
#include <bout/difops.hxx>

#include "../include/parallel_ohms_law.hxx"
#include <bout/constants.hxx>
#include "../include/hermes_utils.hxx"

using bout::globals::mesh;

namespace {
/// Limited free gradient of log of a quantity
/// This ensures that the guard cell values remain positive
/// while also ensuring that the quantity never increases
///
///  fm  fc | fp
///         ^ boundary
///
/// exp( 2*log(fc) - log(fm) )
///
BoutReal limitFree(BoutReal fm, BoutReal fc) {
  if (fm < fc) {
    return fc; // Neumann rather than increasing into boundary
  }
  if (fm < 1e-10) {
    return fc; // Low / no density condition
  }
  BoutReal fp = SQ(fc) / fm;
#if CHECKLEVEL >= 2
  if (!std::isfinite(fp)) {
    throw BoutException("SheathBoundary limitFree: {}, {} -> {}", fm, fc, fp);
  }
#endif

  return fp;
}
} // namespace

ParallelOhmsLaw::ParallelOhmsLaw(std::string name, Options& alloptions, Solver*)
    : name(name) {
  AUTO_TRACE();
  auto& options = alloptions[name];
  diagnose = options["diagnose"]
    .doc("Save additional output diagnostics")
    .withDefault<bool>(false);

  resistivity_floor = options["resistivity_floor"].doc("Minimum resistivity floor").withDefault(1e-5);

  spitzer_resistivity = 
        options["spitzer_resistivity"].doc("Use Spitzer resistivity?").withDefault<bool>(false);


  const Options& units = alloptions["units"];
  // Normalisations
  Tnorm = units["eV"];
  Nnorm = units["inv_meters_cubed"];
  rho_s0 = units["meters"];
  Omega_ci = 1. / units["seconds"].as<BoutReal>();
  Cs0 = rho_s0 * Omega_ci; 

  jpar = 0.0; 
  Ve = 0.0 , NVe = 0.0;
  Ve.setBoundary(std::string("Ve"));
  NVe.setBoundary(std::string("NVe"));

}

void ParallelOhmsLaw::calculateResistivity(Options &electrons, Field3D &Ne_lim) {

  // Normalization constant
  const BoutReal eta_norm = Tnorm / (rho_s0 * Nnorm * Cs0 * SI::qe);

  Field3D resistivity_eta;

  if (!spitzer_resistivity){

    // Get electron collision frequency
    const Field3D nu = softFloor(get<Field3D>(electrons["collision_frequency"]), 1e-10); // Nondimentionalized. Multiply with Omega_ci to get nu in [s^-1]

    // Calculate resistivity eta
    resistivity_eta =  ( nu * Omega_ci) * SI::Me / (Ne_lim * Nnorm) / SQ(SI::qe ); // In [Ohm m]

  } else {

    // Get electron temperature
    const Field3D Te = floor(GET_NOBOUNDARY(Field3D, electrons["temperature"]), 0.0);

    BoutReal Zi = 1.0;      // ion charge number
    BoutReal Zeff = 2.5;    // effective charge?

    Field3D LnLambda = 24.0 - log(sqrt(Zi * Ne_lim * Nnorm / 1.e6) / (Te * Tnorm)); // ln(Lambda)

    BoutReal FZ = (1. + 1.198 * Zi + 0.222 * Zi * Zi) / (1. + 2.996 * Zi + 0.753 * Zi * Zi);  // correction coefficient in Spitzer model

    resistivity_eta = Zeff * FZ * 1.03e-4 * Zi * LnLambda * pow(Te * Tnorm,-1.5); // eta in Ohm-m.

  }

  eta = softFloor( resistivity_eta / eta_norm , resistivity_floor );

  mesh->communicate(eta); // Do we need this communication?

}


void ParallelOhmsLaw::transform(Options &state) {
  AUTO_TRACE();

  if (!IS_SET_NOBOUNDARY(state["fields"]["phi"])) {
    // Here we use the electrostatic potential to calculate the parallel current
    throw BoutException("Parallel Ohm's law requires the electrostatic potential. Otherwise use electron force balance\n");
  }

  if (!state["species"].isSection("e")) {
    throw BoutException("Parallel Ohm's law requires a section for electrons in the input file. \n");
  }

  // Get electrostatic potential, without boundary conditions
  Field3D phi = GET_NOBOUNDARY(Field3D, state["fields"]["phi"]);
  
  // Get the electron temperature and pressure, without boundary conditions
  Options& electrons = state["species"]["e"];
  Field3D Te = floor(GET_NOBOUNDARY(Field3D, electrons["temperature"]), 0.0); // Need boundary to take gradient
  Field3D Pe = floor(GET_NOBOUNDARY(Field3D, electrons["pressure"]), 0.0);
  Field3D Ne = floor(GET_NOBOUNDARY(Field3D, electrons["density"]), 0.0);

  // mesh->communicate(Te, Pe, Ne, phi);
  mesh->communicate(Te, Ne);

  // Note: We need boundary conditions on P, so apply the same
  //       free boundary condition as sheath_boundary.
  // Note: The below calculation requires phi derivatives at the Y boundaries
  //       Setting to free boundaries
  Field3D phi_fa = toFieldAligned(phi);
  Field3D Pe_fa = toFieldAligned(Pe);
  Field3D Te_fa = toFieldAligned(Te);
  for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
    for (int jz = 0; jz < mesh->LocalNz; jz++) {
      phi_fa(r.ind, mesh->ystart - 1, jz) = 2 * phi_fa(r.ind, mesh->ystart, jz) - phi_fa(r.ind, mesh->ystart + 1, jz);
      auto i = indexAt(Pe_fa, r.ind, mesh->ystart, jz);
      auto ip = i.yp();
      auto im = i.ym();
      Pe_fa[im] = limitFree(Pe_fa[ip], Pe_fa[i]);
      Te_fa[im] = limitFree(Te_fa[ip], Te_fa[i]);
    }
  }
  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
    for (int jz = 0; jz < mesh->LocalNz; jz++) {
      phi_fa(r.ind, mesh->yend + 1, jz) = 2 * phi_fa(r.ind, mesh->yend, jz) - phi_fa(r.ind, mesh->yend - 1, jz);
      auto i = indexAt(Pe_fa, r.ind, mesh->yend, jz);
      auto ip = i.yp();
      auto im = i.ym();
      Pe_fa[ip] = limitFree(Pe_fa[im], Pe_fa[i]);
      Te_fa[ip] = limitFree(Te_fa[im], Te_fa[i]);
    }
  }
  phi = fromFieldAligned(phi_fa);
  Pe = fromFieldAligned(Pe_fa);
  Te = fromFieldAligned(Te_fa);

  Field3D Ne_lim = softFloor(Ne, 1e-7);

  const BoutReal AA = get<BoutReal>(electrons["AA"]); // Atomic mass
  ASSERT1(get<BoutReal>(electrons["charge"]) == -1.0);

  calculateResistivity(electrons, Ne_lim); // Sets eta

  // Calculate the contribution of each term 
  Field3D term_phi = -Grad_par(phi) / eta;

  Field3D term_Pe = Grad_par(Pe) / Ne_lim / eta;

  Field3D term_Te = 0.71 * Grad_par(Te) / eta;

  jpar = term_phi + term_Pe + term_Te;

  // Current due to other species
  Field3D current = 0.0;

  // Now calculate current from all other species
  Options& allspecies = state["species"];
  for (auto& kv : allspecies.getChildren()) {
    if (kv.first == "e") {
      continue; // Skip electrons
    }

    Options& species = allspecies[kv.first]; // Note: Need non-const
    if (!(species.isSet("density") and species.isSet("charge"))) {
      continue; // Needs both density and charge to contribute
    }

    if (isSetFinalNoBoundary(species["velocity"], "parallel_ohms_law")) {
      // If velocity is set, update the current
      // Note: Mark final so can't be set later

      const Field3D N = getNoBoundary<Field3D>(species["density"]);
      const BoutReal charge = get<BoutReal>(species["charge"]);
      const Field3D V = getNoBoundary<Field3D>(species["velocity"]);

      if (fabs(charge) < 1e-5) {
        continue; // Not charged
      }

      current += charge * N * V;
    }
  }

  Ve = (jpar - current) / (-1.0 * Ne_lim);
  Ve.applyBoundary();

  NVe = AA * Ne * Ve; // Re-calculate consistent with V and N
  NVe.applyBoundary();

  mesh->communicate(Ve, NVe);

  set(electrons["velocity"], Ve);
  set(electrons["momentum"], NVe);
}

void ParallelOhmsLaw::outputVars(Options& state) {
  AUTO_TRACE();

    set_with_attrs(state["Ve"], Ve,
                   {{"time_dimension", "t"},
                    {"units", "m / s"},
                    {"conversion", Cs0},
                    {"long_name", "e parallel velocity"},
                    {"standard_name", "velocity"},
                    {"species", "e"},
                    {"source", "parallel_ohms_law"}});


    set_with_attrs(state["NVe"], NVe,
                   {{"time_dimension", "t"},
                    {"units", "kg / m^2 / s"},
                    {"conversion", SI::Mp * Nnorm * Cs0},
                    {"standard_name", "momentum"},
                    {"long_name", "e parallel momentum"},
                    {"species", "e"},
                    {"source", "parallel_ohms_law"}});

  if (diagnose) {
    set_with_attrs(state["jpar"], jpar,
                  {{"time_dimension", "t"},
                    {"units", "A / m^2"},
                    {"conversion", SI::qe * Nnorm * Cs0},
                    {"standard_name", "jpar"},
                    {"long_name", "Parallel electric current"},
                    {"source", "parallel_ohms_law"}});


    const BoutReal eta_norm = Tnorm / (rho_s0 * Nnorm * Cs0 * SI::qe);

    set_with_attrs(state["eta_resist"], eta,
                  {{"time_dimension", "t"},
                    {"units", "Ohm m"},
                    {"conversion", eta_norm},
                    {"standard_name", "eta"},
                    {"long_name", "plasma resistivity"},
                    {"source", "parallel_ohms_law"}});
  }
}
