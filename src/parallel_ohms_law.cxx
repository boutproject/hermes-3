
#include <bout/difops.hxx>

#include "../include/parallel_ohms_law.hxx"
#include <bout/constants.hxx>
#include "../include/hermes_utils.hxx"

using bout::globals::mesh;


namespace {
BoutReal clip(BoutReal value, BoutReal min, BoutReal max) {
  if (value < min)
    return min;
  if (value > max)
    return max;
  return value;
}

Ind3D indexAt(const Field3D& f, int x, int y, int z) {
  int ny = f.getNy();
  int nz = f.getNz();
  return Ind3D{(x * ny + y) * nz + z, ny, nz};
}

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

}

ParallelOhmsLaw::ParallelOhmsLaw(std::string name, Options& alloptions, Solver*)
    : name(name) {
  AUTO_TRACE();
  auto& options = alloptions[name];
  diagnose = options["diagnose"]
    .doc("Save additional output diagnostics")
    .withDefault<bool>(false);

  resistivity_floor = options["resistivity_floor"].doc("Minimum resistivity floor").withDefault(1e-4);

  spitzer_resist = 
        options["spitzer_resist"].doc("Use Spitzer resistivity?").withDefault<bool>(false);


  const Options& units = alloptions["units"];
  // Normalisations
  Tnorm = units["eV"];
  Nnorm = units["inv_meters_cubed"];
  rho_s0 = units["meters"];
  Omega_ci = 1. / units["seconds"].as<BoutReal>();
  Cs0 = rho_s0 * Omega_ci; 


  Ve = 0.0 , NVe = 0.0;
  Ve.setBoundary(std::string("Ve"));
  NVe.setBoundary(std::string("NVe"));

  jpar = 0.0; 


  lower_y = options["lower_y"].doc("Boundary on lower y?").withDefault<bool>(true);
  upper_y = options["upper_y"].doc("Boundary on upper y?").withDefault<bool>(true);

  // Read wall voltage, convert to normalised units
  wall_potential = options["wall_potential"]
                       .doc("Voltage of the wall [Volts]")
                       .withDefault(Field3D(0.0))
                   / Tnorm;
  // Convert to field aligned coordinates
  wall_potential = toFieldAligned(wall_potential);

  // Note: wall potential at the last cell before the boundary is used,
  // not the value at the boundary half-way between cells. This is due
  // to how twist-shift boundary conditions and non-aligned inputs are
  // treated; using the cell boundary gives incorrect results.

  floor_potential = options["floor_potential"]
                        .doc("Apply a floor to wall potential when calculating Ve?")
                        .withDefault<bool>(true);

  temperature_floor = options["temperature_floor"].doc("Low temperature scale")
    .withDefault<BoutReal>(0.1) / Tnorm;

  Ge = options["secondary_electron_coef"]
           .doc("Effective secondary electron emission coefficient")
           .withDefault(0.0);

  if ((Ge < 0.0) or (Ge > 1.0)) {
    throw BoutException("Secondary electron emission must be between 0 and 1 ({:e})", Ge);
  }
  
}



void ParallelOhmsLaw::calculateResistivity(Options &electrons, Field3D &Ne_lim) {

  // Normalization constant
  const BoutReal eta_norm = Tnorm / (rho_s0 * Nnorm * Cs0 * SI::qe);

  Field3D resistivity_eta;

  if (!spitzer_resist){

    // Get electron collision frequency
    const Field3D nu = softFloor(get<Field3D>(electrons["collision_frequency"]), 1e-10); // Nondimentionalized. Multiply with Omega_ci to get nu in [s^-1]

    // Calculate resistivity eta
    resistivity_eta =  ( nu * Omega_ci) * SI::Me / (Ne_lim * Nnorm) / SQ(SI::qe ); // In [Ohm m]

  } else {

    // Get electron temperature
    const Field3D Te = floor(GET_VALUE(Field3D, electrons["temperature"]), 0.0);

    BoutReal Zi = 1.0;      // ion charge number
    BoutReal Zeff = 2.5;    // effective charge?

    Field3D LnLambda = 24.0 - log(sqrt(Zi * Ne_lim * Nnorm / 1.e6) / (Te * Tnorm)); // ln(Lambda)

    BoutReal FZ = (1. + 1.198 * Zi + 0.222 * Zi * Zi) / (1. + 2.996 * Zi + 0.753 * Zi * Zi);  // correction coefficient in Spitzer model

    resistivity_eta = Zeff * FZ * 1.03e-4 * Zi * LnLambda * pow(Te * Tnorm,-1.5); // eta in Ohm-m.

  }

  eta = softFloor( resistivity_eta / eta_norm , resistivity_floor );
  // eta = sqrt( SQ(resistivity_eta / eta_norm) + resistivity_floor*resistivity_floor ); // This mitigates discontinuities.

  mesh->communicate(eta); // NOTE(malamast): Do we need this communication?

}


void ParallelOhmsLaw::transform(Options &state) {
  AUTO_TRACE();

  if (!IS_SET(state["fields"]["phi"])) {
    // Here we use the electrostatic potential to calculate the parallel current
    throw BoutException("Parallel Ohm's law requires the electrostatic potential. Otherwise use electron force balance\n");
  }

  if (!state["species"].isSection("e")) {
    throw BoutException("Parallel Ohm's law requires a section for electrons in the input file. \n");
  }

  // Get electrostatic potential, with boundary condition applied
  Field3D phi = get<Field3D>(state["fields"]["phi"]);

  // Get the electron temperature and pressure, with boundary condition applied
  Options& electrons = state["species"]["e"];
  Field3D Te = floor(GET_VALUE(Field3D, electrons["temperature"]), 0.0); // Need boundary to take gradient
  Field3D Pe = floor(GET_VALUE(Field3D, electrons["pressure"]), 0.0);
  Field3D Ne = floor(GET_VALUE(Field3D, electrons["density"]), 0.0);

  // mesh->communicate(Te, Pe, Ne, phi);
  mesh->communicate(Te, Ne); // NOTE(malamast): do we need that?

  Field3D Ne_lim = softFloor(Ne, 1e-7);

  const BoutReal AA = get<BoutReal>(electrons["AA"]); // Atomic mass
  ASSERT1(get<BoutReal>(electrons["charge"]) == -1.0);

  calculateResistivity(electrons, Ne_lim); // Sets eta

  // Calculate the contribution of each term 
  Field3D term_phi = -Grad_par(phi) / eta;

  Field3D term_Pe = Grad_par(Pe) / Ne_lim / eta;

  Field3D term_Te = 0.71 * Grad_par(Te) / eta;


  // NOTE(malamast): Should we include the contribution of the electron momentum source terms? How can we exclude the contribution of collisions there?
  // if (IS_SET(electrons["momentum_source"])) {
  //   // Balance other forces from e.g. collisions
  //   // Note: marked as final so can't be changed later
  //   force_density += GET_VALUE(Field3D, electrons["momentum_source"]);
  // }

  jpar = term_phi + term_Pe + term_Te;

  // Current due to other species
  // Field3D current;
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

      // if (!current.isAllocated()) {
        // Not yet allocated -> Set to the value
        // This avoids having to set to zero initially and add the first time
        // current = charge * N * V;
      // } else {
      current += charge * N * V;
      // }
    }
  }

  Ve = (jpar - current) / (-1.0 * Ne_lim);
  Ve.applyBoundary();

  NVe = AA * Ne * Ve; // Re-calculate consistent with V and N
  NVe.applyBoundary();

  mesh->communicate(Ve, NVe);


  // set(electrons["velocity"], Ve); // # This is causing a problem with option.hasAttribute("final-domain") check in set function.
  electrons["velocity"].force(Ve);

  // set(electrons["momentum"], NVe);
  electrons["momentum"].force(NVe);


  //////////////////////////////////////////////////////////////////
  // sheath BC for Electrons

  // Mass, normalised to proton mass
  const BoutReal Me =
      IS_SET(electrons["AA"]) ? get<BoutReal>(electrons["AA"]) : SI::Me / SI::Mp;

  // Need electron properties
  // Not const because boundary conditions will be set
  Field3D Ne_fa = toFieldAligned(Ne);
  Field3D Te_fa = toFieldAligned(Te);
  Field3D Pe_fa = toFieldAligned(Pe);

  Field3D phi_fa  = toFieldAligned(phi);

  // This is for applying boundary conditions
  Field3D Ve_fa = toFieldAligned(Ve);
  Field3D NVe_fa = toFieldAligned(NVe);


  if (lower_y) {
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        auto i = indexAt(Ne_fa, r.ind, mesh->ystart, jz);
        auto ip = i.yp();
        auto im = i.ym();


        const BoutReal nesheath = 0.5 * (Ne_fa[im] + Ne_fa[i]);
        const BoutReal tesheath = 0.5 * (Te_fa[im] + Te_fa[i]);  // electron temperature
        const BoutReal phi_wall = wall_potential[i];

        const BoutReal phisheath = floor_potential ? floor(
            0.5 * (phi_fa[im] + phi_fa[i]), phi_wall) // Electron saturation at phi = phi_wall
	        : 0.5 * (phi_fa[im] + phi_fa[i]);


        // Electron velocity into sheath (< 0)
        const BoutReal vesheath = (tesheath < 1e-10) ?
          0.0 :
          -sqrt(tesheath / (TWOPI * Me)) * (1. - Ge) * exp(-(phisheath - phi_wall) / tesheath);

        Ve_fa[im] = 2 * vesheath - Ve_fa[i];
        NVe_fa[im] = 2. * Me * nesheath * vesheath - NVe_fa[i];

      }
    }
  }

  if (upper_y) {
    // This is essentially the same as at the lower y boundary
    // except ystart -> yend, ip <-> im
    // 
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        auto i = indexAt(Ne, r.ind, mesh->yend, jz);
        auto ip = i.yp();
        auto im = i.ym();

        const BoutReal nesheath = 0.5 * (Ne_fa[ip] + Ne_fa[i]);
        const BoutReal tesheath = 0.5 * (Te_fa[ip] + Te_fa[i]);  // electron temperature
        const BoutReal phi_wall = wall_potential[i];
        const BoutReal phisheath = floor_potential ? floor(
            0.5 * (phi_fa[ip] + phi_fa[i]), phi_wall) // Electron saturation at phi = phi_wall
          : 0.5 * (phi_fa[ip] + phi_fa[i]);

        // Electron velocity into sheath (> 0)
        const BoutReal vesheath = (tesheath < 1e-10) ?
          0.0 :
          sqrt(tesheath / (TWOPI * Me)) * (1. - Ge) * exp(-(phisheath - phi_wall) / tesheath);

        Ve_fa[ip] = 2 * vesheath - Ve_fa[i];
        NVe_fa[ip] = 2. * Me * nesheath * vesheath - NVe_fa[i];

      }
    }
  }

  Ve_fa.clearParallelSlices();
  setBoundary(electrons["velocity"], fromFieldAligned(Ve_fa));

  NVe_fa.clearParallelSlices();
  setBoundary(electrons["momentum"], fromFieldAligned(NVe_fa));

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
