
#include <bout/difops.hxx>

#include "../include/parallel_ohms_law.hxx"
#include <bout/constants.hxx>

using bout::globals::mesh;


ParallelOhmsLaw::ParallelOhmsLaw(std::string name, Options& alloptions, Solver*)
    : name(name) {
  AUTO_TRACE();
  auto& options = alloptions[name];
  diagnose = options["diagnose"]
    .doc("Save additional output diagnostics")
    .withDefault<bool>(false);

  resistivity_floor = options["resistivity_floor"].doc("Minimum resistivity floor").withDefault(1e-6);


  Ve.setBoundary(std::string("Ve"));
  NVe.setBoundary(std::string("NVe"));

  const Options& units = alloptions["units"];
  // Normalisations
  Tnorm = units["eV"];
  Nnorm = units["inv_meters_cubed"];
  rho_s0 = units["meters"];
  Omega_ci = 1. / units["seconds"].as<BoutReal>();
  Cs0 = rho_s0 * Omega_ci; 
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
  const Field3D Te = GET_VALUE(Field3D, electrons["temperature"]); // Need boundary to take gradient
  const Field3D Pe = GET_VALUE(Field3D, electrons["pressure"]);

  auto Ne = getNoBoundary<Field3D>(electrons["density"]);
  const Field3D Ne_lim = floor(Ne, 1e-7);

  const BoutReal AA = get<BoutReal>(electrons["AA"]); // Atomic mass
  ASSERT1(get<BoutReal>(electrons["charge"]) == -1.0);


  // Get electron collision frequency
  const Field3D nu = GET_VALUE(Field3D, electrons["collision_frequency"]) * Omega_ci; // In s^-1

  // Calculate resistivity eta
  Field3D resistivity_eta =  nu * SI::Me / (Ne_lim * Nnorm) / SQ(SI::qe ); // In [Ohm m]
  const BoutReal eta_norm = Tnorm / (rho_s0 * Nnorm * Cs0 * SI::qe);

  Field3D eta = floor(resistivity_eta / eta_norm,  resistivity_floor);



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
  // set(electrons["velocity"], Ve); # This is causing a problem with option.hasAttribute("final-domain") check in set function.
  electrons["velocity"].force(Ve);

  NVe = AA * Ne * Ve; // Re-calculate consistent with V and N
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

  }
}
