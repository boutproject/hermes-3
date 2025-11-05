#include "classical_diffusion.hxx"

#include <bout/fv_ops.hxx>

ClassicalDiffusion::ClassicalDiffusion(std::string name, Options& alloptions, Solver*)
    : Component({readIfSet("species:{all_species}:{optional}"),
                 readOnly("species:e:{e_vals}"),
                 readWrite("species:{all_species}:{output}")}) {
  AUTO_TRACE();
  Options& options = alloptions[name];

  Bsq = SQ(bout::globals::mesh->getCoordinates()->Bxy);

  diagnose = options["diagnose"].doc("Output additional diagnostics?").withDefault<bool>(false);
  custom_D = options["custom_D"].doc("Custom diffusion coefficient override. -1: Off, calculate D normally").withDefault<BoutReal>(-1);

  state_variable_access.substitute(
      "optional", {"charge", "pressure", "density", "velocity", "temperature"});
  std::vector<std::string> e_vals = {"AA", "density"};
  if (custom_D <= 0.)
    e_vals.push_back("collision_frequency");
  state_variable_access.substitute("e_vals", e_vals);
  // FIXME: momentum and energy sources are only set if velocity and
  // temperature are defined (respectively). Collision frequency is
  // only used if temperature is set. Nothing happens if the charge or
  // density are unset.
  state_variable_access.substitute(
      "output", {"density_source", "momentum_source", "energy_source"});
  if (custom_D < 0.)
    state_variable_access.setAccess(
        readOnly("species:{all_species}:collision_frequency"));
}

void ClassicalDiffusion::transform_impl(GuardedOptions& state) {
  AUTO_TRACE();
  GuardedOptions allspecies = state["species"];
  
  // Particle diffusion coefficient
  // The only term here comes from the resistive drift

  Field3D Ptotal = 0.0;
  for (auto& kv : allspecies.getChildren()) {
    const auto species = kv.second;

    if (!(species.isSet("charge") and species.isSet("pressure"))) {
      continue; // Skip, go to next species
    }
    auto q = get<BoutReal>(species["charge"]);
    if (fabs(q) < 1e-5) {
      continue;
    }
    Ptotal += GET_VALUE(Field3D, species["pressure"]);
  }

  auto electrons = allspecies["e"];
  const auto me = get<BoutReal>(electrons["AA"]);
  const Field3D Ne = GET_VALUE(Field3D, electrons["density"]);

  // Particle diffusion coefficient. Applied to all charged species
  // so that net transport is ambipolar

  if (custom_D > 0) {    // User-set
    Dn = custom_D;   
  } else {                  // Calculated from collisions
    const Field3D nu_e = floor(GET_VALUE(Field3D, electrons["collision_frequency"]), 1e-10);
    Dn = floor(Ptotal, 1e-5) * me * nu_e / (floor(Ne, 1e-5) * Bsq);
  }

  // Set D to zero in all guard cells
  BOUT_FOR(i, Dn.getRegion("RGN_GUARDS")) {
      Dn[i] = 0.0;
    }

  for (auto kv : allspecies.getChildren()) {
    GuardedOptions species = allspecies[kv.first]; // Note: Need non-const

    if (!(species.isSet("charge") and species.isSet("density"))) {
      continue; // Skip, go to next species
    }
    auto q = get<BoutReal>(species["charge"]);
    if (fabs(q) < 1e-5) {
      continue;
    }

    const auto N = GET_VALUE(Field3D, species["density"]);

    add(species["density_source"], FV::Div_a_Grad_perp(Dn, N));

    if (IS_SET(species["velocity"])) {
      const auto V = GET_VALUE(Field3D, species["velocity"]);
      const auto AA = GET_VALUE(BoutReal, species["AA"]);

      add(species["momentum_source"], FV::Div_a_Grad_perp(Dn * AA * V, N));
    }

    if (IS_SET(species["temperature"])) {
      const auto T = GET_VALUE(Field3D, species["temperature"]);
      add(species["energy_source"], FV::Div_a_Grad_perp(Dn * (3. / 2) * T, N));

      // Cross-field heat conduction
      // kappa_perp = 2 * n * nu_ii * rho_i^2

      const auto P = GET_VALUE(Field3D, species["pressure"]);
      const auto AA = GET_VALUE(BoutReal, species["AA"]);

      // TODO: Figure out what to do with the below
      if(custom_D < 0) {
        const Field3D nu = floor(GET_VALUE(Field3D, species["collision_frequency"]), 1e-10);
        add(species["energy_source"], FV::Div_a_Grad_perp(2. * floor(P, 1e-5) * nu * AA / Bsq, T));
      }
    }
  }
}

void ClassicalDiffusion::outputVars(Options &state) {
  AUTO_TRACE();

  if (diagnose) {
    // Normalisations
    auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
    auto rho_s0 = get<BoutReal>(state["rho_s0"]);

    set_with_attrs(state["D_classical"], Dn,
                   {{"time_dimension", "t"},
                    {"units", "m^2 s^-1"},
                    {"conversion", rho_s0 * rho_s0 * Omega_ci},
                    {"standard_name", "Classical particle diffusion"},
                    {"long_name", "Classical cross-field particle diffusion coefficient"},
                    {"source", "classical_diffusion"}});
  }
}
