#include "classical_diffusion.hxx"

#include "../include/div_ops.hxx"
#include <bout/fv_ops.hxx>

ClassicalDiffusion::ClassicalDiffusion(std::string name, Options& alloptions, Solver*)
    : NamedComponent(
          name,
          {readIfSet("species:{all_species}:{optional}"), readOnly("species:e:{e_vals}"),
           // FIXME: This is only read velocity or temperature are set for the species.
           readOnly("species:{all_species}:AA"),
           readWrite("species:{all_species}:{output}")}) {
  Options& options = alloptions[name];

  Bsq = SQ(bout::globals::mesh->getCoordinates()->Bxy());

  diagnose = options["diagnose"]
                 .doc("Output the normalised perpendicular particle diffusivity as "
                      "D_classical with conversion metadata for m^2/s, including "
                      "any custom override?")
                 .withDefault<bool>(false);
  custom_D = options["custom_D"]
                 .doc("Particle and thermal-energy diffusivity override [m^2/s]. "
                      ">0: constant diffusion; <0: classical diffusion and conduction "
                      "(default -1); 0: collision-derived particle and carried-energy "
                      "diffusion without conduction.")
                 .withDefault<BoutReal>(-1);
  if (custom_D > 0.0) {
    const Options& units = alloptions["units"];
    const BoutReal rho_s0 = units["meters"];
    const BoutReal Omega_ci = 1.0 / units["seconds"].as<BoutReal>();
    custom_D /= SQ(rho_s0) * Omega_ci;
  }

  low_density_enhance = options["low_density_enhance"]
                            .doc("Smoothly enhance Dn at low electron density in "
                                 "both custom and classical modes?")
                            .withDefault(false);
  if (low_density_enhance) {
    low_density_floor = options["low_density_floor"]
                            .doc("Electron density below which Dn is enhanced [m^-3]. "
                                 "Must be finite and positive.")
                            .withDefault<BoutReal>(1e14);
    low_density_max_multiplier = options["low_density_max_multiplier"]
                                     .doc("Smooth upper bound on the Dn enhancement "
                                          "multiplier. Must be finite and >= 1; "
                                          "1 leaves Dn unchanged.")
                                     .withDefault<BoutReal>(1000.0);
    if (!(low_density_floor > 0.0) || !std::isfinite(low_density_floor)) {
      throw BoutException("{}: low_density_floor must be finite and positive", name);
    }
    if (!(low_density_max_multiplier >= 1.0)
        || !std::isfinite(low_density_max_multiplier)) {
      throw BoutException("{}: low_density_max_multiplier must be finite and >= 1", name);
    }
    low_density_floor /= get<BoutReal>(alloptions["units"]["inv_meters_cubed"]);
  }

  substitutePermissions("optional",
                        {"charge", "pressure", "density", "velocity", "temperature"});
  std::vector<std::string> e_vals = {"AA", "density"};
  if (custom_D <= 0.) {
    e_vals.push_back("collision_frequency");
  }
  substitutePermissions("e_vals", e_vals);
  // Momentum transport requires velocity. Custom thermal-energy diffusion
  // requires pressure or temperature; classical conduction requires temperature
  // and collision frequency. Nothing happens in classical mode if charge or
  // density are unset.
  substitutePermissions("output", {"density_source", "momentum_source", "energy_source"});
  if (custom_D < 0.) {
    setPermissions(readOnly("species:{all_species}:collision_frequency"));
  }
}

void ClassicalDiffusion::transform_impl(GuardedOptions& state) {
  GuardedOptions allspecies = state["species"];

  // Particle diffusion coefficient
  // The only term here comes from the resistive drift

  Field3D Ptotal = 0.0;
  if (!(custom_D > 0)) {
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
  }

  // Particle diffusion coefficient. Applied to all charged species
  // so that net transport is ambipolar

  if (custom_D > 0) { // User-set
    Dn = custom_D;
  } else { // Calculated from collisions
    auto electrons = allspecies["e"];
    const auto me = get<BoutReal>(electrons["AA"]);
    const Field3D Ne = GET_VALUE(Field3D, electrons["density"]);
    const Field3D nu_e =
        floor(GET_VALUE(Field3D, electrons["collision_frequency"]), 1e-10);
    Dn = floor(Ptotal, 1e-5) * me * nu_e / (floor(Ne, 1e-5) * Bsq);
  }

  if (low_density_enhance && low_density_max_multiplier > 1.0) {
    const Field3D Ne = GET_VALUE(Field3D, allspecies["e"]["density"]);
    const Field3D fraction = clamp(Field3D{Ne / low_density_floor}, 0.0, 1.0);
    const Field3D activation = SQ(1.0 - fraction);
    const BoutReal inverse_boost = 1.0 / (low_density_max_multiplier - 1.0);
    // Flatten the onset at the density floor and saturation at zero density.
    // Retain coefficients in guards using the already communicated state.
    Dn *=
        1.0
        + activation
              / (activation * inverse_boost + SQ(fraction) / (fraction + inverse_boost));
  }

  for (auto kv : allspecies.getChildren()) {
    GuardedOptions species = allspecies[kv.first]; // Note: Need non-const
    if (!(custom_D > 0)) {
      if (!(species.isSet("charge") and species.isSet("density"))) {
        continue; // Skip, go to next species
      }
      auto q = get<BoutReal>(species["charge"]);
      if (fabs(q) < 1e-5) {
        continue;
      }
    }
    const auto N = GET_VALUE(Field3D, species["density"]);

    // add(species["density_source"], FV::Div_a_Grad_perp(Dn, N));
    add(species["density_source"],
        Div_a_Grad_perp_nonorthog(Dn, N, cls_pf_perp_xlow, cls_pf_perp_ylow, false));

    if (IS_SET(species["velocity"])) {
      const auto V = GET_VALUE(Field3D, species["velocity"]);
      const auto AA = GET_VALUE(BoutReal, species["AA"]);

      // add(species["momentum_source"], FV::Div_a_Grad_perp(Dn * AA * V, N));
      add(species["momentum_source"],
          Div_a_Grad_perp_nonorthog(Dn * AA * V, N, cls_mf_perp_xlow, cls_mf_perp_ylow,
                                    false));
    }

    if (custom_D > 0) {
      if (IS_SET(species["pressure"]) || IS_SET(species["temperature"])) {
        const Field3D P = IS_SET(species["pressure"])
                              ? GET_VALUE(Field3D, species["pressure"])
                              : N * GET_VALUE(Field3D, species["temperature"]);
        // Diffuse thermal energy and particles with the same coefficient.
        add(species["energy_source"],
            Div_a_Grad_perp_nonorthog(Dn * (3. / 2), P, cls_nef_perp_xlow,
                                      cls_nef_perp_ylow, false));
      }
    } else if (IS_SET(species["temperature"])) {
      const auto T = GET_VALUE(Field3D, species["temperature"]);
      add(species["energy_source"],
          Div_a_Grad_perp_nonorthog(Dn * (3. / 2) * T, N, cls_nef_perp_xlow,
                                    cls_nef_perp_ylow, false));

      if (custom_D < 0) {
        // Cross-field heat conduction
        // kappa_perp = 2 * n * nu_ii * rho_i^2
        const auto P = GET_VALUE(Field3D, species["pressure"]);
        const auto AA = GET_VALUE(BoutReal, species["AA"]);
        const Field3D nu =
            floor(GET_VALUE(Field3D, species["collision_frequency"]), 1e-10);
        add(species["energy_source"],
            Div_a_Grad_perp_nonorthog(2. * floor(P, 1e-5) * nu * AA / Bsq, T,
                                      cls_tef_perp_xlow, cls_tef_perp_ylow, false));
      }
    }
  }
}

void ClassicalDiffusion::outputVars(Options& state) {

  if (diagnose) {
    // Normalisations
    auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
    auto rho_s0 = get<BoutReal>(state["rho_s0"]);

    set_with_attrs(state["D_classical"], Dn,
                   {{"time_dimension", "t"},
                    {"units", "m^2 s^-1"},
                    {"conversion", rho_s0 * rho_s0 * Omega_ci},
                    {"standard_name", "Perpendicular particle diffusion"},
                    {"long_name", "Perpendicular particle diffusion coefficient"},
                    {"source", "classical_diffusion"}});
  }
}
