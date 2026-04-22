#include "classical_diffusion.hxx"

#include "../include/div_ops.hxx"
#include <bout/fv_ops.hxx>

ClassicalDiffusion::ClassicalDiffusion(std::string name, Options& alloptions, Solver*)
    : Component({readIfSet("species:{all_species}:{optional}"),
                 readOnly("species:{all_species}:AA"), 
                 readOnly("species:e:{e_vals}"),
                 readWrite("species:{all_species}:{output}")}), name(name) {
  Options& options = alloptions[name];

  Bsq = SQ(bout::globals::mesh->getCoordinates()->Bxy);
  zero_BC_transport =
      options["zero_BC_transport"].doc("Zero BC transport zero?").withDefault<bool>(true);
  diagnose = options["diagnose"].doc("Output additional diagnostics?").withDefault<bool>(false);
  custom_D = options["custom_D"].doc("Custom diffusion coefficient override. -1: Off, calculate D normally").withDefault<BoutReal>(-1);
  nonorthogonal_operators =
      options["nonorthogonal_operators"]
          .doc("Use nonorthogonal operators for radial transport? NOTE: may be broken")
          .withDefault<bool>(false);

  substitutePermissions("optional",
                        {"charge", "pressure", "density", "velocity", "temperature"});
  std::vector<std::string> e_vals = {"AA", "density"};
  if (custom_D <= 0.) {
    e_vals.push_back("collision_frequency");
  }
  substitutePermissions("e_vals", e_vals);
  // FIXME: momentum and energy sources are only set if velocity and
  // temperature are defined (respectively). Collision frequency is
  // only used if temperature is set. Nothing happens if the charge or
  // density are unset.
  std::vector<std::string> output_vars;
  output_vars.push_back("density_source");
  output_vars.push_back("cls_particle_flow_xlow");
  output_vars.push_back("cls_particle_flow_ylow");
  output_vars.push_back("energy_source");
  output_vars.push_back("cls_energy_flow_xlow");
  output_vars.push_back("cls_energy_flow_ylow");
  output_vars.push_back("momentum_source");
  output_vars.push_back("cls_momentum_flow_xlow");
  output_vars.push_back("cls_momentum_flow_ylow");
  substitutePermissions("output", output_vars);
  setPermissions(readWrite("species:{all_species}:energy_flow_xlow"));
  setPermissions(readWrite("species:{all_species}:momentum_flow_xlow"));
  setPermissions(readWrite("species:{all_species}:particle_flow_xlow"));
  setPermissions(readWrite("species:{all_species}:energy_flow_ylow"));
  setPermissions(readWrite("species:{all_species}:momentum_flow_ylow"));
  setPermissions(readWrite("species:{all_species}:particle_flow_ylow"));
  if (custom_D < 0.)
    setPermissions(readOnly("species:{all_species}:collision_frequency"));
}

void ClassicalDiffusion::transform_impl(GuardedOptions& state) {
  GuardedOptions allspecies = state["species"];
  GuardedOptions electrons = allspecies["e"];
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

  if (custom_D > 0) {    // User-set
    Dn = custom_D;   
  } else {                  // Calculated from collisions
    auto electrons = allspecies["e"];
    const auto me = get<BoutReal>(electrons["AA"]);
    const Field3D Ne = GET_VALUE(Field3D, electrons["density"]);
    const Field3D nu_e = floor(GET_VALUE(Field3D, electrons["collision_frequency"]), 1e-10);
    Dn = floor(Ptotal, 1e-5) * me * nu_e / (floor(Ne, 1e-5) * Bsq);
  }

  // Set D to zero in all guard cells
  if (zero_BC_transport) {
    BOUT_FOR(i, Dn.getRegion("RGN_GUARDS")) {
      Dn[i] = 0.0;
    }
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

    if (nonorthogonal_operators)
    {
    add(species["density_source"],
        Div_a_Grad_perp_nonorthog(Dn, N, cls_pf_perp_xlow, cls_pf_perp_ylow));
    }
    else
    {
      Div_a_Grad_perp_flows(Dn, N, cls_pf_perp_xlow, cls_pf_perp_ylow);
    } 
      add(species["cls_particle_flow_xlow"], cls_pf_perp_xlow);
      add(species["cls_particle_flow_ylow"], cls_pf_perp_ylow);
      add(species["particle_flow_xlow"], cls_pf_perp_xlow);
      add(species["particle_flow_ylow"], cls_pf_perp_ylow);

      if (IS_SET(species["velocity"])) {
        const auto V = GET_VALUE(Field3D, species["velocity"]);
        const auto AA = GET_VALUE(BoutReal, species["AA"]);

        // add(species["momentum_source"], FV::Div_a_Grad_perp(Dn * AA * V, N));
        if (nonorthogonal_operators) {
          add(species["momentum_source"],
              Div_a_Grad_perp_nonorthog(Dn * AA * V, N, cls_mf_perp_xlow,
                                        cls_mf_perp_ylow));
        } else {
          add(species["momentum_source"],
              Div_a_Grad_perp_flows(Dn * AA * V, N, cls_mf_perp_xlow, cls_mf_perp_ylow));
        }
        add(species["cls_momentum_flow_xlow"], cls_mf_perp_xlow);
        add(species["cls_momentum_flow_ylow"], cls_mf_perp_ylow);
        add(species["momentum_flow_xlow"], cls_mf_perp_xlow);
        add(species["momentum_flow_ylow"], cls_mf_perp_ylow);
    }
    if (IS_SET(species["temperature"])) {
      const auto T = GET_VALUE(Field3D, species["temperature"]);
      // add(species["energy_source"], FV::Div_a_Grad_perp(Dn * (3. / 2) * T, N));
      if (nonorthogonal_operators)
      {
      add(species["energy_source"],
          Div_a_Grad_perp_nonorthog(Dn * (3. / 2) * T, N, cls_nef_perp_xlow,
                                    cls_nef_perp_ylow));
      }
      else{
      add(species["energy_source"],
          Div_a_Grad_perp_flows(Dn * (3. / 2) * T, N, cls_nef_perp_xlow,
                                    cls_nef_perp_ylow));
      }
    add(species["cls_energy_flow_xlow"], cls_nef_perp_xlow);
    add(species["cls_energy_flow_ylow"], cls_nef_perp_ylow);
    add(species["energy_flow_xlow"], cls_nef_perp_xlow);
    add(species["energy_flow_ylow"], cls_nef_perp_ylow);

    // TODO: Figure out what to do with the below
    if (custom_D < 0) {
      // Cross-field heat conduction
      // kappa_perp = 2 * n * nu_ii * rho_i^2
      const auto P = GET_VALUE(Field3D, species["pressure"]);
      const auto AA = GET_VALUE(BoutReal, species["AA"]);
      nu = floor(GET_VALUE(Field3D, species["collision_frequency"]), 1e-10);
      // add(species["energy_source"], FV::Div_a_Grad_perp(2. * floor(P, 1e-5) * nu * AA
      // / Bsq, T));
      Kappa_perp = 2. * floor(P, 1e-5) * nu * AA / Bsq;
      if (zero_BC_transport) {
        BOUT_FOR(i, Kappa_perp.getRegion("RGN_GUARDS")) { Kappa_perp[i] = 0.0; }
      }
      if (nonorthogonal_operators) {
        add(species["energy_source"],
            Div_a_Grad_perp_nonorthog(Kappa_perp, T, cls_tef_perp_xlow,
                                      cls_tef_perp_ylow));
      } else {
        add(species["energy_source"],
            Div_a_Grad_perp_flows(Kappa_perp, T, cls_tef_perp_xlow, cls_tef_perp_ylow));
      }
      add(species["cls_energy_flow_xlow"], cls_tef_perp_xlow);
      add(species["cls_energy_flow_ylow"], cls_tef_perp_ylow);
      add(species["energy_flow_xlow"], cls_tef_perp_xlow);
      add(species["energy_flow_ylow"], cls_tef_perp_ylow);
    }
    }
  }
}

void ClassicalDiffusion::outputVars(Options& state) {

  if (diagnose) {
    // Normalisations
    auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
    auto rho_s0 = get<BoutReal>(state["rho_s0"]);
    auto Nnorm = get<BoutReal>(state["Nnorm"]);

    set_with_attrs(state["D_classical"], Dn,
                   {{"time_dimension", "t"},
                    {"units", "m^2 s^-1"},
                    {"conversion", rho_s0 * rho_s0 * Omega_ci},
                    {"standard_name", "Classical particle diffusion"},
                    {"long_name", "Classical cross-field particle diffusion coefficient"},
                    {"source", "classical_diffusion"}});
    set_with_attrs(state["Kappa_perp_classical"], Kappa_perp,
                   {{"time_dimension", "t"},
                    {"units", "m^2 s^-1"},
                    {"conversion", rho_s0 * rho_s0 * Omega_ci}, //TODO: WRONG
                    {"standard_name", "Classical perpendicular heat conductivity"},
                    {"long_name", "Classical cross-field perpendicular heat conductivity"},
                    {"source", "classical_diffusion"}});
    set_with_attrs(state["nu_classical"], nu,
                   {{"time_dimension", "t"},
                    {"units", "s^-1"},
                    {"conversion", Omega_ci},
                    {"standard_name", "Classical collision frequency"},
                    {"long_name", "Classical cross-field collision frequency"},
                    {"source", "classical_diffusion"} });
    set_with_attrs(state[fmt::format("pf{}_classical_xlow", name)], cls_pf_perp_xlow,
                   {{"time_dimension", "t"},
                    {"units", "s^-1"},
                    {"conversion", rho_s0 * SQ(rho_s0) * Nnorm * Omega_ci},
                    {"standard_name", "classical particle flow"},
                    {"long_name", name + " classical transport particle flow in X. "},
                    {"species", name},
                    {"source", "classical_diffusion"}});
    }

}
