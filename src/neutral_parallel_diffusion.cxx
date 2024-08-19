#include "../include/neutral_parallel_diffusion.hxx"

#include <bout/constants.hxx>
#include <bout/fv_ops.hxx>

using bout::globals::mesh;

void NeutralParallelDiffusion::transform(Options& state) {
  AUTO_TRACE();
  Options& allspecies = state["species"];
  for (auto& kv : allspecies.getChildren()) {
    const auto& species_name = kv.first;

    // Get non-const reference
    auto& species = allspecies[species_name];

    if (species.isSet("charge") and (get<BoutReal>(species["charge"]) != 0.0)) {
      // Skip charged species
      continue;
    }

    auto nu = GET_VALUE(Field3D, species["collision_frequency"]);
    const BoutReal AA = GET_VALUE(BoutReal, species["AA"]); // Atomic mass
    const Field3D Nn = GET_VALUE(Field3D, species["density"]);
    const Field3D Tn = GET_VALUE(Field3D, species["temperature"]);
    const Field3D Pn = IS_SET(species["pressure"]) ?
      GET_VALUE(Field3D, species["pressure"]) : Nn * Tn;

    BoutReal advection_factor = 0;
    BoutReal kappa_factor = 0;

    if (equation_fix) {
      advection_factor = (5. / 2);    // This is equivalent to 5/3 if on pressure basis
      kappa_factor = (5. / 2);
    } else {
      advection_factor = (3. / 2);
      kappa_factor = 1;
    }

    // Pressure-diffusion coefficient
    Field3D Dn = dneut * Tn / (AA * nu);
    Dn.applyBoundary("dirichlet_o2");
    mesh->communicate(Dn);

    // Cross-field diffusion calculated from pressure gradient
    // This is the pressure-diffusion approximation 
    Field3D logPn = log(floor(Pn, 1e-7));
    logPn.applyBoundary("neumann");

    // Particle advection
    Field3D S = FV::Div_par_K_Grad_par(Dn * Nn, logPn);
    add(species["density_source"], S);

    Field3D kappa_n = kappa_factor * Nn * Dn;
    kappa_n.applyBoundary("neumann");

    // Heat transfer
    Field3D E = + FV::Div_par_K_Grad_par(
      Dn * advection_factor * Pn, logPn);        // Pressure advection
    if (thermal_conduction) {
      E += FV::Div_par_K_Grad_par(kappa_n, Tn);   // Conduction
    }
    add(species["energy_source"], E);

    Field3D F = 0.0;
    if (IS_SET(species["velocity"]) and viscosity) {
      // Relationship between heat conduction and viscosity for neutral
      // gas Chapman, Cowling "The Mathematical Theory of Non-Uniform
      // Gases", CUP 1952 Ferziger, Kaper "Mathematical Theory of
      // Transport Processes in Gases", 1972
      //

      Field3D Vn = GET_VALUE(Field3D, species["velocity"]);
      Field3D NVn = GET_VALUE(Field3D, species["momentum"]);

      Field3D eta_n = (2. / 5) * kappa_n;

      // Momentum diffusion
      F = FV::Div_par_K_Grad_par(NVn * Dn, logPn) + FV::Div_par_K_Grad_par(eta_n, Vn);
      add(species["momentum_source"], F);
    }

    if (diagnose) {
      // Find the diagnostics struct for this species
      auto search = diagnostics.find(species_name);
      if (search == diagnostics.end()) {
        // First time, create diagnostic
        diagnostics.emplace(species_name, Diagnostics {Dn, S, E, F});
      } else {
        // Update diagnostic values
        auto& d = search->second;
        d.Dn = Dn;
        d.S = S;
        d.E = E;
        d.F = F;
      }
    }
  }
}

void NeutralParallelDiffusion::outputVars(Options &state) {
  AUTO_TRACE();

  if (diagnose) {
    // Normalisations
    auto Nnorm = get<BoutReal>(state["Nnorm"]);
    auto Tnorm = get<BoutReal>(state["Tnorm"]);
    auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
    auto Cs0 = get<BoutReal>(state["Cs0"]);
    auto rho_s0 = get<BoutReal>(state["rho_s0"]);

    for (const auto& it : diagnostics) {
      const std::string& species_name = it.first;
      const auto& d = it.second;
      set_with_attrs(state[std::string("D") + species_name + std::string("_Dpar")], d.Dn,
                   {{"time_dimension", "t"},
                    {"units", "m^2/s"},
                    {"conversion", rho_s0 * Cs0},
                    {"standard_name", "diffusion coefficient"},
                    {"long_name", species_name + " particle diffusion coefficient"},
                    {"species", species_name},
                    {"source", "neutral_parallel_diffusion"}});

      set_with_attrs(state[std::string("S") + species_name + std::string("_Dpar")], d.S,
                   {{"time_dimension", "t"},
                    {"units", "s^-1"},
                    {"conversion", Nnorm * Omega_ci},
                    {"standard_name", "particle diffusion"},
                    {"long_name", species_name + " particle source due to diffusion"},
                    {"species", species_name},
                    {"source", "neutral_parallel_diffusion"}});

      set_with_attrs(state[std::string("E") + species_name + std::string("_Dpar")], d.E,
                   {{"time_dimension", "t"},
                    {"units", "W / m^3"},
                    {"conversion", SI::qe * Tnorm * Nnorm * Omega_ci},
                    {"standard_name", "energy diffusion"},
                    {"long_name", species_name + " energy source due to diffusion"},
                    {"species", species_name},
                    {"source", "neutral_parallel_diffusion"}});

      set_with_attrs(state[std::string("F") + species_name + std::string("_Dpar")], d.F,
                   {{"time_dimension", "t"},
                    {"units", "kg m^-2 s^-2"},
                    {"conversion", SI::Mp * Nnorm * Cs0 * Omega_ci},
                    {"standard_name", "momentum diffusion"},
                    {"long_name", species_name + " momentum source due to diffusion"},
                    {"species", species_name},
                    {"source", "neutral_parallel_diffusion"}});
    }
  }
}
