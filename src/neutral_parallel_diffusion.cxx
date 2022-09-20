#include "../include/neutral_parallel_diffusion.hxx"

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

    // Diffusion coefficient
    Field3D Dn = dneut * Tn / (AA * nu);
    Dn.applyBoundary("dirichlet_o2");
    mesh->communicate(Dn);

    // Cross-field diffusion calculated from pressure gradient
    Field3D logPn = log(floor(Pn, 1e-7));
    logPn.applyBoundary("neumann");

    // Particle diffusion
    Field3D S = FV::Div_par_K_Grad_par(Dn * Nn, logPn);
    add(species["density_source"], S);

    Field3D kappa_n = Nn * Dn;
    kappa_n.applyBoundary("neumann");

    // Heat conduction
    Field3D E = FV::Div_par_K_Grad_par(kappa_n, Tn) // Parallel
      + FV::Div_par_K_Grad_par(Dn * (3. / 2) * Pn, logPn); // Perpendicular diffusion
    add(species["energy_source"], E);

    Field3D F = 0.0;
    if (IS_SET(species["velocity"])) {
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
        auto it_bool_pair = diagnostics.emplace(species_name, Diagnostics {Dn, S, E, F});
        auto& d = it_bool_pair.first->second;
        bout::globals::dump.addRepeat(d.Dn, std::string("D") + species_name + std::string("_Dpar"));
        bout::globals::dump.addRepeat(d.S, std::string("S") + species_name + std::string("_Dpar"));
        bout::globals::dump.addRepeat(d.E, std::string("E") + species_name + std::string("_Dpar"));
        bout::globals::dump.addRepeat(d.F, std::string("F") + species_name + std::string("_Dpar"));
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
