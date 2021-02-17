#include "../include/neutral_parallel_diffusion.hxx"

#include <bout/fv_ops.hxx>

using bout::globals::mesh;

void NeutralParallelDiffusion::transform(Options& state) {

  Options& allspecies = state["species"];
  for (auto& kv : allspecies.getChildren()) {
    // Get non-const reference
    auto& species = allspecies[kv.first];

    if (species.isSet("charge") and (get<BoutReal>(species["charge"]) != 0.0)) {
      // Skip charged species
      continue;
    }

    auto nu = get<Field3D>(species["collision_frequency"]);
    const BoutReal AA = get<BoutReal>(species["AA"]); // Atomic mass
    Field3D Nn = get<Field3D>(species["density"]);
    Field3D Tn = get<Field3D>(species["temperature"]);
    Field3D Pn = get<Field3D>(species["pressure"]);

    // Diffusion coefficient
    Field3D Dn = dneut * Tn / (AA * nu);
    Dn.applyBoundary("dirichlet_o2");
    mesh->communicate(Dn);

    // Cross-field diffusion calculated from pressure gradient
    Field3D logPn = log(floor(Pn, 1e-7));
    logPn.applyBoundary("neumann");

    // Particle diffusion
    add(species["density_source"], FV::Div_par_K_Grad_par(Dn * Nn, logPn));

    Field3D kappa_n = Nn * Dn;
    kappa_n.applyBoundary("neumann");

    // Heat conduction
    add(species["energy_source"],
        FV::Div_par_K_Grad_par(kappa_n, Tn) // Parallel
            + FV::Div_par_K_Grad_par(Dn * (3. / 2) * Pn,
                                     logPn)); // Perpendicular diffusion

    if (species.isSet("velocity")) {
      // Relationship between heat conduction and viscosity for neutral
      // gas Chapman, Cowling "The Mathematical Theory of Non-Uniform
      // Gases", CUP 1952 Ferziger, Kaper "Mathematical Theory of
      // Transport Processes in Gases", 1972
      //

      Field3D Vn = get<Field3D>(species["velocity"]);
      Field3D NVn = get<Field3D>(species["momentum"]);

      Field3D eta_n = (2. / 5) * kappa_n;

      // Momentum diffusion
      add(species["momentum_source"],
          FV::Div_par_K_Grad_par(NVn * Dn, logPn) + FV::Div_par_K_Grad_par(eta_n, Vn));
    }
  }
}
