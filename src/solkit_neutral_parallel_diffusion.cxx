#include "../include/solkit_neutral_parallel_diffusion.hxx"

#include <bout/fv_ops.hxx>

#include <numeric>

using bout::globals::mesh;

void SOLKITNeutralParallelDiffusion::transform(Options& state) {
  AUTO_TRACE();
  Options& allspecies = state["species"];
  for (auto& kv : allspecies.getChildren()) {
    // Get non-const reference
    auto& species = allspecies[kv.first];

    if (species.isSet("charge") and (get<BoutReal>(species["charge"]) != 0.0)) {
      // Skip charged species
      continue;
    }

    // Sum inverse mean free path with other species
    Field3D inv_mfp = std::accumulate(
      // Iterate through species
      begin(allspecies.getChildren()), end(allspecies.getChildren()),
      // Start with no collisions
      Field3D(0.0),
      [this](Field3D value,
             const std::map<std::string, Options>::value_type &name_species) {
        const Options &species = name_species.second;
        if (name_species.first == "e") {
          // Electrons
          const Field3D Ne = GET_VALUE(Field3D, species["density"]);
          return value + (8.8e-21 / area_norm) * Ne;
        } else if (species.isSet("charge") and
                   (get<BoutReal>(species["charge"]) != 0.0)) {
          // Charged ion species
          const Field3D Ni = GET_VALUE(Field3D, species["density"]);
          return value + (3e-19 / area_norm) * Ni;
        } else {
          // Neutral species, including self
          const Field3D Nn = GET_VALUE(Field3D, species["density"]);
          return value + (8.8e-21 / area_norm) * Nn;
        }
      });
    
    // Diffusion coefficient
    const BoutReal AA = get<BoutReal>(species["AA"]); // Atomic mass
    // Notes:
    //  - neutral_temperature already normalised
    //  - vth,n = sqrt(2kT/m_i)
    Field3D Dn = sqrt(2 * neutral_temperature / AA) / (2 * inv_mfp);
    Dn.applyBoundary("dirichlet_o2");
    mesh->communicate(Dn);

    // Cross-field diffusion calculated from pressure gradient
    const Field3D Nn = GET_VALUE(Field3D, species["density"]);
    const Field3D Pn = Nn * neutral_temperature; // Pressure
    Field3D logPn = log(floor(Pn, 1e-7));
    logPn.applyBoundary("neumann");

    // Particle diffusion
    add(species["density_source"], FV::Div_par_K_Grad_par(Dn * Nn, logPn));
  }
}
