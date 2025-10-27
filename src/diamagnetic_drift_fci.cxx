#include <bout/fv_ops.hxx>
#include <bout/vecops.hxx>
#include <bout/yboundary_regions.hxx>
#include "../include/diamagnetic_drift_fci.hxx"
using bout::globals::mesh;

DiamagneticDriftFCI::DiamagneticDriftFCI(std::string name, Options& alloptions,
                                   Solver* UNUSED(solver)) {

  // Get options for this component
  auto& options = alloptions[name];

  yboundary.init(options);

  // Normalise

  // Get the units
  const auto& units = alloptions["units"];
  BoutReal Bnorm = get<BoutReal>(units["Tesla"]);
  BoutReal Lnorm = get<BoutReal>(units["meters"]);

  const auto coord = mesh->getCoordinates();
  if (mesh->isFci()) {

    bracket_factor = sqrt(coord->g_22.withoutParallelSlices()) / (coord->J.withoutParallelSlices() * coord->Bxy);

    logB = log(coord->Bxy);
    logB.applyBoundary("neumann_o2");
    mesh->communicate(logB);
    logB.applyParallelBoundary("parallel_neumann_o2");

  } else {
    bracket_factor = 1.0;
    logB = log(coord->Bxy);
    mesh->communicate(logB);
    logB.applyBoundary("neumann_o2");
  } 
}

void DiamagneticDriftFCI::transform(Options& state) {
  // Iterate through all subsections
  Options& allspecies = state["species"];

  for (auto& kv : allspecies.getChildren()) {
    Options& species = allspecies[kv.first]; // Note: Need non-const

    if (!(species.isSet("charge") and species.isSet("temperature")))
      continue; // Skip, go to next species

    // Calculate diamagnetic drift velocity for this species
    auto q = get<BoutReal>(species["charge"]);
    if (fabs(q) < 1e-5) {
      continue;
    }
    auto T = GET_VALUE(Field3D, species["temperature"]);

    if (IS_SET(species["density"])) {
      auto N = GET_VALUE(Field3D, species["density"]);

      Field3D jdia_bracket = 2 * bracket(logB, T*N/q, BRACKET_ARAKAWA) * bracket_factor;

      subtract(species["density_source"], jdia_bracket);
    }

    if (IS_SET(species["pressure"])) {
      auto P = get<Field3D>(species["pressure"]);

      Field3D jdia_bracket = 2 * bracket(logB, T*P/q, BRACKET_ARAKAWA) * bracket_factor;
											  
      subtract(species["energy_source"], (5. / 2) * jdia_bracket);
    }

    if (IS_SET(species["momentum"])) {
      auto NV = get<Field3D>(species["momentum"]);
      
      Field3D jdia_bracket = 2 * bracket(logB, T*NV/q, BRACKET_ARAKAWA) * bracket_factor;

      subtract(species["momentum_source"], jdia_bracket);
    }
  }
}
