
#include <difops.hxx>

#include "../include/zero_current.hxx"

void ZeroCurrent::transform(Options &state) {
  // Get the electron pressure
  Options& electrons = state["species"]["e"];
  Field3D Pe = get<Field3D>(electrons["pressure"]);
  Field3D Ne = get<Field3D>(electrons["density"]);

  ASSERT1(get<BoutReal>(electrons["charge"]) == -1.0);

  // Force balance, E = -∇p / n

  Field3D Epar = - Grad_par(Pe) / Ne;

  // Now calculate forces on other species
  Options& allspecies = state["species"];
  for (auto& kv : allspecies.getChildren()) {
    if (kv.first == "e") {
      continue; // Skip electrons
    }
    Options& species = allspecies[kv.first]; // Note: Need non-const
    
    if (!(species.isSet("density") and species.isSet("charge"))) {
      continue; // Needs both density and charge to experience a force
    }
    
    Field3D N = get<Field3D>(species["density"]);
    BoutReal charge = get<BoutReal>(species["charge"]);

    add(species["momentum_source"],
        charge * N * Epar);
  }
}
