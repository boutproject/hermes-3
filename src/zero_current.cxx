
#include <difops.hxx>

#include "../include/zero_current.hxx"

void ZeroCurrent::transform(Options &state) {
  AUTO_TRACE();

  // Get the electron pressure
  // Note: The pressure boundary can be set in sheath boundary condition
  //       which depends on the electron velocity being set here first.
  Options& electrons = state["species"]["e"];
  Field3D Pe = getNoBoundary<Field3D>(electrons["pressure"]);
  Pe.applyBoundary("neumann");

  Field3D Ne = getNoBoundary<Field3D>(electrons["density"]);

  ASSERT1(get<BoutReal>(electrons["charge"]) == -1.0);

  // Force balance, E = (-∇p + F) / n
  Field3D force_density = - Grad_par(Pe);

  if (electrons.isSet("momentum_source")) {
    // Balance other forces from e.g. collisions
    force_density += get<Field3D>(electrons["momentum_source"]);
  }
  const Field3D Epar = force_density / Ne;

  // Current due to ions
  Field3D ion_current;
  
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
    
    const Field3D N = getNoBoundary<Field3D>(species["density"]);
    const BoutReal charge = get<BoutReal>(species["charge"]);

    add(species["momentum_source"],
        charge * N * Epar);

    if (species.isSet("velocity")) {
      // If velocity is set, update the ion current
      
      const Field3D V = getNoBoundary<Field3D>(species["velocity"]);
      
      if (!ion_current.isAllocated()) {
        // Not yet allocated -> Set to the value
        // This avoids having to set to zero initially and add the first time
        ion_current = charge * N * V;
      } else {
        ion_current += charge * N * V;
      }
    }
  }

  if (ion_current.isAllocated()) {
    // If ion current set, set the electron velocity
    set(electrons["velocity"], ion_current / Ne);
  }
}
