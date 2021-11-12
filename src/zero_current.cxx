
#include <difops.hxx>

#include "../include/zero_current.hxx"

ZeroCurrent::ZeroCurrent(std::string name, Options& alloptions, Solver*)
    : name(name) {
  AUTO_TRACE();
  Options &options = alloptions[name];

  charge = options["charge"].doc("Particle charge. electrons = -1");

  // Save the velocity
  bout::globals::dump.addRepeat(velocity, std::string("V") + name);
  velocity = 0.0;

  ASSERT0(charge != 0.0);
}

void ZeroCurrent::transform(Options &state) {
  AUTO_TRACE();

  // Current due to other species
  Field3D current;
  
  // Now calculate forces on other species
  Options& allspecies = state["species"];
  for (auto& kv : allspecies.getChildren()) {
    if (kv.first == name) {
      continue; // Skip self
    }
    Options& species = allspecies[kv.first]; // Note: Need non-const
 
    if (!(species.isSet("density") and species.isSet("charge"))) {
      continue; // Needs both density and charge to contribute
    }

    if (isSetFinalNoBoundary(species["velocity"], "zero_current")) {
      // If velocity is set, update the current
      // Note: Mark final so can't be set later

      const Field3D N = getNoBoundary<Field3D>(species["density"]);
      const BoutReal charge = get<BoutReal>(species["charge"]);
      const Field3D V = getNoBoundary<Field3D>(species["velocity"]);

      if (!current.isAllocated()) {
        // Not yet allocated -> Set to the value
        // This avoids having to set to zero initially and add the first time
        current = charge * N * V;
      } else {
        current += charge * N * V;
      }
    }
  }

  if (!current.isAllocated()) {
    // No currents, probably not what is intended
    throw BoutException("No other species to set velocity from");
  }

  // Get the species density
  Options& species = state["species"][name];
  Field3D N = getNoBoundary<Field3D>(species["density"]);

  velocity = current / (-charge * N);
  set(species["velocity"], velocity);
}
