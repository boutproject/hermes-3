
#include <bout/difops.hxx>
#include <bout/immersed_boundary.hxx>

#include "../include/zero_current.hxx"

ZeroCurrent::ZeroCurrent(std::string name, Options& alloptions, Solver*)
    : name(name) {
  AUTO_TRACE();
  Options &options = alloptions[name];

  charge = options["charge"].doc("Particle charge. electrons = -1");

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

      //IB_TODO: I think * operator is fine here. velocity set with current and the BCs set after so ghosts dont matter here.
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
  if (species["velocity"].isSet()) {
    throw BoutException("Cannot use zero_current in species {} if velocity already set\n", name);
  }
  Field3D N = getNoBoundary<Field3D>(species["density"]);

  Field3D Nfloor = copy(N);
  if (immBndry) {
    //IB_TODO: Just floor, dont need to SetBoundary because velocity will below and loop ignores ghost cells here.
    immBndry->FloorField(Nfloor, 1e-5);
    velocity = current;
    BOUT_FOR(i, velocity.getRegion("RGN_NO_IMM_BNDRY")) {
      velocity[i] /= (-charge * Nfloor[i]);
    }
    velocity.name = std::string("V") + name;
    //IB_TODO: Complex logic for late setup...
    if (!immBndry->CheckFieldSetUp(velocity.name)) {
      immBndry->FieldSetup(velocity);
    }
    immBndry->SetBoundary(velocity);
  } else {
    Field3D Nfloor = floor(N, 1e-5);
    velocity = current / (-charge * Nfloor);
  }
  set(species["velocity"], velocity);
}

void ZeroCurrent::outputVars(Options &state) {
  AUTO_TRACE();
  auto Cs0 = get<BoutReal>(state["Cs0"]);

  // Save the velocity
  set_with_attrs(state[std::string("V") + name], velocity,
                 {{"time_dimension", "t"},
                  {"units", "m / s"},
                  {"conversion", Cs0},
                  {"long_name", name + " parallel velocity"},
                  {"standard_name", "velocity"},
                  {"species", name},
                  {"source", "zero_current"}});
}
