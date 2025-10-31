#pragma once
#ifndef FIXED_VELOCITY_H
#define FIXED_VELOCITY_H

#include "component.hxx"
#include <bout/globals.hxx>

/// Set parallel velocity to a fixed value
///
struct FixedVelocity : public Component {

  FixedVelocity(std::string name, Options& alloptions, Solver* UNUSED(solver))
      : name(name) {
    AUTO_TRACE();

    auto& options = alloptions[name];

    // Normalisation of velocity
    auto& units = alloptions["units"];
    const BoutReal Cs0 = units["meters"].as<BoutReal>() / units["seconds"].as<BoutReal>();

    // Get the velocity and normalise
    // First read from the mesh file e.g. "Ve0"
    if ((bout::globals::mesh->get(V, std::string("V") + name + "0") != 0) and
        !options.isSet("velocity")) {
      throw BoutException("fixed_velocity: Missing mesh V{}0 or option {}:velocity\n", name, name);
    }
    // Option overrides mesh value
    // so use mesh value (if any) as default value.
    V = options["velocity"].withDefault(V) / Cs0;
  }

  void outputVars(Options& state) override {
    AUTO_TRACE();
    auto Cs0 = get<BoutReal>(state["Cs0"]);

    // Save the density, not time dependent
    set_with_attrs(state[std::string("V") + name], V,
                   {{"units", "m / s"},
                    {"conversion", Cs0},
                    {"long_name", name + " parallel velocity"},
                    {"standard_name", "velocity"},
                    {"species", name},
                    {"source", "fixed_velocity"}});
  }

private:
  std::string name; ///< Short name of species e.g "e"

  Field3D V; ///< Species velocity (normalised)

  /// This sets in the state
  /// - species
  ///   - <name>
  ///     - velocity
  ///     - momentum
  void transform(GuardedOptions& state) override {
    AUTO_TRACE();
    auto& species = state["species"][name];
    set(species["velocity"], V);

    // If density is set, also set momentum
    if (isSetFinalNoBoundary(species["density"])) {
      const Field3D N = getNoBoundary<Field3D>(species["density"]);
      const BoutReal AA = get<BoutReal>(species["AA"]); // Atomic mass

      set(species["momentum"], AA * N * V);
    }
  }
};

namespace {
RegisterComponent<FixedVelocity> registercomponentfixedvelocity("fixed_velocity");
}

#endif // FIXED_VELOCITY_H
