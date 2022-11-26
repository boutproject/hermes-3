#pragma once
#ifndef FIXED_VELOCITY_H
#define FIXED_VELOCITY_H

#include "component.hxx"

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
    V = options["velocity"].as<Field3D>() / Cs0;
  }

  /// This sets in the state
  /// - species
  ///   - <name>
  ///     - velocity
  ///     - momentum
  void transform(Options& state) override {
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

  void outputVars(Options& state) override {
    AUTO_TRACE();
    auto Cs0 = state["Cs0"].as<BoutReal>();

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
};

namespace {
RegisterComponent<FixedVelocity> registercomponentfixedvelocity("fixed_velocity");
}

#endif // FIXED_VELOCITY_H
