#pragma once
#ifndef FIXED_VELOCITY_H
#define FIXED_VELOCITY_H

#include "component.hxx"
#include "guarded_options.hxx"
#include "permissions.hxx"

#include <bout/assert.hxx>
#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/field3d.hxx>
#include <bout/globals.hxx>
#include <bout/mesh.hxx>
#include <bout/options.hxx>
#include <bout/unused.hxx>

#include <string>

/// Set parallel velocity to a fixed value
///
struct FixedVelocity : public NamedComponent<FixedVelocity> {

  FixedVelocity(std::string name, Options& alloptions, Solver* UNUSED(solver))
      : NamedComponent(name, {readIfSet("species:{name}:density", Regions::Interior),
                              // FIXME: AA is only read if density is set
                              readOnly("species:{name}:AA"),
                              readWrite("species:{name}:{output}")}),
        name(name) {

    auto& options = alloptions[name];

    // Normalisation of velocity
    auto& units = alloptions["units"];
    const BoutReal Cs0 = units["meters"].as<BoutReal>() / units["seconds"].as<BoutReal>();

    // Get the velocity and normalise
    // First read from the mesh file e.g. "Ve0"
    if ((bout::globals::mesh->get(V, std::string("V") + name + "0") != 0)
        and !options.isSet("velocity")) {
      throw BoutException("fixed_velocity: Missing mesh V{}0 or option {}:velocity\n",
                          name, name);
    }
    // Option overrides mesh value
    // so use mesh value (if any) as default value.
    V = options["velocity"].withDefault(V) / Cs0;

    bout::globals::mesh->communicate(V);
    if (V.isFci()) {
      V.applyParallelBoundary("parallel_neumann_o2");
      ASSERT2(V.hasParallelSlices());
    }

    substitutePermissions("name", {name});
    // FIXME: Momentum is only written if density is set
    substitutePermissions("output", {"velocity", "momentum"});
  }

  void outputVars(Options& state) override {
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

  static constexpr auto type = "fixed_velocity";

private:
  std::string name; ///< Short name of species e.g "e"

  Field3D V; ///< Species velocity (normalised)

  /// This sets in the state
  /// - species
  ///   - <name>
  ///     - velocity
  ///     - momentum
  void transform_impl(GuardedOptions& state) override {
    auto species = state["species"][name];
    set(species["velocity"], V);

    // If density is set, also set momentum
    if (isSetFinalNoBoundary(species["density"])) {
      const Field3DParallel N = getNoBoundary<Field3D>(species["density"]);
      const BoutReal AA = get<BoutReal>(species["AA"]); // Atomic mass

      set(species["momentum"], Field3DParallel{AA * N * V});
    }
  }
};

namespace {
RegisterComponent<FixedVelocity> registercomponentfixedvelocity;
}

#endif // FIXED_VELOCITY_H
