#pragma once
#ifndef FIXED_DENSITY_H
#define FIXED_DENSITY_H

#include "component.hxx"

/// Set ion density to a fixed value
///
struct FixedDensity : public NamedComponent<FixedDensity> {
  /// Inputs
  /// - <name>
  ///   - AA
  ///   - charge
  ///   - density   value (expression) in units of m^-3
  FixedDensity(std::string name, Options& alloptions, Solver* UNUSED(solver))
      : NamedComponent(name, {readWrite("species:{name}:{vars}")}) {

    auto& options = alloptions[name];

    // Charge and mass
    charge = options["charge"].doc("Particle charge. electrons = -1");
    AA = options["AA"].doc("Particle atomic mass. Proton = 1");

    // Normalisation of density
    const BoutReal Nnorm = alloptions["units"]["inv_meters_cubed"];

    // Get the density and normalise
    N = options["density"].as<Field3D>() / Nnorm;
    substitutePermissions("name", {name});
    substitutePermissions("vars", {"AA", "charge", "density"});
  }

  void outputVars(Options& state) override {
    auto Nnorm = get<BoutReal>(state["Nnorm"]);
    const auto& name = objectName();

    // Save the density, not time dependent
    set_with_attrs(state[std::string("N") + name], N,
                   {{"units", "m^-3"},
                    {"conversion", Nnorm},
                    {"standard_name", "density"},
                    {"long_name", name + " number density"},
                    {"species", name},
                    {"source", "fixed_density"}});
  }

  static constexpr auto type = "fixed_density";

private:
  BoutReal charge; ///< Species charge e.g. electron = -1
  BoutReal AA;     ///< Atomic mass e.g. proton = 1

  Field3D N; ///< Species density (normalised)

  /// Sets in the state the density, mass and charge of the species
  ///
  /// - species
  ///   - <name>
  ///     - AA
  ///     - charge
  ///     - density
  void transform_impl(GuardedOptions& state) override {
    auto species = state["species"][objectName()];
    if (charge != 0.0) { // Don't set charge for neutral species
      set(species["charge"], charge);
    }
    set(species["AA"], AA); // Atomic mass
    set(species["density"], N);
  }
};

namespace {
RegisterComponent<FixedDensity> registercomponentfixeddensity;
}

#endif // FIXED_DENSITY_H
