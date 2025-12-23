#pragma once
#ifndef FIXED_DENSITY_H
#define FIXED_DENSITY_H

#include "component.hxx"

using bout::globals::mesh;

/// Set ion density to a fixed value
///
struct FixedDensity : public Component {
  /// Inputs
  /// - <name>
  ///   - AA
  ///   - charge
  ///   - density   value (expression) in units of m^-3
  FixedDensity(std::string name, Options& alloptions, Solver* UNUSED(solver))
      : name(name) {
    AUTO_TRACE();

    auto& options = alloptions[name];

    // Charge and mass
    charge = options["charge"].doc("Particle charge. electrons = -1");
    AA = options["AA"].doc("Particle atomic mass. Proton = 1");

    // Normalisation of density
    const BoutReal Nnorm = alloptions["units"]["inv_meters_cubed"];

    // Get the density and normalise
    N = options["density"].as<Field3D>() / Nnorm;

    // Initilize the Field3D N from 2D profiles stored in the mesh file.
    initialize_from_mesh = options["initialize_from_mesh"]
      .doc("Initilize field from 2D profiles stored in the mesh file?")
      .withDefault<bool>(false);

    if (initialize_from_mesh) {
      // Try to read the initial field from the mesh
      mesh->get(N, std::string("N") + name + "_init"); // Units: [m^-3/s]
      N /= Nnorm; // Normalization
    }    

  }

  /// Sets in the state the density, mass and charge of the species
  ///
  /// - species
  ///   - <name>
  ///     - AA
  ///     - charge
  ///     - density
  void transform(Options& state) override {
    AUTO_TRACE();
    auto& species = state["species"][name];
    if (charge != 0.0) { // Don't set charge for neutral species
      set(species["charge"], charge);
    }
    set(species["AA"], AA); // Atomic mass
    set(species["density"], N);
  }

  void outputVars(Options& state) override {
    AUTO_TRACE();
    auto Nnorm = get<BoutReal>(state["Nnorm"]);

    // Save the density, not time dependent
    set_with_attrs(state[std::string("N") + name], N,
                   {{"units", "m^-3"},
                    {"conversion", Nnorm},
                    {"standard_name", "density"},
                    {"long_name", name + " number density"},
                    {"species", name},
                    {"source", "fixed_density"}});
  }
private:
  std::string name; ///< Short name of species e.g "e"

  BoutReal charge; ///< Species charge e.g. electron = -1
  BoutReal AA;     ///< Atomic mass e.g. proton = 1

  Field3D N; ///< Species density (normalised)

  bool initialize_from_mesh;  ///< Initilize the Field3D N from 2D profiles stored in the mesh file. 

};

namespace {
RegisterComponent<FixedDensity> registercomponentfixeddensity("fixed_density");
}

#endif // FIXED_DENSITY_H
