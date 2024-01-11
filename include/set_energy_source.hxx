#pragma once
#ifndef SET_ENERGY_SOURCE_H
#define SET_ENERGY_SOURCE_H

#include "component.hxx"
#include <bout/constants.hxx>

/// Set species energy source to the energy source of another species
///
/// # Example
///
/// [hermes]
/// components = e, d, ...
///
/// [e]
/// type = ... // energy source should be set
///
/// [d+]
/// type = set_energy_source, ...
///
/// energy_source_from = e // Set Pd+ = Pe
struct SetEnergySource : public Component {
  /// Inputs
  /// - <name>
  ///   - energy_source_from   name of species
  SetEnergySource(std::string name, Options& alloptions, Solver* UNUSED(solver))
      : name(name) {
    AUTO_TRACE();

    auto& options = alloptions[name];

    energy_source_from = options["energy_source_from"]
                           .doc("Name of species to take energy source from (e.g 'e')")
                           .as<std::string>();
    
    energy_source_ratio = options["energy_source_ratio"]
                               .doc("Scalar factor to apply to the energy source.")
                               .withDefault(1.0);

    diagnose = options["diagnose"]
                   .doc("Save additional output diagnostics")
                   .withDefault<bool>(false);
  }

  ///
  /// Inputs
  /// - species
  ///   - <energy_source_from>
  ///     - energy_source
  ///
  /// Sets in the state:
  /// - species
  ///   - <name>
  ///     - energy_source
  ///
  void transform(Options& state) override {
    AUTO_TRACE();

    // Get the energy_source
    Field3D E_src = getNoBoundary<Field3D>(state["species"][energy_source_from]["energy_source"]);

    // Set energy_source
    auto& species = state["species"][name];
    add(species["energy_source"], energy_source_ratio * E_src);

  }

  void outputVars(Options& state) override {
    AUTO_TRACE();

    if (diagnose) {
      auto Tnorm = get<BoutReal>(state["Tnorm"]);
      auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
      auto Nnorm = get<BoutReal>(state["Nnorm"]);
      BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation
      BoutReal SPnorm = Pnorm * Omega_ci; // Pressure-source normalisation [Pa/s] or [W/m^3] if converted to energy

    set_with_attrs(state[std::string("SP") + name], E_src,
                      {{"time_dimension", "t"},
                      {"units", "Pa / s"},
                      {"conversion", SPnorm},
                      {"standard_name", "energy source"},
                      {"long_name", name + "energy source"},
                      {"source", "set_energy_source"}});
    }
  }

private:
  std::string name;               ///< Short name of species e.g "e"
  std::string energy_source_from; ///< The species that the energy_source is taken from
  BoutReal energy_source_ratio;   ///< Constant factor applied to the energy_source
  Field3D E_src;                  ///< The energy_source
  bool diagnose;                  ///< Output diagnostics?
};

namespace {
RegisterComponent<SetEnergySource> register_set_energy_source("set_energy_source");
}

#endif // SET_ENERGY_SOURCE_H
