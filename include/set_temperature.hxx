#pragma once
#ifndef SET_TEMPERATURE_H
#define SET_TEMPERATURE_H

#include "component.hxx"

/// Set species temperature to the temperature of another species
///
/// # Example
///
/// [hermes]
/// components = e, d, ...
///
/// [e]
/// type = ... // Evolve Te
///
/// [d]
/// type = set_temperature, ...
///
/// temperature_from = e // Set Td = Te
struct SetTemperature : public Component {
  /// Inputs
  /// - <name>
  ///   - temperature_from   name of species
  SetTemperature(std::string name, Options& alloptions, Solver* UNUSED(solver))
      : name(name) {
    AUTO_TRACE();

    auto& options = alloptions[name];

    temperature_from = options["temperature_from"]
      .doc("Name of species to take temperature from (e.g 'e')")
      .as<std::string>();

    if (options["diagnose"]
            .doc("Save additional output diagnostics")
            .withDefault<bool>(false)) {
      // Save temperature as time-varying field
      bout::globals::dump.addRepeat(T, std::string("T") + name);
      T = 0.0;
    }
  }

  ///
  /// Inputs
  /// - species
  ///   - <temperature_from>
  ///     - temperature
  ///
  /// Sets in the state:
  /// - species
  ///   - <name>
  ///     - temperature
  ///     - pressure (if density is set)
  ///
  void transform(Options& state) override {
    AUTO_TRACE();

    // Get the temperature
    T = GET_NOBOUNDARY(Field3D, state["species"][temperature_from]["temperature"]);

    // Set temperature
    auto& species = state["species"][name];
    set(species["temperature"], T);

    if (isSetFinalNoBoundary(species["density"])) {
      // Note: The boundary of N may not be set yet
      auto N = GET_NOBOUNDARY(Field3D, species["density"]);
      set(species["pressure"], N * T);
    }
  }

private:
  std::string name;             ///< Short name of species e.g "e"
  std::string temperature_from; //< The species that the temperature is taken from
  Field3D T;                    //< The temperature
};

namespace {
RegisterComponent<SetTemperature> registercomponentsettemperature("set_temperature");
}

#endif // SET_TEMPERATURE_H
