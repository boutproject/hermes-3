#pragma once
#ifndef FIXED_TEMPERATURE_H
#define FIXED_TEMPERATURE_H

#include "component.hxx"

/// Set species temperature to a fixed value
///
struct FixedTemperature : public Component {
  /// Inputs
  /// - <name>
  ///   - temperature   value (expression) in units of eV
  FixedTemperature(std::string name, Options& alloptions, Solver* UNUSED(solver))
      : name(name) {
    AUTO_TRACE();

    auto& options = alloptions[name];

    // Normalisation of temperature
    auto Tnorm = get<BoutReal>(alloptions["units"]["eV"]);

    // Get the temperature and normalise
    T = options["temperature"].doc("Constant temperature [eV]").as<Field3D>()
        / Tnorm; // Normalise

    diagnose = options["diagnose"]
      .doc("Save additional output diagnostics")
      .withDefault<bool>(false);
  }

  /// Sets in the state the temperature and pressure of the species
  ///
  /// Inputs
  /// - species
  ///   - <name>
  ///     - density (optional)
  ///
  /// Sets in the state
  ///
  /// - species
  ///   - <name>
  ///     - temperature
  ///     - pressure (if density is set)
  void transform(Options& state) override {
    AUTO_TRACE();
    auto& species = state["species"][name];

    set(species["temperature"], T);

    // If density is set, also set pressure
    if (isSetFinalNoBoundary(species["density"])) {
      // Note: The boundary of N may not be set yet
      auto N = GET_NOBOUNDARY(Field3D, species["density"]);
      P = N * T;
      set(species["pressure"], P);
    }
  }

  void outputVars(Options& state) override {
    AUTO_TRACE();
    auto Tnorm = get<BoutReal>(state["Tnorm"]);

    // Save temperature to output files
    set_with_attrs(state[std::string("T") + name], T,
                   {{"units", "eV"},
                    {"conversion", Tnorm},
                    {"standard_name", "temperature"},
                    {"long_name", name + " temperature"},
                    {"species", name},
                    {"source", "fixed_temperature"}});

    if (diagnose) {
      auto Nnorm = get<BoutReal>(state["Nnorm"]);
      BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

      // Save pressure as time-varying field
      set_with_attrs(state[std::string("P") + name], P,
                     {{"time_dimension", "t"},
                      {"units", "Pa"},
                      {"conversion", Pnorm},
                      {"standard_name", "pressure"},
                      {"long_name", name + " pressure"},
                      {"species", name},
                      {"source", "fixed_temperature"}});
    }
  }

private:
  std::string name; ///< Short name of species e.g "e"

  Field3D T; ///< Species temperature (normalised)
  Field3D P; ///< Species pressure (normalised)

  bool diagnose; ///< Output additional fields
};

namespace {
RegisterComponent<FixedTemperature> registercomponentfixedtemperature("fixed_temperature");
}

#endif // FIXED_TEMPERATURE_H
