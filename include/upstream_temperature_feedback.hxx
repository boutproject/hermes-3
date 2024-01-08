#pragma once
#ifndef UPSTREAM_TEMPERATURE_FEEDBACK_H
#define UPSTREAM_TEMPERATURE_FEEDBACK_H

#include "component.hxx"

/// Adds a time-varying temperature source, depending on the difference
/// between the upstream temperature at y=0 and the specified value
struct UpstreamTemperatureFeedback : public Component {

  /// Inputs
  ///  - <name> (e.g. "d+")
  ///    - temperature_upstream        Upstream temperature (y=0) in eV
  ///    - temperature_controller_p    Feedback proportional to error
  ///    - temperature_controller_i    Feedback proportional to error integral
  ///    - temperature_integral_positive  Force integral term to be positive? (default: false)
  ///    - temperature_source_positive    Force temperature source to be positive? (default: true)
  ///    - diagnose           Output diagnostic information?
  ///
  ///  - N<name>  (e.g. "Nd+")
  ///    - source_shape  The initial source that is scaled by a time-varying factor
  ///
  UpstreamTemperatureFeedback(std::string name, Options& alloptions, Solver*) : name(name) {
    const auto& units = alloptions["units"];
    BoutReal Tnorm = get<BoutReal>(units["eV"]);
    BoutReal FreqNorm = 1. / get<BoutReal>(units["seconds"]);

    Options& options = alloptions[name];

    temperature_upstream =
        options["temperature_upstream"].doc("Upstream temperature (at y=0) [eV]").as<BoutReal>()
        / Tnorm;

    temperature_controller_p = options["temperature_controller_p"]
                               .doc("Feedback controller proportional (p) parameter")
                               .withDefault(1e-2);
    temperature_controller_i = options["temperature_controller_i"]
                               .doc("Feedback controller integral (i) parameter")
                               .withDefault(1e-3);

    temperature_integral_positive = options["temperature_integral_positive"]
                                    .doc("Force integral term to be positive?")
                                    .withDefault<bool>(false);
    temperature_source_positive = options["temperature_source_positive"]
                                  .doc("Force source to be positive?")
                                  .withDefault<bool>(true);

    // NOTE: temperature_error_integral should be set from the restart file here.
    // There doesn't seem to be a way to get the restart Options,
    // so for now this is set in restartVars

    // Source shape the same as used in EvolvePressure
    pressure_source_shape =
      alloptions[std::string("P") + name]["source_shape"]
            .doc("Source term in ddt(P" + name + std::string("). Units [W/s]"))
            .withDefault(Field3D(0.0))
        / (Tnorm * FreqNorm);

    diagnose = options["diagnose"]
                   .doc("Output additional diagnostics?")
                   .withDefault<bool>(false);
  }

  /// Inputs
  ///  - <name>
  ///    - temperature
  ///
  /// Outputs
  ///
  ///  - <name>
  ///    - temperature_source
  ///
  void transform(Options& state) override;

  void outputVars(Options& state) override {
    AUTO_TRACE();
    if (diagnose) {
      // Normalisations
      auto Tnorm = get<BoutReal>(state["Tnorm"]);
      auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
      auto Nnorm = get<BoutReal>(state["Nnorm"]);
      BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

      // Shape is not time-dependent and has units of [W s-1]
      set_with_attrs(
          state[std::string("temperature_feedback_src_shape_") + name], temperature_source_shape,
          {{"units", "W s^-1"},
           {"conversion", Pnorm * Omega_ci},
           {"long_name", name + " temperature source shape"},
           {"source", "upstream_temperature_feedback"}});

      // The source multiplier is time-dependent, but dimensionless
      // because all the units are attached to the shape
      set_with_attrs(state[std::string("temperature_feedback_src_mult_") + name],
                     source_multiplier,
                     {{"time_dimension", "t"},
                      {"long_name", name + " temperature source multiplier"},
                      {"source", "upstream_temperature_feedback"}});

      set_with_attrs(state[std::string("S") + name + std::string("_feedback")], temperature_source_shape * source_multiplier,
                      {{"time_dimension", "t"},
                      {"units", "W s^-1"},
                    {"conversion", Pnorm * Omega_ci},
                    {"standard_name", "temperature source"},
                    {"long_name", name + "upstream temperature feedback controller source"},
                    {"source", "upstream_temperature_feedback"}});

      // Save proportional and integral component of source for diagnostics/tuning
      // Multiplier = proportional term + integral term
      set_with_attrs(
          state[std::string("temperature_feedback_src_p_") + name], proportional_term,
          {{"time_dimension", "t"},
          {"long_name", name + " proportional feedback term"},
          {"source", "upstream_temperature_feedback"}});

      
      set_with_attrs(
          state[std::string("temperature_feedback_src_i_") + name], integral_term,
          {{"time_dimension", "t"},
           {"long_name", name + " integral feedback term"},
           {"source", "upstream_temperature_feedback"}});

    }
  }

  void restartVars(Options& state) override {
    AUTO_TRACE();

    // NOTE: This is a hack because we know that the loaded restart file
    //       is passed into restartVars in PhysicsModel::postInit
    // The restart value should be used in init() rather than here
    static bool first = true;
    if (first and state.isSet(name + "_temperature_error_integral")) {
      first = false;
      temperature_error_integral = state[name + "_temperature_error_integral"].as<BoutReal>();
    }

    // Save the temperature error integral
    set_with_attrs(state[name + "_temperature_error_integral"], temperature_error_integral,
                   {{"long_name", name + " temperature error integral"},
                    {"source", "upstream_temperature_feedback"}});
  }

private:
  std::string name; ///< The species name

  BoutReal temperature_upstream;                           ///< Normalised upstream temperature
  BoutReal temperature_controller_p, temperature_controller_i; ///< PI controller parameters
  BoutReal error;
  BoutReal temperature_error_integral{0.0}; ///< Time integral of the error

  bool temperature_integral_positive; ///< Force integral term to be positive?
  bool temperature_source_positive;   ///< Force source to be positive?

  // Terms used in Trapezium rule integration of error
  BoutReal temperature_error_lasttime{-1.0};
  BoutReal temperature_error_last{0.0};

  Field3D temperature_source_shape; ///< This shape source is scaled up and down

  BoutReal source_multiplier; ///< Factor to multiply source

  BoutReal proportional_term, integral_term; ///< Components of resulting source for diagnostics

  bool diagnose; ///< Output diagnostic information?
};

namespace {
RegisterComponent<UpstreamTemperatureFeedback> register_uds("upstream_temperature_feedback");
}

#endif // UPSTREAM_TEMPERATURE_FEEDBACK_H
