#pragma once
#ifndef temperature_feedback_H
#define temperature_feedback_H

#include "component.hxx"
#include <bout/constants.hxx>

/// Adds a time-varying temperature source, depending on the difference
/// between the upstream temperature at y=0 and the specified value
struct TemperatureFeedback : public Component {

  /// Inputs
  ///  - <name> (e.g. "d+")
  ///    - temperature_setpoint        Desired temperature in eV
  ///    - control_target_temperature  Adjust the pressure source to match the upstream (if false) or target (if true) temperature
  ///    - temperature_controller_p    Feedback proportional to error
  ///    - temperature_controller_i    Feedback proportional to error integral
  ///    - temperature_integral_positive  Force integral term to be positive? (default: false)
  ///    - temperature_source_positive    Force temperature source to be positive? (default: true)
  ///    - diagnose           Output diagnostic information?
  ///
  ///  - T<name>  (e.g. "Td+")
  ///    - source_shape  The initial source that is scaled by a time-varying factor
  ///
  TemperatureFeedback(std::string name, Options& alloptions, Solver*) : name(name) {

    Options& options = alloptions[name];
    const auto& units = alloptions["units"];
    BoutReal Tnorm = get<BoutReal>(units["eV"]);
    BoutReal Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
    BoutReal Omega_ci = 1. / get<BoutReal>(units["seconds"]);

    BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation
    BoutReal SPnorm = Pnorm * Omega_ci; // Pressure-source normalisation [Pa/s] or [W/m^3] if converted to energy

    temperature_setpoint =
        options["temperature_setpoint"].doc("Desired temperature in eV, upstream (y=0) unless control_target_temperature=true").as<BoutReal>()
        / Tnorm;
    
    control_target_temperature = options["control_target_temperature"]
                                    .doc("Adjust the pressure source to match the upstream (if false) or target (if true) temperature")
                                    .withDefault<bool>(false);

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
    
    species_list = strsplit(options["species_for_temperature_feedback"]
                                  .doc("Comma-separated list of species to apply the PI-controlled source to")
                                  .as<std::string>(),
                            ',');
    scaling_factors_list = strsplit(options["scaling_factors_for_temperature_feedback"]
                                  .doc("Comma-separated list of scaling factors to apply to the PI-controlled source, 1 for each species")
                                  .as<std::string>(),
                            ',');

    if (species_list.size() != scaling_factors_list.size()) {
      throw BoutException("TemperatureFeedback: species_list length doesn't match scaling_factors length. Need 1 scaling factor per species.");
    }

    // NOTE: temperature_error_integral should be set from the restart file here.
    // There doesn't seem to be a way to get the restart Options,
    // so for now this is set in restartVars

    // Source shape the same as used in EvolvePressure
    source_shape =
      (alloptions[std::string("P") + name]["source_shape"]
        .doc("Source term in ddt(P" + name + std::string("). Units [Pa/s], note P = 2/3 E."))
        .as<BoutReal>()
      )  / (SPnorm);

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
      BoutReal SPnorm = Pnorm * Omega_ci; // Pressure-source normalisation [Pa/s] or [W/m^3] if converted to energy

      // Shape is not time-dependent and has units of [Pa s-1]
      set_with_attrs(
          state[std::string("temperature_feedback_src_shape_") + name], source_shape,
          {{"units", "Pa / s"},
           {"conversion", SPnorm},
           {"long_name", name + " temperature source shape"},
           {"source", "temperature_feedback"}});

      // The source multiplier is time-dependent, but dimensionless
      // because all the units are attached to the shape
      set_with_attrs(state[std::string("temperature_feedback_src_mult_") + name],
                     source_multiplier,
                     {{"time_dimension", "t"},
                      {"long_name", name + " temperature source multiplier"},
                      {"source", "temperature_feedback"}});

      set_with_attrs(state[std::string("SP") + name + std::string("_feedback")], source_shape * source_multiplier,
                      {{"time_dimension", "t"},
                      {"units", "Pa / s"},
                    {"conversion", SPnorm},
                    {"standard_name", "temperature source"},
                    {"long_name", name + "upstream temperature feedback controller source"},
                    {"source", "temperature_feedback"}});

      // Save proportional and integral component of source for diagnostics/tuning
      // Multiplier = proportional term + integral term
      set_with_attrs(
          state[std::string("temperature_feedback_src_p_") + name], proportional_term,
          {{"time_dimension", "t"},
          {"long_name", name + " proportional feedback term"},
          {"source", "temperature_feedback"}});

      
      set_with_attrs(
          state[std::string("temperature_feedback_src_i_") + name], integral_term,
          {{"time_dimension", "t"},
           {"long_name", name + " integral feedback term"},
           {"source", "temperature_feedback"}});

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
                    {"source", "temperature_feedback"}});
  }

private:
  std::string name; ///< The species name

  std::list<std::string> species_list; ///< Which species to apply the factor to
  std::list<std::string> scaling_factors_list; ///< Factor to apply

  BoutReal temperature_setpoint;                           ///< Normalised setpoint temperature
  BoutReal temperature_controller_p, temperature_controller_i; ///< PI controller parameters
  BoutReal error;
  BoutReal temperature_error_integral{0.0}; ///< Time integral of the error

  bool control_target_temperature; ///<Adjust the pressure source to match the upstream (if false) or target (if true) temperature
  bool temperature_integral_positive; ///< Force integral term to be positive?
  bool temperature_source_positive;   ///< Force source to be positive?

  // Terms used in Trapezium rule integration of error
  BoutReal temperature_error_lasttime{-1.0};
  BoutReal temperature_error_last{0.0};

  Field3D source_shape; ///< This shape source is scaled up and down

  BoutReal source_multiplier; ///< Factor to multiply source

  BoutReal proportional_term, integral_term; ///< Components of resulting source for diagnostics

  bool diagnose; ///< Output diagnostic information?
};

namespace {
RegisterComponent<TemperatureFeedback> register_uts("temperature_feedback");
}

#endif // temperature_feedback_H
