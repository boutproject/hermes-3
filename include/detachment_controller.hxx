#pragma once
#ifndef detachment_controller_H
#define detachment_controller_H

#include "component.hxx"
#include <bout/constants.hxx>

struct DetachmentController : public Component {

  DetachmentController(std::string, Options& options, Solver*) {

    Options& detachment_controller_options = options["detachment_controller"];

    const auto& units = options["units"];
    BoutReal Tnorm = get<BoutReal>(units["eV"]);
    BoutReal Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
    BoutReal Omega_ci = 1. / get<BoutReal>(units["seconds"]);

    BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

    connection_length = options["mesh"]["length"].as<BoutReal>();

    detachment_front_desired_location = 
      detachment_controller_options["detachment_front_desired_location"]
      .doc("Desired location of the detachment front relative to the divertor target in [m].")
      .as<BoutReal>();
    
    // Convert from distance to target into the y coordinate.
    detachment_front_desired_location = connection_length - detachment_front_desired_location;

    error_on_log_distance = 
      detachment_controller_options["error_on_log_distance"]
      .doc("If true, calculate the difference in log-space instead of linear space. Equivalent to taking log10(detachment_front_desired_location/detachment_front_location).")
      .withDefault<bool>(false);
    
    minval_for_log = 
      detachment_controller_options["minval_for_log"]
      .doc("A floor value used instead of 0 when evaluating error_on_log_distance.")
      .withDefault<BoutReal>(1E-12);

    exponential_control = 
      detachment_controller_options["exponential_control"]
      .doc("Raise the PI correction to the power of 10.")
      .withDefault<bool>(false);
    
    species_for_source_shape =
      detachment_controller_options["species_for_source_shape"]
      .doc("Which species to select the source_shape from?")
      .withDefault<std::string>("e");
    
    neutral_species =
      detachment_controller_options["neutral_species"]
      .doc("Which is the main neutral species?")
      .as<std::string>();
    
    control_power = detachment_controller_options["control_power"]
                    .doc("Control power source (if true) or particle source (if false)?")
                    .withDefault<bool>(false);
    
    controller_p = detachment_controller_options["detachment_controller_p"]
                               .doc("Feedback controller proportional (p) parameter")
                               .withDefault(1e-2);
    controller_i = detachment_controller_options["detachment_controller_i"]
                               .doc("Feedback controller integral (i) parameter")
                               .withDefault(1e-3);

    force_integral_positive = detachment_controller_options["force_integral_positive"]
                                    .doc("Force integral term to be positive?")
                                    .withDefault<bool>(false);
    force_source_positive = detachment_controller_options["force_source_positive"]
                                  .doc("Force source to be positive?")
                                  .withDefault<bool>(true);
    
    species_list = strsplit(detachment_controller_options["species_for_detachment_feedback"]
                                  .doc("Comma-separated list of species to apply the PI-controlled source to")
                                  .as<std::string>(),
                            ',');
    scaling_factors_list = strsplit(detachment_controller_options["scaling_factors_for_detachment_feedback"]
                                  .doc("Comma-separated list of scaling factors to apply to the PI-controlled source, 1 for each species")
                                  .as<std::string>(),
                            ',');

    if (species_list.size() != scaling_factors_list.size()) {
      throw BoutException("DetachmentController: species_list length doesn't match scaling_factors length. Need 1 scaling factor per species.");
    }

    // NOTE: detachment_control_error_integral should be set from the restart file here.
    // There doesn't seem to be a way to get the restart Options,
    // so for now this is set in restartVars

    // Source shape the same as used in EvolvePressure
    if (control_power) {
      source_shape =
        (options[std::string("P") + species_for_source_shape]["source_shape"]
          .doc("Source term in ddt(P" + species_for_source_shape + std::string("). Units [Pa/s], note P = 2/3 E."))
          .withDefault(Field3D(0.0))
        )  / (Pnorm * Omega_ci);

        source_units = "Pa / s";
        source_conversion = Pnorm * Omega_ci;

    } else {
      source_shape = 
        (options[std::string("N") + species_for_source_shape]["source_shape"]
          .doc("Source term in ddt(N" + species_for_source_shape + std::string("). Units [m^-3/s]"))
          .withDefault(Field3D(0.0))
        ) / (Nnorm * Omega_ci);

        source_units = "m^-3 / s";
        source_conversion = Nnorm * Omega_ci;

    }

    diagnose = detachment_controller_options["diagnose"]
                  .doc("Output additional diagnostics?")
                  .withDefault<bool>(false);
    
  }

  void transform(Options& state) override;

  void outputVars(Options& state) override {
    AUTO_TRACE();
    if (diagnose) {
      set_with_attrs(
        state[std::string("detachment_front_location")], detachment_front_location,
        {{"units", "m"},
         {"time_dimension", "t"},
         {"long_name", "detachment front position"},
         {"source", "detachment_controller"}});
      
      // Shape is not time-dependent and has units
      set_with_attrs(
          state[std::string("detachment_control_src_shape")], source_shape,
          {{"units", source_units},
           {"conversion", source_conversion},
           {"long_name", "detachment control source shape"},
           {"source", "detachment_controller"}});

      // The source multiplier is time-dependent, but dimensionless
      // because all the units are attached to the shape
      set_with_attrs(state[std::string("detachment_control_src_mult")],
                     source_multiplier,
                     {{"time_dimension", "t"},
                      {"long_name", "detachment control source multiplier"},
                      {"source", "detachment_controller"}});

      set_with_attrs(state[std::string("detachment_source_feedback")], detachment_source_feedback,
                      {{"time_dimension", "t"},
                      {"units", source_units},
                      {"conversion", source_conversion},
                      {"standard_name", "detachment control source"},
                      {"long_name", "detachment control source"},
                      {"source", "detachment_controller"}});

      // Save proportional and integral component of source for diagnostics/tuning
      // Multiplier = proportional term + integral term
      set_with_attrs(
          state[std::string("detachment_control_src_p")], proportional_term,
          {{"time_dimension", "t"},
          {"long_name", "proportional feedback term"},
          {"source", "detachment_controller"}});

      
      set_with_attrs(
          state[std::string("detachment_control_src_i")], integral_term,
          {{"time_dimension", "t"},
           {"long_name", "integral feedback term"},
           {"source", "detachment_controller"}});

    }
  }

  void restartVars(Options& state) override {
    AUTO_TRACE();

    // NOTE: This is a hack because we know that the loaded restart file
    //       is passed into restartVars in PhysicsModel::postInit
    // The restart value should be used in init() rather than here
    static bool first = true;
    if (first and state.isSet("detachment_control_error_integral")) {
      first = false;
      error_integral = state["detachment_control_error_integral"].as<BoutReal>();
    }

    // Save the detachment control error integral
    set_with_attrs(state["detachment_control_error_integral"], error_integral,
                   {{"long_name", "detachment control error integral"},
                    {"source", "detachment_controller"}});
  }

private:
  std::list<std::string> species_list; ///< Which species to apply the factor to
  std::list<std::string> scaling_factors_list; ///< Factor to apply

  BoutReal detachment_front_desired_location;
  bool error_on_log_distance;
  bool exponential_control;
  bool control_power;
  BoutReal minval_for_log;

  std::string species_for_source_shape;
  std::string neutral_species;
  
  bool force_integral_positive; ///< Force integral term to be positive?
  bool force_source_positive;   ///< Force source to be positive?

  BoutReal connection_length;
  BoutReal controller_p, controller_i; ///< PI controller parameters
  BoutReal error;
  BoutReal error_integral{0.0}; ///< Time integral of the error

  // Terms used in Trapezium rule integration of error
  BoutReal error_lasttime{-1.0};
  BoutReal error_last{0.0};

  Field3D source_shape; ///< This shape source is scaled up and down
  Field3D detachment_source_feedback;
  std::string source_units;
  BoutReal source_conversion;

  BoutReal detachment_front_location{0.0};
  int detachment_front_index{0};
  BoutReal source_multiplier; ///< Factor to multiply source
  BoutReal proportional_term, integral_term; ///< Components of resulting source for diagnostics

  bool diagnose; ///< Output diagnostic information?
};

namespace {
  RegisterComponent<DetachmentController> register_detachment_controller("detachment_controller");
}

#endif // detachment_controller_H
