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
    BoutReal Omega_ci = 1.0 / get<BoutReal>(units["seconds"]);
    time_normalisation = 1.0 / Omega_ci;

    BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

    connection_length = options["mesh"]["length"].as<BoutReal>();

    detachment_front_setpoint = 
      detachment_controller_options["detachment_front_setpoint"]
      .doc("How far from the divertor target should the detachment front be? (in m).")
      .as<BoutReal>();
    
    min_time_for_change = 
      detachment_controller_options["min_time_for_change"]
      .doc("Minimum time change before changing the control signal.")
      .withDefault<BoutReal>(1E-12);
    
    min_error_for_change = 
      detachment_controller_options["min_error_for_change"]
      .doc("Minimum error change before changing the control signal.")
      .withDefault<BoutReal>(0.0);
    
    minval_for_source_multiplier = 
      detachment_controller_options["minval_for_source_multiplier"]
      .doc("Minimum value for the control signal.")
      .withDefault<BoutReal>(-INFINITY);
    
    maxval_for_source_multiplier = 
      detachment_controller_options["maxval_for_source_multiplier"]
      .doc("Maximum value for the control signal.")
      .withDefault<BoutReal>(INFINITY);
    
    species_for_source_shape =
      detachment_controller_options["species_for_source_shape"]
      .doc("Which species to select the source_shape from?")
      .withDefault<std::string>("e");
    
    neutral_species =
      detachment_controller_options["neutral_species"]
      .doc("Which is the main neutral species?")
      .as<std::string>();
    
    actuator = 
      detachment_controller_options["actuator"]
      .doc("What should we adjust to control the detachment front position (options are 'power' or 'particles' (either main species or impurity))?")
      .as<std::string>(); 
    if (actuator == "power") {
      control_mode = control_power;
      response_sign = -1.0;
    } else if (actuator == "particles") {
      control_mode = control_particles;
      response_sign = 1.0;
    } else {
      // Invalid control mode
      ASSERT2(false);
    }
    
    control = detachment_controller_options["initial_source_multiplier"]
      .doc("Initial source multiplier for controller.")
      .withDefault<BoutReal>(1.0);
    previous_control = control;

    ignore_restart = detachment_controller_options["ignore_restart"]
                               .doc("Ignore the restart file (mainly useful for development).")
                               .withDefault(false);

    controller_gain = detachment_controller_options["controller_gain"]
                               .doc("Detachment controller gain (Kc parameter)")
                               .withDefault(0.0);
    integral_time = detachment_controller_options["integral_time"]
                               .doc("Detachment controller integral time")
                               .withDefault(INFINITY);
    derivative_time = detachment_controller_options["derivative_time"]
                               .doc("Detachment controller detachment time")
                               .withDefault(0.0);
    
    alpha_de = detachment_controller_options["alpha_de"]
                               .doc("Smoothing factor applied to the change_in_error (low pass filter, must be between 0 and 1)")
                               .withDefault(0.0);
    ASSERT2(((alpha_de >= 0.0) && (alpha_de <= 1.0)));
    
    alpha_d2e = detachment_controller_options["alpha_d2e"]
                               .doc("Smoothing factor applied to the change_in_derivative (low pass filter, must be between 0 and 1)")
                               .withDefault(0.0);
    ASSERT2(((alpha_d2e >= 0.0) && (alpha_d2e <= 1.0)));

    species_list = strsplit(detachment_controller_options["species_list"]
                                  .doc("Comma-separated list of species to apply the PI-controlled source to")
                                  .as<std::string>(),
                            ',');
    scaling_factors_list = strsplit(detachment_controller_options["scaling_factors_list"]
                                  .doc("Comma-separated list of scaling factors to apply to the PI-controlled source, 1 for each species")
                                  .as<std::string>(),
                            ',');

    if (species_list.size() != scaling_factors_list.size()) {
      throw BoutException("DetachmentController: species_list length doesn't match scaling_factors length. Need 1 scaling factor per species.");
    }

    if (control_mode == control_power) {
      source_shape =
        (options[std::string("P") + species_for_source_shape]["source_shape"]
          .doc("Source term in ddt(P" + species_for_source_shape + std::string("). Units [Pa/s], note P = 2/3 E."))
          .withDefault(Field3D(0.0))
        )  / (Pnorm * Omega_ci);

      source_units = "Pa / s";
      source_conversion = Pnorm * Omega_ci;
    } else if (control_mode == control_particles) {
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
    
    debug =
      detachment_controller_options["debug"]
      .doc("Print debugging information to the screen (0 for none, 1 for basic, 2 for extensive).")
      .withDefault<int>(0);
    
  };

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
      set_with_attrs(state[std::string("detachment_control_src_mult")], control,
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
      
      set_with_attrs(state[std::string("detachment_control_change_in_error")], change_in_error,
                    {{"time_dimension", "t"},
                    {"standard_name", "change_in_error"},
                    {"long_name", "detachment control change_in_error"},
                    {"source", "detachment_controller"}});
      set_with_attrs(state[std::string("detachment_control_change_in_time")], change_in_time,
                    {{"time_dimension", "t"},
                    {"standard_name", "change_in_time"},
                    {"long_name", "detachment control change_in_time"},
                    {"source", "detachment_controller"}});
      set_with_attrs(state[std::string("detachment_control_error")], error,
                    {{"time_dimension", "t"},
                    {"standard_name", "error"},
                    {"long_name", "detachment control error"},
                    {"source", "detachment_controller"}});
      set_with_attrs(state[std::string("detachment_control_derivative")], derivative,
                    {{"time_dimension", "t"},
                    {"standard_name", "derivative"},
                    {"long_name", "detachment control derivative"},
                    {"source", "detachment_controller"}});
      set_with_attrs(state[std::string("detachment_control_change_in_derivative")], change_in_derivative,
                    {{"time_dimension", "t"},
                    {"standard_name", "change_in_derivative"},
                    {"long_name", "detachment control change_in_derivative"},
                    {"source", "detachment_controller"}});
      set_with_attrs(state[std::string("detachment_control_change_in_control")], change_in_control,
                    {{"time_dimension", "t"},
                    {"standard_name", "change_in_control"},
                    {"long_name", "detachment control change_in_control"},
                    {"source", "detachment_controller"}});
      set_with_attrs(state[std::string("detachment_control_control")], control,
                     {{"time_dimension", "t"},
                      {"long_name", "detachment control response"},
                      {"source", "detachment_controller"}});
  }}

  void restartVars(Options& state) override {
    AUTO_TRACE();
    
    if ((initialise) && (not ignore_restart)) {
      if (state.isSet("detachment_control_src_mult")) {
        control = state["detachment_control_src_mult"].as<BoutReal>();
      }

      if (state.isSet("detachment_control_previous_control")) {
        previous_control = state["detachment_control_previous_control"].as<BoutReal>();
      }
      if (state.isSet("detachment_control_previous_time")) {
        previous_time = state["detachment_control_previous_time"].as<BoutReal>();
      }
      if (state.isSet("detachment_control_previous_error")) {
        previous_error = state["detachment_control_previous_error"].as<BoutReal>();
      }
      if (state.isSet("detachment_control_previous_derivative")) {
        previous_derivative = state["detachment_control_previous_derivative"].as<BoutReal>();
      }
      if (state.isSet("detachment_control_previous_change_in_error")) {
        previous_change_in_error = state["detachment_control_previous_change_in_error"].as<BoutReal>();
      }
      if (state.isSet("detachment_control_previous_change_in_derivative")) {
        previous_change_in_derivative = state["detachment_control_previous_change_in_derivative"].as<BoutReal>();
      }

      initialise = false;
    }
    
    set_with_attrs(state["detachment_control_src_mult"], control, {{"source", "detachment_controller"}});
    set_with_attrs(state["detachment_control_previous_control"], previous_control, {{"source", "detachment_controller"}});
    set_with_attrs(state["detachment_control_previous_time"], previous_time, {{"source", "detachment_controller"}});
    set_with_attrs(state["detachment_control_previous_error"], previous_error, {{"source", "detachment_controller"}});
    set_with_attrs(state["detachment_control_previous_derivative"], previous_derivative, {{"source", "detachment_controller"}});
    set_with_attrs(state["detachment_control_previous_change_in_error"], previous_change_in_error, {{"source", "detachment_controller"}});
    set_with_attrs(state["detachment_control_previous_change_in_derivative"], previous_change_in_derivative, {{"source", "detachment_controller"}});
  }

  private:
    // Inputs variables
    BoutReal connection_length;
    BoutReal detachment_front_setpoint;
    BoutReal min_time_for_change;
    BoutReal min_error_for_change;
    std::string species_for_source_shape;
    std::string neutral_species;
    std::string actuator;
    bool ignore_restart;
    BoutReal response_sign;
    BoutReal controller_gain;
    BoutReal integral_time;
    BoutReal derivative_time;
    std::list<std::string> species_list; ///< Which species to apply the factor to
    std::list<std::string> scaling_factors_list; ///< Factor to apply
    Field3D source_shape;
    std::string source_units;
    BoutReal source_conversion;
    bool diagnose;
    BoutReal minval_for_source_multiplier;
    BoutReal maxval_for_source_multiplier;
    BoutReal alpha_de;
    BoutReal alpha_d2e;
    int debug;

    int control_mode;
    int control_power{0};
    int control_particles{1};

    // System state variables for output
    Field3D detachment_source_feedback;
    BoutReal detachment_front_location;
    BoutReal control;
    BoutReal change_in_error;
    BoutReal change_in_time;

    BoutReal time;
    BoutReal error;
    BoutReal derivative;
    BoutReal change_in_derivative;
    BoutReal change_in_control;

    // Private system state variables for calculations
    BoutReal previous_control;
    BoutReal previous_time{0.0};
    BoutReal previous_error{0.0};
    BoutReal previous_derivative{0.0};
    BoutReal previous_change_in_error{0.0};
    BoutReal previous_change_in_derivative{0.0};
    BoutReal time_normalisation;
    bool initialise{true};
    bool evaluate_derivatives{false};

};

namespace {
  RegisterComponent<DetachmentController> register_detachment_controller("detachment_controller");
}

#endif // detachment_controller_H
