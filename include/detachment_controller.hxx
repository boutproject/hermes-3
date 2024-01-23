#pragma once
#ifndef detachment_controller_H
#define detachment_controller_H

#include "component.hxx"
#include <bout/constants.hxx>

struct DetachmentController : public Component {

  DetachmentController(std::string, Options& options, Solver*) {

    Options& detachment_controller_options = options["detachment_controller"]

    const auto& units = options["units"];
    BoutReal Tnorm = get<BoutReal>(units["eV"]);
    BoutReal Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
    BoutReal Omega_ci = 1. / get<BoutReal>(units["seconds"]);

    BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

    connection_length = get<BoutReal>(options["mesh"]["length"]);

    detachment_front_location = 
      detachment_controller_options["detachment_front_location"]
      .doc("Desired location of the detachment front in m, relative to y=0 unless set_location_relative_to_target=true").as<BoutReal>;
    
    set_location_relative_to_target = 
      detachment_controller_options["set_location_relative_to_target"]
      .doc("Set detachment_front_location relative to target (y=-1) if true, else relative to upstream (y=0)")
      .withDefault<bool>(true);
    
    species_for_source_shape =
      detachment_controller_options["species_for_source_shape"]
      .doc("Which species to select the source_shape from?")
      .withDefault<std::string>("e");
    
    neutral_species =
      detachment_controller_options["neutral_species"]
      .doc("Which is the main neutral species?")
      .as<std::string>();
    
    control_power = options["control_power"]
                    .doc("Control power source (if true) or particle source (if false)?")
                    .withDefault<bool>(false);
    
    detachment_controller_p = options["detachment_controller_p"]
                               .doc("Feedback controller proportional (p) parameter")
                               .withDefault(1e-2);
    detachment_controller_i = options["detachment_controller_i"]
                               .doc("Feedback controller integral (i) parameter")
                               .withDefault(1e-3);

    force_integral_positive = options["force_integral_positive"]
                                    .doc("Force integral term to be positive?")
                                    .withDefault<bool>(false);
    force_source_positive = options["force_source_positive"]
                                  .doc("Force source to be positive?")
                                  .withDefault<bool>(true);
    
    species_list = strsplit(options["species_for_detachment_feedback"]
                                  .doc("Comma-separated list of species to apply the PI-controlled source to")
                                  .as<std::string>(),
                            ',');
    scaling_factors_list = strsplit(options["scaling_factors_for_detachment_feedback"]
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
    if control_power {
      source_shape =
        (options[std::string("P") + species_for_source_shape]["source_shape"]
          .doc("Source term in ddt(P" + species_for_source_shape + std::string("). Units [Pa/s], note P = 2/3 E."))
          .as<BoutReal>()
        )  / (Pnorm * Omega_ci);

        source_units = "Pa / s";
        source_conversion = Pnorm * Omega_ci;

    } else {
      source_shape = 
        (options[std::string("N") + species_for_source_shape]["source_shape"]
          .doc("Source term in ddt(N" + species_for_source_shape + std::string("). Units [m^-3/s]"))
          .as<BoutReal>()
        ) / (Nnorm * Omega_ci);

        source_units = "m^-3 / s";
        source_conversion = Pnorm * Omega_ci;

    }

    diagnose = options["diagnose"]
                  .doc("Output additional diagnostics?")
                  .withDefault<bool>(false);
    
  }

  void transform(Options& state) override;

  void outputVars(Options& state) override {
    AUTO_TRACE();
    if (diagnose) {
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

      set_with_attrs(state[std::string("detachment_source_feedback")], source_shape * source_multiplier,
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
      detachment_control_error_integral = state["detachment_control_error_integral"].as<BoutReal>();
    }

    // Save the detachment control error integral
    set_with_attrs(state["detachment_control_error_integral"], detachment_control_error_integral,
                   {{"long_name", "detachment control error integral"},
                    {"source", "detachment_controller"}});
  }

private:
  std::list<std::string> species_list; ///< Which species to apply the factor to
  std::list<std::string> scaling_factors_list; ///< Factor to apply

  BoutReal detachment_front_location;
  bool set_location_relative_to_target;
  bool control_power;
  
  std::string species_for_source_shape;
  std::string neutral_species;
  
  BoutReal detachment_controller_p;
  BoutReal detachment_controller_i;
  bool force_integral_positive; ///< Force integral term to be positive?
  bool force_source_positive;   ///< Force source to be positive?
  BoutReal connection_length;

  BoutReal detachment_control_controller_p, detachment_control_controller_i; ///< PI controller parameters
  BoutReal error;
  BoutReal detachment_control_error_integral{0.0}; ///< Time integral of the error

  // Terms used in Trapezium rule integration of error
  BoutReal detachment_control_error_lasttime{-1.0};
  BoutReal detachment_control_error_last{0.0};

  Field3D source_shape; ///< This shape source is scaled up and down
  std::string source_units;
  BoutReal source_conversion;

  BoutReal source_multiplier; ///< Factor to multiply source

  BoutReal proportional_term, integral_term; ///< Components of resulting source for diagnostics

  bool diagnose; ///< Output diagnostic information?
};

namespace {
RegisterComponent<DetachmentController> register_detachment_controller("detachment_controller");
}

#endif // detachment_controller_H
