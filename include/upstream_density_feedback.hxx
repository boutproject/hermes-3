#pragma once
#ifndef UPSTREAM_DENSITY_FEEDBACK_H
#define UPSTREAM_DENSITY_FEEDBACK_H

#include "component.hxx"

struct UpstreamDensityFeedback : public Component {
  UpstreamDensityFeedback(std::string name, Options& alloptions, Solver*) : name(name) {
    const auto& units = alloptions["units"];
    BoutReal Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
    BoutReal FreqNorm = 1. / get<BoutReal>(units["seconds"]);

    Options& options = alloptions[name];

    density_upstream =
        options["density_upstream"].doc("Upstream density (at y=0) [m^-3]").as<BoutReal>()
        / Nnorm;

    density_controller_p = options["density_controller_p"]
                               .doc("Feedback controller proportional (p) parameter")
                               .withDefault(1e-2);
    density_controller_i = options["density_controller_i"]
                               .doc("Feedback controller integral (i) parameter")
                               .withDefault(1e-3);

    density_integral_positive = options["density_integral_positive"]
                                    .doc("Force integral term to be positive?")
                                    .withDefault<bool>(false);
    density_source_positive = options["density_source_positive"]
                                  .doc("Force source to be positive?")
                                  .withDefault<bool>(true);

    // NOTE: density_error_integral should be set from the restart file here.
    // There doesn't seem to be a way to get the restart Options,
    // so for now this is set in restartVars

    // Source shape the same as used in EvolveDensity
    density_source_shape =
        alloptions[std::string("N") + name]["source"]
            .doc("Source term in ddt(N" + name + std::string("). Units [m^-3/s]"))
            .withDefault(Field3D(0.0))
        / (Nnorm * FreqNorm);

    diagnose = options["diagnose"]
                   .doc("Output additional diagnostics?")
                   .withDefault<bool>(false);
  }

  /// Inputs
  ///  - <name>
  ///    - density
  ///
  /// Outputs
  ///
  ///  - <name>
  ///    - density_source
  ///
  void transform(Options& state) override;

  void outputVars(Options& state) override {
    AUTO_TRACE();
    if (diagnose) {
      // Normalisations
      auto Nnorm = get<BoutReal>(state["Nnorm"]);
      auto Omega_ci = get<BoutReal>(state["Omega_ci"]);

      // Shape is not time-dependent
      set_with_attrs(
          state[std::string("density_source_shape_") + name], density_source_shape,
          {{"units", "m^-3 s^-1"},
           {"conversion", Nnorm * Omega_ci},
           {"long_name", name + " density source shape"},
           {"source", "upstream_density_feedback"}});

      // The source multiplier is time-dependent, but dimensionless
      // because all the units are attached to the shape
      set_with_attrs(state[std::string("density_source_multiplier_") + name],
                     source_multiplier,
                     {{"time_dimension", "t"},
                      {"long_name", name + " density source multiplier"},
                      {"source", "upstream_density_feedback"}});
    }
  }

  void restartVars(Options& state) override {
    AUTO_TRACE();

    // NOTE: This is a hack because we know that the loaded restart file
    //       is passed into restartVars in PhysicsModel::postInit
    // The restart value should be used in init() rather than here
    static bool first = true;
    if (first and state.isSet(name + "_density_error_integral")) {
      first = false;
      density_error_integral = state[name + "_density_error_integral"].as<BoutReal>();
    }

    // Save the density error integral
    set_with_attrs(state[name + "_density_error_integral"], density_error_integral,
                   {{"long_name", name + " density error integral"},
                    {"source", "upstream_density_feedback"}});
  }

private:
  std::string name; ///< The species name

  BoutReal density_upstream;                           ///< Normalised upstream density
  BoutReal density_controller_p, density_controller_i; ///< PI controller parameters

  BoutReal density_error_integral{0.0}; ///< Time integral of the error

  bool density_integral_positive; ///< Force integral term to be positive?
  bool density_source_positive;   ///< Force source to be positive?

  // Terms used in Trapezium rule integration of error
  BoutReal density_error_lasttime{-1.0};
  BoutReal density_error_last{0.0};

  Field3D density_source_shape; ///< This shape source is scaled up and down

  BoutReal source_multiplier; ///< Factor to multiply source

  bool diagnose; ///< Output diagnostic information?
};

namespace {
RegisterComponent<UpstreamDensityFeedback> register_uds("upstream_density_feedback");
}

#endif // UPSTREAM_DENSITY_FEEDBACK_H
