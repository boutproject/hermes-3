#pragma once
#ifndef simple_pump_H
#define simple_pump_H

#include "component.hxx"
#include <bout/constants.hxx>

struct SimplePump : public Component {

  SimplePump(std::string name, Options& alloptions, Solver*) : name(name) {

    Options& options = alloptions[name];

    const auto& units = alloptions["units"];
    Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
    Omega_ci = 1.0 / get<BoutReal>(units["seconds"]);

    residence_time = (options["residence_time"]
        .doc("Pumping time-constant. Units [s]")
        .as<BoutReal>()
    ) * Omega_ci;

    sink_shape = (options["sink_shape"]
        .doc("Shape of pumping sink.")
        .withDefault(Field3D(0.0))
    );

    diagnose = options["diagnose"]
                  .doc("Output additional diagnostics?")
                  .withDefault<bool>(false);

    };

    void transform(Options& state) override {

        Field3D species_density = getNoBoundary<Field3D>(state["species"][name]["density"]);

        pumping_sink = (sink_shape * species_density) * (-1.0 / residence_time);

        add(state["species"][name]["density_source"], pumping_sink);

    };

    void outputVars(Options& state) override {
    AUTO_TRACE();
    if (diagnose) {

      set_with_attrs(
          state[std::string("simple_pump_src_shape_" + name)], sink_shape,
          {{"long_name", "simple pump source shape"},
           {"source", "simple_pump"}});
    
      set_with_attrs(
          state[std::string("simple_pump_sink_" + name)], pumping_sink,
          {{"time_dimension", "t"},
           {"units", "m^-3 / s"},
           {"conversion", Nnorm * Omega_ci},
           {"long_name", "simple pump source shape"},
           {"source", "simple_pump"}});

    }}

  private:
    std::string name; ///< The species name
    Field3D sink_shape;
    Field3D pumping_sink;
    BoutReal Nnorm;
    BoutReal Omega_ci;
    BoutReal residence_time;
    bool diagnose;

};

namespace {
  RegisterComponent<SimplePump> register_simple_pump("simple_pump");
}

#endif // simple_pump_H
