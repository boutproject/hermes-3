#pragma once
#ifndef SCALE_TIMEDERIVS_H
#define SCALE_TIMEDERIVS_H

#include "component.hxx"
#include <globals.hxx>

/// Scale time derivatives of the system
///
/// This is intended for steady-state calculations
/// where the aim is to reach ddt -> 0
///
struct ScaleTimeDerivs : public Component {
  ScaleTimeDerivs(std::string, Options&, Solver*) {}

  /// Sets in the state
  ///
  /// - scale_timederivs
  ///
  void transform(Options &state) override {

    auto* coord = bout::globals::mesh->getCoordinates();
    Field2D dl2 = coord->g_22 * SQ(coord->dy);

    // Scale by parallel heat conduction CFL timescale
    auto Te = get<Field3D>(state["species"]["e"]["temperature"]);
    Field3D dt = dl2 / pow(floor(Te, 1e-5), 5./2);
    scaling = dt / max(dt, true); // Saved for output

    state["scale_timederivs"] = scaling;
  }

  void outputVars(Options& state) override {
    set_with_attrs(
        state["scale_timederivs"], scaling,
        {{"time_dimension", "t"},
         {"long_name", "Scaling factor applied to all time derivatives"},
         {"source", "scale_timederivs"}});
  }
private:
  Field3D scaling; // The scaling factor applied to each cell
};

namespace {
RegisterComponent<ScaleTimeDerivs> registercomponentscaletimederivs("scale_timederivs");
}

#endif // SCALE_TIMEDERIVS_H

