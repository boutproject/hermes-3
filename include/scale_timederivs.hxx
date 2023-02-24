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
  ScaleTimeDerivs(std::string name, Options &alloptions, Solver*) {
    Options &options = alloptions[name];

  }

  /// Sets in the state
  ///
  /// - solver
  ///   - scale_timederivs
  ///
  void transform(Options &state) override {
    
    auto* coord = bout::globals::mesh->getCoordinates();
    Field2D dl2 = coord->g_22 * SQ(coord->dy);
    
    // Scale by parallel heat conduction CFL timescale
    auto Te = get<Field3D>(state["species"]["e"]["temperature"]);
    Field3D dt = dl2 / pow(Te, 5./2);

    state["scale_timederivs"] = dt / max(dt, true);
  }

private:
};

namespace {
RegisterComponent<ScaleTimeDerivs> registercomponentscaletimederivs("scale_timederivs");
}

#endif // SCALE_TIMEDERIVS_H

