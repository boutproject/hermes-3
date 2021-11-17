#pragma once
#ifndef ION_VISCOSITY_H
#define ION_VISCOSITY_H

#include "component.hxx"

/// Ion viscosity terms
///
/// Adds a viscosity to all species which are not electrons
///
/// Needs to be calculated after collisions, because collision
/// frequency is used to calculate parallel viscosity
struct IonViscosity : public Component {
  IonViscosity(std::string name, Options& alloptions, Solver*);

  /// Inputs
  /// - species
  ///   - <name>   (skips "e")
  ///     - pressure  (skips if not present)
  ///     - velocity  (skips if not present)
  ///     - collision_frequency
  ///
  /// Sets in the state
  /// - species
  ///   - <name>
  ///     - momentum_source
  ///
  void transform(Options &state) override;
private:

};

namespace {
RegisterComponent<IonViscosity> registercomponentionviscosity("ion_viscosity");
}

#endif
