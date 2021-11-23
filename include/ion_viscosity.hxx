#pragma once
#ifndef ION_VISCOSITY_H
#define ION_VISCOSITY_H

#include "component.hxx"

/// Ion viscosity terms
///
/// Adds a viscosity to all species which are not electrons
///
/// Uses Braginskii collisional form, combined with a
/// SOLPS-like flux limiter.
///
/// Needs to be calculated after collisions, because collision
/// frequency is used to calculate parallel viscosity
struct IonViscosity : public Component {
  /// Inputs
  /// - <name>
  ///   - eta_limit_alpha    Flux limiter coefficient
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
  BoutReal eta_limit_alpha; ///< Flux limit coefficient
};

namespace {
RegisterComponent<IonViscosity> registercomponentionviscosity("ion_viscosity");
}

#endif
