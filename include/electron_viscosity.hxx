#pragma once
#ifndef ELECTRON_VISCOSITY_H
#define ELECTRON_VISCOSITY_H

#include "component.hxx"

/// Electron viscosity 
///
/// Adds Braginskii parallel electron viscosity, with SOLPS-style
/// viscosity flux limiter
///
/// Needs to be calculated after collisions, because collision
/// frequency is used to calculate parallel viscosity
///
/// References
///  - https://farside.ph.utexas.edu/teaching/plasma/lectures1/node35.html
///
struct ElectronViscosity : public Component {
  /// Inputs
  /// - <name>
  ///   - diagnose: bool, default false
  ///     Output diagnostic SNVe_viscosity?
  ///   - eta_limit_alpha: float, default -1.0
  ///     Flux limiter coefficient. < 0 means no limiter
  ElectronViscosity(std::string name, Options& alloptions, Solver*);

  /// Inputs
  /// - species
  ///   - e
  ///     - pressure  (skips if not present)
  ///     - velocity  (skips if not present)
  ///     - collision_frequency
  ///
  /// Sets in the state
  /// - species
  ///   - e
  ///     - momentum_source
  ///
  void transform(Options &state) override;

  void outputVars(Options &state) override;
private:
  BoutReal eta_limit_alpha; ///< Flux limit coefficient
  bool diagnose; ///< Output viscosity diagnostic?
  Field3D viscosity; ///< The viscosity momentum source
};

namespace {
RegisterComponent<ElectronViscosity> registercomponentelectronviscosity("electron_viscosity");
}

#endif
