#pragma once
#ifndef BRAGINSKII_FRICTION_H
#define BRAGINSKII_FRICTION_H

#include <bout/field3d.hxx>

#include "component.hxx"

/// Calculates the friction force applied to each species due to collisions.
/// 
struct BraginskiiFriction : public Component {
  ///
  /// @param alloptions Settings, which should include:
  ///    - units
  ///      - eV
  ///      - inv_meters_cubed
  ///      - meters
  ///      - seconds
  ///
  /// There are switches for other terms:
  ///
  ///   - frictional_heating    Include R dot v heating term as energy source? (includes Ohmic heating)
  ///
  BraginskiiFriction(std::string name, Options& alloptions, Solver*);

  void transform(Options &state) override;

  /// Add extra fields for output, or set attributes e.g docstrings
  void outputVars(Options &state) override;

private:
  BoutReal Tnorm; // Temperature normalisation [eV]
  BoutReal Nnorm; // Density normalisation [m^-3]
  BoutReal rho_s0;  // Length normalisation [m]
  BoutReal Omega_ci; // Frequency normalisation [s^-1]

  /// Include frictional heating term?
  bool frictional_heating;

  /// Calculated friction heating and momentum rates saved for post-processing and use by other components
  /// Saved in options, the BOUT++ dictionary-like object
  Options friction_energy_channels, momentum_channels;

  /// Save more diagnostics?
  bool diagnose;
};

namespace {
RegisterComponent<BraginskiiFriction> registercomponentcollisions("braginskii_friction");
}

#endif // BRAGINSKII_FRICTION_H
