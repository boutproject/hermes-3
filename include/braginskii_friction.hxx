#pragma once
#ifndef BRAGINSKII_FRICTION_H
#define BRAGINSKII_FRICTION_H

#include <bout/field3d.hxx>

#include "component.hxx"

/// Calculates the friction force applied to each species due to collisions.
/// 
struct BraginskiiFriction : public Component {
  ///
  /// @param alloptions Settings, which has switches for additional terms:
  ///
  ///   - frictional_heating    Include R dot v heating term as energy source? (includes Ohmic heating)
  ///
  BraginskiiFriction(std::string name, Options& alloptions, Solver*);


  /// Calculate transfer of momentum and energy between species due to
  /// friction arising from collisions.
  ///
  /// Uses
  ///   - species
  ///     - <name>
  ///       - AA
  ///       - charge
  ///       - collision_frequencies
  ///       - density
  ///       - velocity
  ///
  /// Modifies
  ///   - species
  ///     - <name>
  ///       - momentum_source   if species1 or species2 velocity is set
  ///       - energy_source     if velocity is set and frictional_heating
  ///
  void transform(Options &state) override;

  /// Add extra fields for output, or set attributes e.g docstrings
  void outputVars(Options &state) override;

private:

  /// Include frictional heating term?
  bool frictional_heating;

  /// Calculated friction heating and momentum rates saved for post-processing and use by other components
  /// Saved in options, the BOUT++ dictionary-like object
  Options friction_energy_channels, momentum_channels;

  /// Save more diagnostics?
  bool diagnose;
};

namespace {
RegisterComponent<BraginskiiFriction> registercomponentbraginskiifriction("braginskii_friction");
}

#endif // BRAGINSKII_FRICTION_H
