#pragma once
#ifndef BRAGINSKII_HEAT_EXCHANGE_H
#define BRAGINSKII_HEAT_EXCHANGE_H

#include <bout/field3d.hxx>

#include "component.hxx"

/// Calculates the heat exchange between species due to collisions
///
struct BraginskiiHeatExchange : public Component {
  ///
  /// @param alloptions Settings. There is nothing to be configured.
  ///
  BraginskiiHeatExchange(std::string name, Options& alloptions, Solver*);

  /// Calculate thermal energy exchange between species due to collisions.
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
  void transform(Options& state) override;

  /// Add extra fields for output, or set attributes e.g docstrings
  void outputVars(Options& state) override;

private:
  /// Calculated energy transfer for post-processing and use by other components
  /// Saved in options, the BOUT++ dictionary-like object
  Options energy_channels;

  /// Save more diagnostics?
  bool diagnose;
};

namespace {
RegisterComponent<BraginskiiHeatExchange>
    registercomponentbraginskiiheatexchange("braginskii_heat_exchange");
}

#endif // BRAGINSKII_HEAT_EXCHANGE_H
