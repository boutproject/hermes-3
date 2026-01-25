#pragma once
#ifndef BRAGINSKII_HEAT_EXCHANGE_H
#define BRAGINSKII_HEAT_EXCHANGE_H

#include <string>

#include <bout/options.hxx>

#include "component.hxx"

/// Calculates the heat exchange between species due to collisions
///
struct BraginskiiHeatExchange : public Component {
  ///
  /// @param alloptions Settings. There is nothing to be configured.
  ///
  BraginskiiHeatExchange(const std::string& name, Options& alloptions, Solver*);

  /// Add extra fields for output, or set attributes e.g docstrings
  void outputVars(Options& state) override;

private:
  /// Calculated energy transfer for post-processing and use by other components
  /// Saved in options, the BOUT++ dictionary-like object
  Options energy_channels;

  /// Save more diagnostics?
  bool diagnose;

  /// Calculate thermal energy exchange between species due to collisions.
  ///
  /// Uses
  ///   - species
  ///     - <name>
  ///       - AA
  ///       - charge (if set)
  ///       - collision_frequencies (if section)
  ///       - density
  ///       - temperature (if set)
  ///
  /// Modifies
  ///   - species
  ///     - <name>
  ///       - momentum_source   if species1 or species2 temperature is set
  ///       - energy_source     if temperature is set
  ///
  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<BraginskiiHeatExchange>
    registercomponentbraginskiiheatexchange("braginskii_heat_exchange");
}

#endif // BRAGINSKII_HEAT_EXCHANGE_H
