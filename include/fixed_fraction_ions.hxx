#pragma once
#ifndef FIXED_FRACTION_IONS_H
#define FIXED_FRACTION_IONS_H

#include "component.hxx"

/// Set ion densities from electron densities
///
struct FixedFractionIons : public NamedComponent<FixedFractionIons> {
  /// Inputs
  /// - <name>
  ///   - fractions   A comma-separated list of pairs separated by @
  ///                 e.g. 'd+ @ 0.5, t+ @ 0.5'
  FixedFractionIons(std::string name, Options& options, Solver* UNUSED(solver));

  static constexpr auto type = "fixed_fraction_ions";

private:
  std::vector<std::pair<std::string, BoutReal>> fractions;

  /// Required inputs
  ///
  /// - species
  ///   - e
  ///     - density
  ///
  /// Sets in the state the density of each species
  ///
  /// - species
  ///   - <species1>
  ///     - density  = <fraction1> * electron density
  ///   - ...
  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<FixedFractionIons> registercomponentfixedfractionions;
}

#endif // FIXED_FRACTION_IONS_H
