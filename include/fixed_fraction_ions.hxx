#pragma once
#ifndef FIXED_FRACTION_IONS_H
#define FIXED_FRACTION_IONS_H

#include "component.hxx"

/// Set ion densities from electron densities
/// 
struct FixedFractionIons : public Component {
  /// Inputs
  /// - <name>
  ///   - fractions   A comma-separated list of pairs separated by @
  ///                 e.g. 'd+ @ 0.5, t+ @ 0.5'
  FixedFractionIons(std::string name, Options &options, Solver *UNUSED(solver));

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
  void transform(Options &state) override;

 private:
  std::vector<std::pair<std::string, BoutReal>> fractions;
};

namespace {
RegisterComponent<FixedFractionIons> registercomponentfixedfractionions("fixed_fraction_ions");
}

#endif // FIXED_FRACTION_IONS_H
