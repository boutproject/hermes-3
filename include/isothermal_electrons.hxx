#pragma once
#ifndef ISOTHERMAL_ELECTRONS_H
#define ISOTHERMAL_ELECTRONS_H

#include "component.hxx"

/// Set electron temperature to a fixed value
///
struct IsothermalElectrons : public Component {
  IsothermalElectrons(std::string name, Options &options, Solver *);

  /// Inputs
  /// - species
  ///   - e
  ///     - density
  ///
  /// Sets in the state
  ///
  /// - species
  ///   - e
  ///     - temperature
  ///     - pressure
  ///
  void transform(Options &state) override;

private:
  BoutReal Te;
};

namespace {
RegisterComponent<IsothermalElectrons>
    registercomponentisothermalelectrons("isothermal_electrons");
}

#endif // ISOTHERMAL_ELECTRONS_H
