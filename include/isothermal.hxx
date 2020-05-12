#pragma once
#ifndef ISOTHERMAL_H
#define ISOTHERMAL_H

#include "component.hxx"

/// Set temperature to a fixed value
///
struct Isothermal : public Component {
  Isothermal(std::string name, Options &options, Solver *);

  /// Inputs
  /// - species
  ///   - <name>
  ///     - density
  ///
  /// Sets in the state
  ///
  /// - species
  ///   - <name>
  ///     - temperature
  ///     - pressure
  ///
  void transform(Options &state) override;

private:
  std::string name; // Species name
  
  BoutReal T; // The normalised tempertature
};

namespace {
RegisterComponent<Isothermal> registercomponentisothermal("isothermal");
}

#endif // ISOTHERMAL_H
