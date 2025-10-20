#pragma once
#ifndef ISOTHERMAL_H
#define ISOTHERMAL_H

#include "component.hxx"

/// Set temperature to a fixed value
///
struct Isothermal : public Component {
  Isothermal(std::string name, Options& options, Solver*);

  /// Inputs
  /// - species
  ///   - <name>
  ///     - density (optional)
  ///
  /// Sets in the state
  ///
  /// - species
  ///   - <name>
  ///     - temperature
  ///     - pressure (if density is set)
  ///
  void transform(Options& state) override;

  void outputVars(Options& state) override;

private:
  std::string name; // Species name

  Field3D T; ///< The normalised temperature
  Field3D P; ///< The normalised pressure

  bool diagnose; ///< Output additional diagnostics?
};

namespace {
RegisterComponent<Isothermal> registercomponentisothermal("isothermal");
}

#endif // ISOTHERMAL_H
