#pragma once
#ifndef ISOTHERMAL_H
#define ISOTHERMAL_H

#include "component.hxx"

/// Set temperature to a fixed value
///
struct Isothermal : public NamedComponent<Isothermal> {
  Isothermal(std::string name, Options& options, Solver*);

  void outputVars(Options& state) override;

  static constexpr auto type = "isothermal";

private:
  std::string name; // Species name
  bool dipole_scaling; ///< Apply dipole scaling ~ B to the temperature
  BoutReal T; ///< The normalised temperature
  Field3D P; ///< The normalised pressure
  Field2D T2D, B2D, B2D_edge; ///< The normalised temperature, averaged to 2D (for dipole scaling)
  bool temperature_3D; ///< Should the temperature be a 3D field (true) or a constant (false)?
  bool diagnose; ///< Output additional diagnostics?

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
  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<Isothermal> registercomponentisothermal;
}

#endif // ISOTHERMAL_H
