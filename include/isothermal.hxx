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
  BoutReal T; ///< The normalised temperature
  Field3D P;  ///< The normalised pressure

  bool initialize_from_mesh; ///< Initilize the Field3D T from 2D profiles stored in the
                             ///< mesh file.

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
