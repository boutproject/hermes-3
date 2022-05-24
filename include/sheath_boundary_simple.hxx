#pragma once
#ifndef SHEATH_BOUNDARY_SIMPLE_H
#define SHEATH_BOUNDARY_SIMPLE_H

#include "component.hxx"

/// Boundary condition at the wall in Y
///
/// This is a collective component, because it couples all charged species
///
/// This implements a simple boundary condition, where each species
/// goes to their own sound velocity at the sheath entrance.
///
/// Notes:
///   - It is recommended to use SheathBoundary rather than SheathBoundarySimple;
///     this is here for comparison to that more complete model.
///
struct SheathBoundarySimple : public Component {
  SheathBoundarySimple(std::string name, Options &options, Solver *);

  ///
  /// Inputs
  /// - species
  ///   - e
  ///     - density
  ///     - temperature
  ///     - pressure    Optional
  ///     - velocity    Optional
  ///     - mass        Optional
  ///     - adiabatic   Optional. Ratio of specific heats, default 5/3.
  ///   - <ions>  if charge is set (i.e. not neutrals)
  ///     - charge
  ///     - mass
  ///     - density
  ///     - temperature
  ///     - pressure     Optional
  ///     - velocity     Optional. Default 0
  ///     - momentum     Optional. Default mass * density * velocity
  ///     - adiabatic    Optional. Ratio of specific heats, default 5/3.
  /// - fields
  ///   - phi    Optional. If not set, calculated at boundary (see note below)
  ///
  /// Outputs
  /// - species
  ///   - e
  ///     - density      Sets boundary
  ///     - temperature  Sets boundary
  ///     - velocity     Sets boundary
  ///     - energy_source
  ///   - <ions>
  ///     - density      Sets boundary
  ///     - temperature  Sets boundary
  ///     - velocity     Sets boundary
  ///     - momentum     Sets boundary
  ///     - energy_source
  /// - fields
  ///   - phi   Sets boundary
  ///
  /// If the field phi is set, then this is used in the boundary condition.
  /// If not set, phi at the boundary is calculated and stored in the state.
  /// Note that phi in the domain will not be set, so will be invalid data.
  ///
  ///
  void transform(Options &state) override;
private:
  BoutReal Ge; // Secondary electron emission coefficient
  BoutReal sin_alpha; // sin of angle between magnetic field and wall.

  BoutReal gamma_e; ///< Electron sheath heat transmission
  BoutReal gamma_i; ///< Ion sheath heat transmission
  BoutReal sheath_ion_polytropic; ///< Polytropic coefficient in sheat velocity
  
  bool lower_y; // Boundary on lower y?
  bool upper_y; // Boundary on upper y?

  bool always_set_phi; ///< Set phi field?
};

namespace {
RegisterComponent<SheathBoundarySimple>
    registercomponentsheathboundarysimple("sheath_boundary_simple");
}

#endif // SHEATH_BOUNDARY_SIMPLE_H
