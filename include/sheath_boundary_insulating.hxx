#pragma once
#ifndef SHEATH_BOUNDARY_INSULATING_H
#define SHEATH_BOUNDARY_INSULATING_H

#include "component.hxx"

/// Insulating sheath boundary condition at the wall in Y
///
/// This is a collective component, because it couples all charged species
///
/// Adapted from the `sheath_boundary` component, but always sets the current
/// density to zero
struct SheathBoundaryInsulating : public Component {
  SheathBoundaryInsulating(std::string name, Options &options, Solver *);

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
  
  bool lower_y; // Boundary on lower y?
  bool upper_y; // Boundary on upper y?

  BoutReal gamma_e; ///< Electron sheath heat transmission
};

namespace {
RegisterComponent<SheathBoundaryInsulating>
    registercomponentsheathboundaryinsulating("sheath_boundary_insulating");
}

#endif // SHEATH_BOUNDARY_INSULATING_H
