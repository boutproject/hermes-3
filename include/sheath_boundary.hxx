#pragma once
#ifndef SHEATH_BOUNDARY_H
#define SHEATH_BOUNDARY_H

#include "component.hxx"

/// Boundary condition at the wall in Y
///
/// These are based on
/// "Boundary conditions for the multi-ion magnetized plasma-wall transition"
///  by D.Tskhakaya, S.Kuhn. JNM 337-339 (2005), 405-409
///
/// Note: The approximation used here is for ions having similar
/// gyro-orbit sizes
struct SheathBoundary : public Component {
  SheathBoundary(std::string name, Options &options, Solver *);

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
};

namespace {
RegisterComponent<SheathBoundary>
    registercomponentsheathboundary("sheath_boundary");
}

#endif // SHEATH_BOUNDARY_H
