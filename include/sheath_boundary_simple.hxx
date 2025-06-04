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
  /// # Input options
  /// - <name>  e.g. "sheath_boundary_simple"
  ///   - lower_y                  Boundary on lower y?
  ///   - upper_y                  Boundary on upper y?
  ///   - gamma_e                  Electron sheath heat transmission coefficient
  ///   - gamma_i                  Ion sheath heat transmission coefficient
  ///   - sheath_ion_polytropic    Ion polytropic coefficient in Bohm sound speed. Default 1.
  ///   - wall_potential           Voltage of the wall [Volts]
  ///   - secondary_electron_coef  Effective secondary electron emission coefficient
  ///   - sin_alpha                Sine of the angle between magnetic field line and wall surface (0 to 1)
  ///   - always_set_phi           Always set phi field? Default is to only modify if already set
  SheathBoundarySimple(std::string name, Options &options, Solver *);

  ///
  /// # Inputs
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
  /// # Outputs
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
  void outputVars(Options &state) override;
private:
  BoutReal Ge; // Secondary electron emission coefficient
  BoutReal sin_alpha; // sin of angle between magnetic field and wall.

  BoutReal gamma_e; ///< Electron sheath heat transmission
  BoutReal gamma_i; ///< Ion sheath heat transmission
  BoutReal sheath_ion_polytropic; ///< Polytropic coefficient in sheat velocity
  
  bool lower_y; // Boundary on lower y?
  bool upper_y; // Boundary on upper y?

  bool always_set_phi; ///< Set phi field?

  Field3D wall_potential; ///< Voltage of the wall. Normalised units.

  Field3D hflux_e;  // Electron heat flux through sheath
  Field3D phi; // Phi at sheath
  Field3D ion_sum; // Sum of ion current at sheath

  bool diagnose; // Save diagnostic variables?
  Options diagnostics;   // Options object to store diagnostic fields like a dict

  bool no_flow; ///< No flow speed, only remove energy

  BoutReal density_boundary_mode, pressure_boundary_mode, temperature_boundary_mode; ///< BC mode: 0=LimitFree, 1=ExponentialFree, 2=LinearFree
};

namespace {
RegisterComponent<SheathBoundarySimple>
    registercomponentsheathboundarysimple("sheath_boundary_simple");
}

#endif // SHEATH_BOUNDARY_SIMPLE_H