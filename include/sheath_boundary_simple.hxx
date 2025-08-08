#pragma once
#ifndef SHEATH_BOUNDARY_SIMPLE_H
#define SHEATH_BOUNDARY_SIMPLE_H

#include "component.hxx"

#include <bout/bout_enum_class.hxx>

#include <cstdint>

BOUT_ENUM_CLASS(SheathLimitMode, limit_free, exponential_free, linear_free);

namespace hermes {
enum class SheathKind : std::uint8_t {
  normal,     // Uses a model for electron heat transmission
  simple,     // User-set value for electron heat transmission
  insulating, // Insulating sheath, free potential BC
};
}

/// Generalised sheath boundary condition at the wall in Y
///
/// This is a collective component, because it couples all charged species
///
///
template <hermes::SheathKind kind>
struct SheathBoundaryBase : public Component {
  /// # Input options
  /// - <name>  e.g. "sheath_boundary_simple"
  ///   - lower_y                  Boundary on lower y?
  ///   - upper_y                  Boundary on upper y?
  ///   - gamma_e                  Electron sheath heat transmission coefficient
  ///   - gamma_i                  Ion sheath heat transmission coefficient
  ///   - sheath_ion_polytropic    Ion polytropic coefficient in Bohm sound speed.
  ///   Default 1.
  ///   - wall_potential           Voltage of the wall [Volts]
  ///   - secondary_electron_coef  Effective secondary electron emission coefficient
  ///   - sin_alpha                Sine of the angle between magnetic field line and wall
  ///   surface (0 to 1)
  ///   - always_set_phi           Always set phi field? Default is to only modify if
  ///   already set
  SheathBoundaryBase(Options& options, BoutReal Tnorm);

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
  void transform(Options& state) override;
  void outputVars(Options& state) override;

private:
  static constexpr hermes::SheathKind sheath_kind = kind;
  static constexpr bool is_insulating() {
    return sheath_kind == hermes::SheathKind::insulating;
  }

  BoutReal Ge;        //< Secondary electron emission coefficient
  BoutReal sin_alpha; //< sin of angle between magnetic field and wall.

  BoutReal gamma_e;               ///< Electron sheath heat transmission
  BoutReal gamma_i;               ///< Ion sheath heat transmission
  BoutReal sheath_ion_polytropic; ///< Polytropic coefficient in sheat velocity

  bool lower_y; // Boundary on lower y?
  bool upper_y; // Boundary on upper y?

  bool always_set_phi; ///< Set phi field?

  Field3D wall_potential; ///< Voltage of the wall. Normalised units.
  bool floor_potential;   ///< Apply floor to sheath potential?

  Field3D hflux_e; // Electron heat flux through sheath
  Field3D phi;     // Phi at sheath

  bool diagnose;       // Save diagnostic variables?
  Options diagnostics; // Options object to store diagnostic fields like a dict

  bool no_flow; ///< No flow speed, only remove energy

  SheathLimitMode density_boundary_mode;     ///< BC for density
  SheathLimitMode pressure_boundary_mode;    ///< BC for pressure
  SheathLimitMode temperature_boundary_mode; ///< BC for temperature
};

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
struct SheathBoundarySimple : public SheathBoundaryBase<hermes::SheathKind::simple> {
  /// # Input options
  /// - <name>  e.g. "sheath_boundary_simple"
  ///   - lower_y                  Boundary on lower y?
  ///   - upper_y                  Boundary on upper y?
  ///   - gamma_e                  Electron sheath heat transmission coefficient
  ///   - gamma_i                  Ion sheath heat transmission coefficient
  ///   - sheath_ion_polytropic    Ion polytropic coefficient in Bohm sound speed.
  ///   Default 1.
  ///   - wall_potential           Voltage of the wall [Volts]
  ///   - secondary_electron_coef  Effective secondary electron emission coefficient
  ///   - sin_alpha                Sine of the angle between magnetic field line and wall
  ///   surface (0 to 1)
  ///   - always_set_phi           Always set phi field? Default is to only modify if
  ///   already set
  SheathBoundarySimple(const std::string& name, Options& options,
                       [[maybe_unused]] Solver* solver)
      : SheathBoundaryBase(options[name], options["units"]["eV"]) {}
};

/// Boundary condition at the wall in Y
///
/// This is a collective component, because it couples all charged species
///
/// These are based on
/// "Boundary conditions for the multi-ion magnetized plasma-wall transition"
///  by D.Tskhakaya, S.Kuhn. JNM 337-339 (2005), 405-409
///
/// Notes:
///   - The approximation used here is for ions having similar
///     gyro-orbit sizes
///   - No boundary condition is applied to neutral species
///   - Boundary conditions are applied to field-aligned fields
///     using to/fromFieldAligned
struct SheathBoundary : public SheathBoundaryBase<hermes::SheathKind::normal> {
  /// # Input options
  /// - <name>  e.g. "sheath_boundary_simple"
  ///   - lower_y                  Boundary on lower y?
  ///   - upper_y                  Boundary on upper y?
  ///   - gamma_e                  Electron sheath heat transmission coefficient
  ///   - wall_potential           Voltage of the wall [Volts]
  ///   - floor_potential          Apply floor to sheath potential?
  ///   - secondary_electron_coef  Effective secondary electron emission coefficient
  ///   - sin_alpha                Sine of the angle between magnetic field line and wall
  ///                              surface (0 to 1)
  ///   - always_set_phi           Always set phi field? Default is to only modify if
  ///                              already set
  SheathBoundary(const std::string& name, Options& options,
                 [[maybe_unused]] Solver* solver)
      : SheathBoundaryBase(options[name], options["units"]["eV"]) {}
};

/// Insulating sheath boundary condition at the wall in Y
///
/// This is a collective component, because it couples all charged species
///
/// Adapted from the `sheath_boundary` component, but always sets the current
/// density to zero
struct SheathBoundaryInsulating
    : public SheathBoundaryBase<hermes::SheathKind::insulating> {
  /// # Input options
  /// - <name>  e.g. "sheath_boundary_simple"
  ///   - lower_y                  Boundary on lower y?
  ///   - upper_y                  Boundary on upper y?
  ///   - gamma_e                  Electron sheath heat transmission coefficient
  ///   - gamma_i                  Ion sheath heat transmission coefficient
  ///   - sheath_ion_polytropic    Ion polytropic coefficient in Bohm sound speed.
  ///                              Default 1.
  ///   - wall_potential           Voltage of the wall [Volts]
  ///   - secondary_electron_coef  Effective secondary electron emission coefficient
  ///   - sin_alpha                Sine of the angle between magnetic field line and wall
  ///                              surface (0 to 1)
  SheathBoundaryInsulating(const std::string& name, Options& options,
                           [[maybe_unused]] Solver* solver)
      : SheathBoundaryBase(options[name], options["units"]["eV"]) {}
};

template class SheathBoundaryBase<hermes::SheathKind::normal>;
template class SheathBoundaryBase<hermes::SheathKind::simple>;
template class SheathBoundaryBase<hermes::SheathKind::insulating>;

namespace {
const RegisterComponent<SheathBoundarySimple>
    registercomponentsheathboundarysimple("sheath_boundary_simple");
const RegisterComponent<SheathBoundary>
    registercomponentsheathboundary("sheath_boundary");
const RegisterComponent<SheathBoundaryInsulating>
    registercomponentsheathboundaryinsulating("sheath_boundary_insulating");
} // namespace

#endif // SHEATH_BOUNDARY_SIMPLE_H
