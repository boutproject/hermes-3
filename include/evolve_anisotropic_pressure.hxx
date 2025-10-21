#pragma once
#ifndef EVOLVE_ANISOTROPIC_PRESSURE_H
#define EVOLVE_ANISOTROPIC_PRESSURE_H

#include "../include/hermes_utils.hxx"
#include "component.hxx"
#include <bout/field3d.hxx>

/// Evolves species anisotropic pressure (P_parallel & P_perp) in time.
/// This formulation evolves
///   E    Energy, combining thermal and parallel kinetic
///   PA   Pressure anisotropy, PA = P_par - P_perp
///
/// # Mesh inputs
///
/// P<name>_src   A source of pressure, in Pascals per second
///               This can be over-ridden by the `source` option setting.
///
struct EvolveAnisotropicPressure : public Component {
  ///
  /// # Inputs
  ///
  /// - <name>
  ///   - bndry_flux           Allow flows through radial boundaries? Default is true
  ///   - density_floor        Minimum density floor. Default 1e-5 normalised units.
  ///   - diagnose             Output additional diagnostic fields?
  ///   - poloidal_flows       Include poloidal ExB flows? Default is true
  ///   - precon               Enable preconditioner? Note: solver may not use it even if
  ///   enabled.
  ///   - p_div_v              Use p * Div(v) form? Default is v * Grad(p) form
  ///   - thermal_conduction   Include parallel heat conduction? Default is true
  ///
  /// - P<name>  e.g. "Pe", "Pd+"
  ///   - source     Source of pressure [Pa / s].
  ///                NOTE: This overrides mesh input P<name>_src
  ///   - source_only_in_core         Zero the source outside the closed field-line
  ///   region?
  ///   - neumann_boundary_average_z  Apply Neumann boundaries with Z average?
  ///
  EvolveAnisotropicPressure(std::string name, Options& options, Solver* solver);

  /// Inputs
  /// - species
  ///   - <name>
  ///     - density
  ///
  /// Sets
  /// - species
  ///   - <name>
  ///     - pressure       <- This is (2pressure_perp + pressure_par)/3
  ///     - temperature   Requires density
  ///     - pressure_par
  ///     - pressure_perp
  ///
  void transform(Options& state) override;

  ///
  /// Optional inputs
  ///
  /// - species
  ///   - <name>
  ///     - velocity. Must have sound_speed or temperature
  ///     - energy_source
  ///     - collision_rate  (needed if thermal_conduction on)
  /// - fields
  ///   - phi      Electrostatic potential -> ExB drift
  ///
  void finally(const Options& state) override;

  void outputVars(Options& state) override;

private:
  std::string name; ///< Short name of the species e.g. h+

  BoutReal adiabatic_index; ///< Ratio of specific heats, γ = Cp / Cv
  BoutReal Cv; ///< /// Heat capacity at constant volume (3/2 for ideal monatomic gas)

  Field3D E;  ///< Energy (normalised): P + 1/2 m n v^2
  Field3D PA; ///< Pressure anisotropy: P_par - P_perp

  Field3D P;             ///< Pressure (normalised): P = (2 P_perp + P_par) / 3
  Field3D P_perp, P_par; ///< Perpendicular and parallel pressure (normalised)

  Field3D T, N; ///< Temperature, density

  bool bndry_flux;
  bool neumann_boundary_average_z; ///< Apply neumann boundary with Z average?
  bool poloidal_flows;

  BoutReal density_floor;     ///< Minimum density for calculating T
  BoutReal temperature_floor; ///< Low temperature scale for low_T_diffuse_perp
  BoutReal pressure_floor;    ///< When non-zero pressure is needed

  Field3D source, final_source; ///< External pressure source
  Field3D Sp;                   ///< Total pressure source
  FieldGeneratorPtr source_prefactor_function;

  bool diagnose;                 ///< Output additional diagnostics?
  BoutReal source_normalisation; ///< Normalisation factor [Pa/s]
  BoutReal time_normalisation;   ///< Normalisation factor [s]
  bool source_time_dependent;    ///< Is the input source time dependent?

  bool fix_momentum_boundary_flux; ///< Fix momentum flux to boundary condition?
};

namespace {
RegisterComponent<EvolveAnisotropicPressure>
    registercomponentevolveanisotropicpressure("evolve_anisotropic_pressure");
}

#endif // EVOLVE_ANISOTROPIC_PRESSURE_H
