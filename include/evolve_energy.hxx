#pragma once
#ifndef EVOLVE_ENERGY_H
#define EVOLVE_ENERGY_H

#include <bout/field3d.hxx>

#include "component.hxx"

/// Evolves species internal energy in time
///
/// # Mesh inputs
///
/// P<name>_src   A source of pressure, in Pascals per second
///               This can be over-ridden by the `source` option setting.
///
struct EvolveEnergy : public Component {
  ///
  /// # Inputs
  ///
  /// - <name>
  ///   - bndry_flux           Allow flows through radial boundaries? Default is true
  ///   - density_floor        Minimum density floor. Default 1e-5 normalised units.
  ///   - diagnose             Output additional diagnostic fields?
  ///   - evolve_log           Evolve logarithm of pressure? Default is false
  ///   - hyper_z              Hyper-diffusion in Z
  ///   - kappa_coefficient    Heat conduction constant. Default is 3.16 for
  ///   electrons, 3.9 otherwise
  ///   - kappa_limit_alpha    Flux limiter, off by default.
  ///   - poloidal_flows       Include poloidal ExB flows? Default is true
  ///   - precon               Enable preconditioner? Note: solver may not use it even if
  ///   enabled.
  ///   - thermal_conduction   Include parallel heat conduction? Default is true
  ///
  /// - E<name>  e.g. "Ee", "Ed+"
  ///   - source     Source of energy [W / s].
  ///                NOTE: This overrides mesh input P<name>_src
  ///   - source_only_in_core         Zero the source outside the closed field-line
  ///                                 region?
  ///   - neumann_boundary_average_z  Apply Neumann boundaries with Z average?
  ///
  EvolveEnergy(std::string name, Options& options, Solver* solver);

  /// Inputs
  /// - species
  ///   - <name>
  ///     - density
  ///     - velocity
  ///
  /// Sets
  /// - species
  ///   - <name>
  ///     - pressure
  ///     - temperature
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

  /// Preconditioner
  ///
  void precon(const Options& UNUSED(state), BoutReal gamma) override;

private:
  std::string name; ///< Short name of the species e.g. h+

  Field3D E;    ///< Energy (normalised): P + 1/2 m n v^2
  Field3D P;    ///< Pressure (normalised)
  Field3D T, N; ///< Temperature, density

  BoutReal adiabatic_index; ///< Ratio of specific heats, γ = Cp / Cv
  BoutReal Cv; /// Heat capacity at constant volume (3/2 for ideal monatomic gas)
  bool bndry_flux;
  bool neumann_boundary_average_z; ///< Apply neumann boundary with Z average?
  bool poloidal_flows;
  bool thermal_conduction;    ///< Include thermal conduction?
  std::vector<std::string> collision_names; ///< Collisions used for collisionality
  std::string conduction_collisions_mode;  ///< Collision selection, either legacy or braginskii
  Field3D nu;   ///< Collision frequency for conduction
  BoutReal kappa_coefficient; ///< Leading numerical coefficient in parallel heat flux
                              ///< calculation
  BoutReal kappa_limit_alpha; ///< Flux limit if >0

  bool evolve_log; ///< Evolve logarithm of E?
  Field3D logE;    ///< Natural logarithm of E

  BoutReal density_floor; ///< Minimum density for calculating T
  Field3D kappa_par;      ///< Parallel heat conduction coefficient

  Field3D source; ///< External power source
  Field3D Se;     ///< Total energy source

  BoutReal hyper_z; ///< Hyper-diffusion

  bool diagnose;      ///< Output additional diagnostics?
  bool enable_precon; ///< Enable preconditioner?
  Field3D flow_xlow, flow_ylow; ///< Energy flow diagnostics
};

namespace {
RegisterComponent<EvolveEnergy> registercomponentevolveenergy("evolve_energy");
}

#endif // EVOLVE_ENERGY_H
