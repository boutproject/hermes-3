#pragma once
#ifndef EVOLVE_PRESSURE_H
#define EVOLVE_PRESSURE_H

#include <field3d.hxx>

#include "component.hxx"

/// Evolves species pressure in time
///
/// # Mesh inputs
///
/// P<name>_src   A source of pressure, in Pascals per second
///               This can be over-ridden by the `source` option setting.
///
struct EvolvePressure : public Component {
  ///
  /// # Inputs
  ///
  /// - <name>
  ///   - bndry_flux           Allow flows through radial boundaries? Default is true
  ///   - density_floor        Minimum density floor. Default 1e-5 normalised units.
  ///   - diagnose             Output additional diagnostic fields?
  ///   - evolve_log           Evolve logarithm of pressure? Default is false
  ///   - hyper_z              Hyper-diffusion in Z
  ///   - kappa_coefficient    Heat conduction constant. Default is 3.16 for electrons, 3.9 otherwise
  ///   - kappa_limit_alpha    Flux limiter, off by default.
  ///   - poloidal_flows       Include poloidal ExB flows? Default is true
  ///   - precon               Enable preconditioner? Note: solver may not use it even if enabled.
  ///   - p_div_v              Use p * Div(v) form? Default is v * Grad(p) form
  ///   - thermal_conduction   Include parallel heat conduction? Default is true
  ///
  /// - P<name>  e.g. "Pe", "Pd+"
  ///   - source     Source of pressure [Pa / s].
  ///                NOTE: This overrides mesh input P<name>_src
  ///
  EvolvePressure(std::string name, Options& options, Solver* solver);

  /// Inputs
  /// - species
  ///   - <name>
  ///     - density
  ///
  /// Sets
  /// - species
  ///   - <name>
  ///     - pressure
  ///     - temperature   Requires density
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
  void precon(const Options &UNUSED(state), BoutReal gamma) override;
private:
  std::string name; ///< Short name of the species e.g. h+

  Field3D P;    ///< Pressure (normalised)
  Field3D T, N; ///< Temperature, density

  bool bndry_flux;
  bool poloidal_flows;
  bool thermal_conduction;    ///< Include thermal conduction?
  BoutReal kappa_coefficient; ///< Leading numerical coefficient in parallel heat flux calculation
  BoutReal kappa_limit_alpha; ///< Flux limit if >0

  bool p_div_v; ///< Use p*Div(v) form? False -> v * Grad(p)

  bool evolve_log; ///< Evolve logarithm of P?
  Field3D logP;    ///< Natural logarithm of P

  BoutReal density_floor; ///< Minimum density for calculating T
  Field3D kappa_par; ///< Parallel heat conduction coefficient

  Field3D source; ///< External pressure source
  Field3D Sp;     ///< Total pressure source

  BoutReal hyper_z; ///< Hyper-diffusion

  bool diagnose; ///< Output additional diagnostics?
  bool enable_precon; ///< Enable preconditioner?
};

namespace {
RegisterComponent<EvolvePressure> registercomponentevolvepressure("evolve_pressure");
}

#endif // EVOLVE_PRESSURE_H
