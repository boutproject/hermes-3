#pragma once
#ifndef EVOLVE_PRESSURE_H
#define EVOLVE_PRESSURE_H

#include "../include/hermes_utils.hxx"
#include "component.hxx"
#include <bout/field3d.hxx>

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
  ///   - poloidal_flows       Include poloidal ExB flows? Default is true
  ///   - precon               Enable preconditioner? Note: solver may not use it even if
  ///                          enabled.
  ///   - p_div_v              Use p * Div(v) form? Default is v * Grad(p) form
  ///   - thermal_conduction   Include parallel heat conduction? Default is true
  ///
  /// - P<name>  e.g. "Pe", "Pd+"
  ///   - source     Source of pressure [Pa / s].
  ///                NOTE: This overrides mesh input P<name>_src
  ///   - source_only_in_core         Zero the source outside the closed field-line
  ///                                 region?
  ///   - neumann_boundary_average_z  Apply Neumann boundaries with Z average?
  ///
  EvolvePressure(std::string name, Options& options, Solver* solver);

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

  Field3D P;    ///< Pressure (normalised)
  Field3D T, N; ///< Temperature, density

  bool bndry_flux;
  bool neumann_boundary_average_z; ///< Apply neumann boundary with Z average?
  bool poloidal_flows;
  bool thermal_conduction; ///< Include thermal conduction?

  bool p_div_v; ///< Use p*Div(v) form? False -> v * Grad(p)

  bool evolve_log; ///< Evolve logarithm of P?
  Field3D logP;    ///< Natural logarithm of P

  BoutReal density_floor;     ///< Minimum density for calculating T
  bool low_n_diffuse_perp;    ///< Cross-field diffusion at low density?
  BoutReal temperature_floor; ///< Low temperature scale for low_T_diffuse_perp
  bool low_T_diffuse_perp;    ///< Add cross-field diffusion at low temperature?
  BoutReal pressure_floor;    ///< When non-zero pressure is needed
  bool low_p_diffuse_perp;    ///< Add artificial cross-field diffusion at low electron
                              ///< pressure?
  bool damp_p_nt;             ///< Damp P - N*T. Active when P < 0 or N < density_floor

  Field3D source, final_source; ///< External pressure source
  Field3D Sp;                   ///< Total pressure source
  FieldGeneratorPtr source_prefactor_function;

  BoutReal hyper_z;   ///< Hyper-diffusion
  BoutReal hyper_z_T; ///< 4th-order dissipation in T

  bool diagnose;                 ///< Output additional diagnostics?
  bool enable_precon;            ///< Enable preconditioner?
  BoutReal source_normalisation; ///< Normalisation factor [Pa/s]
  BoutReal time_normalisation;   ///< Normalisation factor [s]
  bool source_time_dependent;    ///< Is the input source time dependent?
  Field3D flow_xlow, flow_ylow;  ///< Energy flow diagnostics
  Field3D flow_ylow_advection;   ///< Advection energy flow diagnostics
  Field3D
      flow_ylow_viscous_heating; ///< Flow of kinetic energy due to numerical viscosity

  bool numerical_viscous_heating;  ///< Include heating due to numerical viscosity?
  bool fix_momentum_boundary_flux; ///< Fix momentum flux to boundary condition?
  Field3D Sp_nvh;                  ///< Pressure source due to artificial viscosity
  Field3D E_PdivV,
      E_VgradP; ///< Diagnostic energy source terms for p*Div(V) and V*Grad(P)

  ////// Diagnostics for pressure equation terms

  bool diagnose_terms;

  // Physics terms
  Field3D ddtP_ExB;
  Field3D ddtP_advection;
  Field3D ddtP_PdivV;
  Field3D ddtP_VgradP;
  Field3D ddtP_flutter1, ddtP_flutter2, ddtP_flutter3;
  Field3D ddtP_viscousheat;
  Field3D ddtP_cond;

  // Diffusion in low density/pressure/temperature regions
  Field3D ddtP_low_N_diff;
  Field3D ddtP_low_N_diff_perp;
  Field3D ddtP_low_T_diff_perp;
  Field3D ddtP_low_P_diff_perp;

  // Hyperdiffusion
  Field3D ddtP_hyperz;
  Field3D ddtP_T_hyperz;

  // Term to enforce P = N*T
  Field3D ddtP_damp_NT;

  // User source and source from elsewhere in code (e.g. reactions)
  Field3D ddtP_user_source;
  Field3D ddtP_component_source;

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
  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<EvolvePressure> registercomponentevolvepressure("evolve_pressure");
}

#endif // EVOLVE_PRESSURE_H
