#pragma once
#ifndef EVOLVE_VORTICITY_ION_ENERGY_H
#define EVOLVE_VORTICITY_ION_ENERGY_H

#include <bout/bout_types.hxx>
#include <bout/coordinates.hxx>
#include <bout/field3d.hxx>
#include <bout/invert_laplace.hxx>
#include <bout/vectormetric.hxx>

#include "component.hxx"
#include "guarded_options.hxx"

#include <memory>
#include <string>

/// Evolves vorticity and ion total energy
///
/// These are evolved together because the ion perpendicular kinetic energy
/// depends on the ion pressure through the diamagnetic drift.
/// This formulation means that dissipation in the vorticity equation
/// is automatically converted to ion heating.
///
/// # Mesh inputs
///
/// P<name>_src   A source of pressure, in Pascals per second
///               This can be over-ridden by the `source` option setting.
///
struct EvolveVorticityIonEnergy : public NamedComponent<EvolveVorticityIonEnergy> {
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
  ///   enabled.
  ///
  /// - E<name>  e.g. "Ee", "Ed+"
  ///   - source     Source of energy [W / s].
  ///                NOTE: This overrides mesh input P<name>_src
  ///   - source_only_in_core         Zero the source outside the closed field-line
  ///                                 region?
  ///   - neumann_boundary_average_z  Apply Neumann boundaries with Z average?
  ///
  EvolveVorticityIonEnergy(std::string name, Options& options, Solver* solver);

  ///
  /// Optional inputs
  ///
  /// - species
  ///   - <name>
  ///     - velocity. Must have sound_speed or temperature
  ///     - energy_source
  /// - fields
  ///   - phi      Electrostatic potential -> ExB drift
  ///
  void finally(const Options& state) override;

  // Save and restore potential phi
  void restartVars(Options& state) override {

    // NOTE: This is a hack because we know that the loaded restart file
    //       is passed into restartVars in PhysicsModel::postInit
    // The restart value should be used in init() rather than here
    static bool first = true;
    if (first and state.isSet("phi")) {
      first = false;
      phi = state["phi"].as<Field3D>();
    }

    // Save the potential
    set_with_attrs(state["phi"], phi,
                   {{"long_name", "plasma potential"}, {"source", "vorticity"}});
  }

  void outputVars(Options& state) override;

  static constexpr auto type = "evolve_vorticity_ion_energy";

private:
  Field3D Vort; ///< Evolving vorticity
  Field3D Eion; ///< Evolving ion energy

  Field3D ion_KE_parallel; ///< Parallel ion kinetic energy K_||i
  Field3D ion_KE_perp;     ///< Perpendicular ion kinetic energy K⊥

  BoutReal average_atomic_mass;    ///< Weighted average atomic mass
  bool poloidal_flows;             ///< Include poloidal (Y) flows?
  bool neumann_boundary_average_z; ///< Set pressure constant in Z, Neumann in X
  BoutReal density_floor;          ///< Used to calculate Ti

  Field3D phi{0.0};    ///< Electrostatic potential ϕ
  Field3D dPhidt{0.0}; ///< Generalised potential time derivative ∂tΦ
  Field3D Pi_sum;      ///< Total ion pressure
  Field3D psi;         ///< Pressure potential Σ A_s p_s / (Abar Z_s)
  Field3D Ti;          ///< Common ion temperature
  Coordinates::FieldMetric boussinesq_coeff; ///< c = mass / B^2
  std::unique_ptr<Laplacian> phiSolver;      ///< Laplacian solver in X-Z

  bool phi_boundary_relax;           ///< Relax boundary to zero-gradient
  BoutReal phi_boundary_timescale;   ///< Relaxation timescale [normalised]
  BoutReal phi_boundary_last_update; ///< Time when last updated
  bool phi_core_averagey;            ///< Average phi core boundary in Y?

  int phi_max_iterations; ///< Maximum number of Vort,E -> phi, Pi iterations
  BoutReal phi_atol;      ///< Absolute tolerance on phi
  BoutReal phi_rtol;      ///< Relative tolerance on phi

  VectorMetric Curlb_B; // Curvature vector Curl(b/B)

  // Divergence of currents
  Field3D DivJdia;      ///< Divergence of diamagnetic current
  Field3D DivJReynolds; ///< Divergence of Reynolds stress current

  // Divergence of energy fluxes
  Field3D DivEExB;      ///< Divergence of ExB energy flux
  Field3D DivEReynolds; ///< Divergence of Reynolds stress energy flux
  Field3D DivEelectric; ///< Divergence of ion electrical work Div(phi (Jpar + Jdia))
  Field3D DivEparallel; ///< Divergence of parallel ion enthalpy and kinetic energy flux
  Field3D DivEpol;      ///< Divergence of polarisation drift energy flux

  Field3D phi_Div_Je; ///< Energy transferred to/from electrons

  Field3D Eion_src; ///< Ion energy sources from other components

  bool diagnose; ///< Output additional diagnostics?

  /// Inputs
  /// - species
  ///   - <name>
  ///     - AA
  ///     - density
  ///     - velocity
  ///
  /// Sets
  /// - species
  ///   - <name>
  ///     - pressure
  ///     - temperature
  ///
  void transform_impl(GuardedOptions& state) override;

  /// Apply a boundary condition to given pressure field `P`
  void applyPressureBoundary(Field3D& P) const;

  /// Sum thermal sources and kinetic energy changes from momentum and particle sources
  Field3D sumIonEnergySources(const Options& allspecies) const;

  /// Calculate ion perpendicular kinetic energy, K⊥ = 1/2 c |∇⊥Φ|^2
  ///
  /// This should be consistent with the vorticity discretisation,
  /// so that the integration by parts of Vort * Phi over the domain
  /// gives the kinetic perpendicular energy plus boundary terms.
  ///
  Field3D calculateKEperp(const Field3D& phi, const Field3D& psi) const;
};

namespace {
RegisterComponent<EvolveVorticityIonEnergy> registercomponentevolvevorticityionenergy;
}

#endif // EVOLVE_VORTICITY_ION_ENERGY_H
