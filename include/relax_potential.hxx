#pragma once
#ifndef RELAX_POTENTIAL_H
#define RELAX_POTENTIAL_H

#include "../include/div_ops.hxx"
#include "../include/guarded_options.hxx"

#include <bout/bout_types.hxx>
#include <bout/coordinates.hxx>
#include <bout/options.hxx>
#include <bout/vector2d.hxx>
#include <bout/vectormetric.hxx>

#include "component.hxx"

#include <memory>
#include <string>

/// Evolve vorticity and potential in time.
///
/// Uses a relaxation method for the potential, which is valid for
/// steady state, but not for timescales shorter than the relaxation
/// timescale.
///
struct RelaxPotential : public NamedComponent<RelaxPotential> {
  /// Options
  ///
  /// - <name>
  ///   - boussinesq: bool, default true
  ///     Use the Boussinesq approximation?
  ///   - average_atomic_mass: float, default 2.0
  ///     Weighted average ion atomic mass for polarisation current
  ///   - bndry_flux: bool, default true
  ///     Allow flows through radial (X) boundaries?
  ///   - collisional_friction: bool, default false
  ///     Damp vorticity based on mass-weighted collision frequency?
  ///   - diagnose: bool, false
  ///     Output additional diagnostics?
  ///   - diamagnetic: bool, default true
  ///     Include diamagnetic current, using curvature vector?
  ///   - diamagnetic_polarisation: bool, default true
  ///     Include ion diamagnetic drift in polarisation current?
  ///   - exb_advection: bool, default true
  ///     Include ExB advection (nonlinear term)?
  ///   - exb_advection_simplified: bool, default true
  ///     Use the simplified ExB advection form?
  ///   - evolve_vorticity: bool, default true
  ///     Evolve vorticity, rather than deriving it from phi?
  ///   - hyper_z: float, default -1.0
  ///     Hyper-viscosity in Z. Negative values disable it.
  ///   - lambda_1: float, default 3e8 [SI]
  ///     Relaxation rate for the phi1 evolution equation
  ///   - lambda_2: float, default 1.0
  ///     Scaling between phi1 and phi, with phi1 = lambda_2 * phi
  ///   - phi_boundary_relax: bool, default false
  ///     Relax radial phi boundaries towards zero-gradient?
  ///   - phi_boundary_timescale: float, 1e-4
  ///     Timescale for phi boundary relaxation [seconds]
  ///   - phi_core_averagey: bool, default false
  ///     Average phi core boundary in Y? (if phi_boundary_relax)
  ///   - phi_dissipation: bool, default true
  ///     Parallel dissipation of potential (Recommended)
  ///   - phi_sheath_dissipation: bool, default false
  ///     Add dissipation when phi < 0 at the sheath?
  ///   - poloidal_flows: bool, default true
  ///     Include poloidal ExB flow?
  ///   - sheath_boundary: bool, default false
  ///     If phi_boundary_relax is false, set the radial boundary to the sheath potential?
  ///   - viscosity_perp: Field2D, default 0.0
  ///     Kinematic viscosity in perpendicular diraction  [m^2/s]
  ///   - viscosity_par: Field2D, default 0.0
  ///     Kinematic viscosity in parallel diraction  [m^2/s]
  ///   - vort_dissipation: bool, default false
  ///     Parallel dissipation of vorticity?
  ///   - damp_core_vorticity: bool, default false
  ///     Damp axisymmetric component of vorticity in cell next to core boundary
  ///

  RelaxPotential(std::string name, Options& options, Solver* solver);

  /// Optional inputs
  /// - fields
  ///   - DivJextra    Divergence of current, including parallel current
  ///                  Not including diamagnetic or polarisation currents
  ///
  void finally(const Options& state) override;

  void outputVars(Options& state) override;

  // The following are public functions for unit testing

  /// Diamagnetic term in vorticity, weighted by mass/charge
  Field3D calculatePihat(GuardedOptions allspecies);

  /// Apply the configured radial boundary condition to phi
  void applyPhiBoundary(Field3DParallel& phi, GuardedOptions state);

  /// Calculates Div(Jdia) and sets energy_source for all
  /// charged species with pressure.
  /// The parallel boundary of phi is set using extrapolation
  Field3D calculateDivJdia(Field3DParallel& phi, GuardedOptions allspecies);

  /// Calculates the collisional friction current divergence.
  Field3D calculateDivJcol(const Field3DParallel& phi, const Field3D& pi_hat,
                           GuardedOptions allspecies);

  /// Calculate the ExB advection contribution to ddt(Vort).
  Field3D calculateExBAdvectionSource(const Field3D& vort, const Field3DParallel& phi,
                                      const Field3D& pi_hat);

  /// Calculate the current-divergence source terms added in finally().
  Field3D calculateParallelCurrentSource(const Options& state);

  /// Calculate viscosity and dissipation contributions added in finally().
  Field3D calculateDissipationSource(const Options& state, const Field3D& vort,
                                     const Field3DParallel& phi);

  /// Calculate ddt(phi1) from the chosen evolution mode.
  Field3D calculatePhi1Source(const Field3D& vort, const Field3D& vort_rhs,
                              const Field3D& vort_from_phi) const;

  /// Calculate vorticity from potential and species
  Field3D vorticity(const Field3DParallel& phi, GuardedOptions& allspecies);

  static constexpr auto type = "relax_potential";

private:
  bool evolve_vorticity; ///< Evolve vorticity?
  Field3D Vort;          // Evolving vorticity
  Field3D Vort_from_phi; ///< Vort calculated from phi

  Field3D phi1;        // Scaled electrostatic potential, evolving in time ϕ_1 = λ_2 ϕ
  Field3DParallel phi; // Electrostatic potential

  Field3D Pi_hat; ///< Contribution from ion pressure, weighted by atomic mass / charge

  bool exb_advection;            ///< Include nonlinear ExB advection?
  BoutReal scale_ExB;            ///< Multiply ExB advection by this factor
  bool exb_advection_simplified; ///< Simplify nonlinear ExB advection form?
  bool diamagnetic;              ///< Include diamagnetic current?
  bool diamagnetic_polarisation; ///< Include diamagnetic drift in polarisation current?
  bool diamagnetic_bracketform;  ///< Use Arakawa form of diamagnetic current?
  bool boussinesq;               ///< Use the Boussinesq approximation?
  BoutReal average_atomic_mass; ///< Weighted average atomic mass, for polarisaion current
                                // (Boussinesq approximation)
  BoutReal permittivity;        ///< Vacuum permittivity
  bool poloidal_flows;          ///< Include poloidal ExB flow?
  bool bndry_flux;              ///< Allow flows through radial boundaries?

  bool collisional_friction; ///< Damping of vorticity due to collisional friction

  bool sheath_boundary; ///< Set outer boundary to j=0?

  bool vort_dissipation;       ///< Parallel dissipation of vorticity
  bool phi_dissipation;        ///< Parallel dissipation of potential
  bool phi_sheath_dissipation; ///< Dissipation at the sheath if phi < 0
  bool damp_core_vorticity;    ///< Damp axisymmetric component of vorticity

  bool phi_boundary_relax;           ///< Relax boundary to zero-gradient
  BoutReal phi_boundary_timescale;   ///< Relaxation timescale [normalised]
  BoutReal phi_boundary_last_update; ///< Time when last updated
  bool phi_core_averagey;            ///< Average phi core boundary in Y?

  bout::FieldMetricParallel Bsq; ///< SQ(coord->Bxy())
  VectorMetric Curlb_B;          ///< Curvature vector Curl(b/B)
  bout::FieldMetricParallel logB;
  BoutReal hyper_z;      ///< Hyper-viscosity in Z
  Field3D viscosity;     ///< Perpendicular Kinematic viscosity
  Field3D viscosity_par; ///< Parallel Kinematic viscosity
  bool include_viscosity, include_viscosity_par;
  // Relax-potential related variables
  BoutReal lambda_1; ///< Relaxation parameters.  NOTE: lambda_1 has dimensions!
  BoutReal lambda_2; ///< Relaxation parameters

  // Diagnostic outputs
  Field3D DivJdia, DivJcol; // Divergence of diamagnetic and collisional current
  Field3D Vort_rhs;         // Rate of change of vorticity, even if not evolved

  bool diagnose; ///< Output additional diagnostics?

  std::shared_ptr<FCI::dagp_fv> dagp_op; ///< Perp. diffusion operator for Fci
  Coordinates::FieldMetric
      bracket_factor; ///< For non-Clebsch coordinate systems (e.g. FCI)

  /// Optional inputs
  ///
  /// - species
  ///   - pressure and charge => Calculates diamagnetic terms [if diamagnetic=true]
  ///   - pressure, charge and mass => Calculates polarisation current terms
  ///     [if diamagnetic_polarisation=true]
  ///
  /// Sets in the state
  /// - species
  ///   - [if has pressure and charge]
  ///     - energy_source
  /// - fields
  ///   - vorticity
  ///   - phi         Electrostatic potential
  ///   - DivJdia     Divergence of diamagnetic current [if diamagnetic=true]
  ///
  /// Note: Diamagnetic current calculated here, but could be moved
  ///       to a component with the diamagnetic drift advection terms
  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<RelaxPotential> registercomponentrelaxpotential;
}

#endif // RELAX_POTENTIAL_H
