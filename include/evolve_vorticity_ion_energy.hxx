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

namespace bout::ion_energy {
struct PerpendicularFaceFluxes {
  BoutReal xp, xm, zp, zm;
};
/// Integrated fluxes on the four faces of a physical cell, including g13.
PerpendicularFaceFluxes perpendicularFaceFluxes(const Field3D& a, const Field3D& f,
                                                const Ind3D& i, bool bndry_flux = true);
/// Conservative X-Z operator, including the cross metric g13.
Field3D Div_a_Grad_perp_XZ(const Field3D& a, const Field3D& f, bool bndry_flux = true);
/// Face-based energy density consistent with Div_a_Grad_perp_XZ.
Field3D kineticEnergyPerp(const Field3D& a, const Field3D& f);
} // namespace bout::ion_energy

/// Evolves vorticity and ion total energy
///
/// These are evolved together because the ion perpendicular kinetic energy
/// depends on the ion pressure through the diamagnetic drift.
/// This formulation means that dissipation in the vorticity equation
/// is automatically converted to ion heating.
/// All species energy_source fields supply thermal energy. Electron pressure
/// must be evolved by evolve_pressure; electron parallel electrical work is
/// accounted for by the electron momentum equation.
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
  ///   - viscosity            Kinematic viscosity [m^2/s], default 0.0
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
  ///     - energy_source. Thermal energy, excluding parallel kinetic energy
  /// - fields
  ///   - phi      Electrostatic potential -> ExB drift
  ///
  void finally(const Options& state) override;

  // Save and restore potential phi
  void restartVars(Options& state) override;

  void outputVars(Options& state) override;

  static constexpr auto type = "evolve_vorticity_ion_energy";

private:
  bool restore_phi_on_restart{true}; ///< Read restart phi only on the first hook call
  Field3D Vort;                      ///< Solver vorticity; radial boundaries store phi
  Field3D vorticity{0.0};            ///< Physical vorticity, including ordinary guards
  Field3D Eion;                      ///< Evolving ion energy

  Field3D ion_KE_parallel; ///< Parallel ion kinetic energy K_||i
  Field3D ion_KE_perp;     ///< Perpendicular ion kinetic energy K⊥

  BoutReal average_atomic_mass;    ///< Weighted average atomic mass
  bool poloidal_flows;             ///< Include poloidal (Y) flows?
  bool neumann_boundary_average_z; ///< Set pressure constant in Z, Neumann in X
  BoutReal density_floor;          ///< Used to calculate Ti

  bool include_viscosity;                ///< Was kinematic viscosity specified?
  Coordinates::FieldMetric viscosity;    ///< Normalised kinematic viscosity
  Field3D viscous_heating{0.0};          ///< Face-based -nu grad(Vort) dot grad(Phi)
  Field3D electron_viscous_heating{0.0}; ///< Electron mass fraction of viscous heating

  Field3D phi{0.0}; ///< Electrostatic potential ϕ
  Field3D chi{0.0}; ///< Auxiliary current potential, L(chi) = physical ∂tVort
  Field3D Pi_sum;   ///< Total ion pressure
  Field3D psi;      ///< Pressure potential Σ A_s p_s / (Abar Z_s)
  Field3D Ti;       ///< Common ion temperature
  Coordinates::FieldMetric boussinesq_coeff; ///< c = mass / B^2
  std::unique_ptr<Laplacian> phiSolver;      ///< Laplacian solver in X-Z

  bool phi_boundary_relax;         ///< Relax boundary to zero-gradient
  BoutReal phi_boundary_timescale; ///< Relaxation timescale [normalised]
  bool phi_core_averagey;          ///< Average phi core boundary in Y?

  int phi_max_iterations; ///< Maximum number of Vort,E -> phi, Pi iterations
  BoutReal phi_atol;      ///< Absolute tolerance on phi
  BoutReal phi_rtol;      ///< Relative tolerance on phi
  BoutReal vort_atol;     ///< Absolute tolerance on forward vorticity residual

  // Internal transport sources are excluded from conversion of external
  // sources to total energy: their energy flux is already in DivEpol.
  Options polarisation_sources;

  VectorMetric Curlb_B; // Curvature vector Curl(b/B)

  // Divergence of currents
  Field3D DivJdia;      ///< Divergence of diamagnetic current
  Field3D DivJReynolds; ///< Divergence of Reynolds stress current
  Field3D DivJpol;      ///< Polarisation transport divergence with closed physical faces
  Field3D DivJviscous{0.0}; ///< Divergence of nu grad_perp(Vort)

  // Divergence of energy fluxes
  Field3D DivEExB;      ///< Divergence of ExB energy flux
  Field3D DivEReynolds; ///< Divergence of Reynolds stress energy flux
  Field3D DivEelectric; ///< Divergence of ion electrical work Div(phi (Jpar + Jdia))
  Field3D DivEparallel; ///< Divergence of parallel ion enthalpy and kinetic energy flux
  Field3D DivEpol;      ///< Divergence of polarisation drift energy flux
  Field3D DivEviscous{0.0}; ///< Divergence of Phi nu grad_perp(Vort)

  Field3D phi_Div_Je;        ///< Energy transferred to/from electrons
  Field3D electron_exb_work; ///< Thermal work -Jdia_e dot Grad(phi)

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

  /// Initialise encoded boundary midpoints from phi, including old restarts.
  void initialiseBoundaryState();
  /// Reimpose phi midpoints from the solver's boundary state during recovery.
  void applyPotentialBoundary();
  /// Relax encoded midpoints towards averages of the current recovered phi.
  void boundaryRHS();

  /// Shared viscous current, energy flux and heating, using physical vorticity guards.
  void calculateViscosity();

  /// Apply a boundary condition to given pressure field `P`
  void applyPressureBoundary(Field3D& P) const;

  /// Sum thermal sources and kinetic energy changes from momentum and particle sources
  Field3D sumIonEnergySources(const Options& allspecies) const;

  /// Use the Laplacian as a preconditioner for a finite-volume residual step.
  Field3D correctPotential(const Field3D& rhs, const Field3D& potential) const;

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
