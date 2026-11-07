#pragma once
#ifndef EVOLVE_VORTICITY_ION_ENERGY_H
#define EVOLVE_VORTICITY_ION_ENERGY_H

#include <bout/bout_enum_class.hxx>
#include <bout/bout_types.hxx>
#include <bout/coordinates.hxx>
#include <bout/field3d.hxx>
#include <bout/invert_laplace.hxx>
#include <bout/vectormetric.hxx>

#include "component.hxx"
#include "guarded_options.hxx"

#include <memory>
#include <string>
#include <vector>

BOUT_ENUM_CLASS_NS(bout::ion_energy, RadialBoundary, fixed, free_slip, no_slip);
BOUT_ENUM_CLASS_NS(bout::ion_energy, PolarisationACBoundary, neumann, dirichlet);

namespace bout::ion_energy {
/// All radial settings are derived once from the two boundary enums.
struct RadialBoundaryModel {
  bool walls{false};
  RadialBoundary mode{RadialBoundary::fixed};
  PolarisationACBoundary chi_ac{PolarisationACBoundary::neumann};
  bool evolve_wall_potential{false};
  bool vorticity_wall_zero{false};
  bool viscous_wall_flux{false};
  bool chi_wall_flux{false};
  bool core_flux_surface_average{false};
  BoutReal phi_boundary_velocity{0.0}; ///< Normalised velocity, only for no_slip
  int phi_flags{INVERT_SET};
  int chi_flags{INVERT_AC_GRAD};
};

/// Mean radial Neumann and zero AC midpoint, copied to all radial guards.
/// Physical Y rows and X walls only; interior and periodic X are preserved.
void applyWallPressureBoundary(Field3D& f);
/// Neumann vorticity, or zero midpoint for a stress-free wall; all X guards.
void applyVorticityWallBoundary(Field3D& w, bool midpoint_zero);

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
/// Estimate the inward wall derivative from an encoded arithmetic midpoint.
///
/// Reflection sets phi_guard = 2 * phi_midpoint - phi0, so this is the centred
/// face derivative (phi0 - phi_guard) / dx0. With the guard centre mirrored
/// across the face, it is second order and exact for quadratic profiles.
///
/// @param phi_midpoint Arithmetic mean of guard and interior potential, normally
///                     averaged over Z; this is the value stored in Vort guards.
/// @param phi0 Potential in the first interior cell, normally averaged over Z.
/// @param dx0 Positive x-coordinate spacing of the first interior cell.
/// @returns Potential derivative per unit x-coordinate distance, directed into
///          the domain at either wall. The physical normal metric is not applied.
BoutReal wallInwardDerivative(BoutReal phi_midpoint, BoutReal phi0, BoutReal dx0);

/// Impose DC midpoint Dirichlet and the selected AC radial condition on chi.
///
/// Each radial guard layer reflects the corresponding interior layer's Z mean
/// about the supplied DC value. AC fluctuations are retained for Neumann
/// or reflected for zero-midpoint Dirichlet.
/// Interior cells are preserved. All local Y rows are processed at physical X
/// boundaries; periodic X leaves the field unchanged.
///
/// @param[in,out] chi Auxiliary potential whose radial guard cells are updated.
/// @param ac Radial condition on Z-dependent fluctuations.
/// @param dc_inner Inner wall DC values in chi's units, indexed by local Y.
///                 If supplied, the array must contain LocalNy entries;
///                 nullptr specifies a zero DC midpoint.
/// @param dc_outer Outer wall DC values, with the same convention as dc_inner.
/// @throws BoutException If a nonperiodic X mesh lacks complete periodic Z lines.
void applyPolarisationBoundary(Field3D& chi, PolarisationACBoundary ac,
                               const std::vector<BoutReal>* dc_inner = nullptr,
                               const std::vector<BoutReal>* dc_outer = nullptr);

/// Approximate the auxiliary potential with a fixed sequence of spectral solves.
///
/// Starting from zero, perform one initial solve plus the requested number of
/// finite-volume defect corrections. DC wall values enter through the FV residual;
/// spectral corrections use homogeneous boundaries. The result depends on these inputs
/// and the mesh, with no initial guess retained from earlier RHS evaluations.
///
/// @param coefficient Boussinesq coefficient c in L(chi) = div(c grad_perp chi).
/// @param rhs Physical vorticity rate to reproduce with L(chi), excluding encoded
///            wall-state rates. Only physical interior cells supply the residual.
/// @param[in,out] inverse Spectral approximate inverse configured for coefficient
///                       and the appropriate homogeneous correction boundaries.
/// @param corrections Nonnegative number of FV corrections after the first solve.
/// @param ac Radial condition on Z-dependent fluctuations.
/// @param dc_inner Optional inner wall DC values, as in applyPolarisationBoundary.
/// @param dc_outer Optional outer wall DC values, as in applyPolarisationBoundary.
/// @returns A new chi field with the requested radial guards imposed. A finite
///          correction count can leave a nonzero FV residual; coefficient and rhs
///          fields are preserved.
/// @throws BoutException If corrections is negative or the boundary mesh is invalid.
Field3D fixedPolarisationPotential(const Field3D& coefficient, const Field3D& rhs,
                                   Laplacian& inverse, int corrections,
                                   PolarisationACBoundary ac,
                                   const std::vector<BoutReal>* dc_inner = nullptr,
                                   const std::vector<BoutReal>* dc_outer = nullptr);
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
struct EvolveVorticityIonEnergy : public NamedComponent<EvolveVorticityIonEnergy> {
  ///
  /// # Inputs
  ///
  /// - <name>
  ///   - radial_boundary      fixed (default), free_slip or no_slip
  ///   - phi_boundary_velocity  Required positive velocity [m/s] for no_slip
  ///   - polarisation_ac_boundary  neumann (default) or dirichlet
  ///   - average_atomic_mass  Boussinesq reference atomic mass, default 2.0
  ///   - density_floor        Minimum normalised total ion density, default 1e-7
  ///   - diagnose             Output additional diagnostic fields? Default false
  ///   - poloidal_flows       Include poloidal ExB flows? Default false
  ///   - viscosity            Kinematic viscosity [m^2/s], default 0.0
  ///   - phi_max_iterations   Recovery iteration limit, default 100
  ///   - phi_atol             Potential convergence scale [eV], default 1.0
  ///   - phi_rtol             Recovery relative tolerance, default 1e-3
  ///   - vort_atol            Normalised vorticity residual tolerance, default 1e-10
  ///   - polarisation_corrections  FV auxiliary defect corrections, default 2
  ///   - laplacian            Approximate potential inverse options
  ///   - polarisation_laplacian  Auxiliary inverse options; type defaults to laplacian
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

  BoutReal average_atomic_mass; ///< Weighted average atomic mass
  bool poloidal_flows;          ///< Include poloidal (Y) flows?
  BoutReal density_floor;       ///< Used to calculate Ti

  bool include_viscosity;                ///< Was kinematic viscosity specified?
  Coordinates::FieldMetric viscosity;    ///< Normalised kinematic viscosity
  Field3D viscous_heating{0.0};          ///< Face-based -nu grad(Vort) dot grad(Phi)
  Field3D electron_viscous_heating{0.0}; ///< Electron mass fraction of viscous heating

  Field3D phi{0.0}; ///< Electrostatic potential ϕ
  Field3D chi{0.0}; ///< Auxiliary current potential with selected AC boundaries
  Field3D Pi_sum;   ///< Total ion pressure

  // Starting guess for potential recovery, updated only at outputs
  Field3D phi_guess{0.0};         ///< Potential at the last output state
  Field3D ion_KE_perp_guess{0.0}; ///< K⊥ at the last output state
  bool have_solution{false};      ///< phi and ion_KE_perp hold a converged solution
  int phi_iterations{0};          ///< Iterations used in the last potential recovery
  Field3D psi;                    ///< Pressure potential Σ A_s p_s / (Abar Z_s)
  Field3D Ti;                     ///< Common ion temperature
  Coordinates::FieldMetric boussinesq_coeff; ///< c = mass / B^2
  std::unique_ptr<Laplacian> phiSolver;      ///< Laplacian solver in X-Z
  std::unique_ptr<Laplacian> chiSolver;      ///< Mixed-boundary auxiliary inverse
  int chi_iterations{0};           ///< Spectral solves in the last chi calculation
  BoutReal chi_residual{0.0};      ///< Max FV residual of the last chi calculation
  int polarisation_corrections{2}; ///< FV defect corrections after the spectral solve

  bout::ion_energy::RadialBoundaryModel radial_model;
  std::vector<BoutReal> wall_rate_inner; ///< Current inner wall dphi/dt, by local y
  std::vector<BoutReal> wall_rate_outer; ///< Current outer wall dphi/dt, by local y

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
  Field3D DivJpol; ///< Polarisation transport divergence with the selected wall closure
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

  /// Encode the initial or legacy-restart phi wall values in the solver state.
  ///
  /// Reads phi's interior and adjacent radial guard values and averages their
  /// midpoints over Z (zero on a fresh start). The inner closed-core midpoint is also
  /// averaged over the connected Y surface using wall-face J dy weights.
  /// Writes that midpoint to every radial Vort guard column and Z cell at each
  /// physical Y row, after zeroing all Vort boundary cells. Interior Vort is
  /// preserved; Y boundary storage remains zero. Used for every wall mode.
  void initialiseBoundaryState();

  /// Restore phi's radial midpoint conditions during potential recovery.
  ///
  /// Reads the Z mean of Vort's first radial guard column as the stored wall
  /// value, using its J dy-weighted surface average at an averaged core wall.
  /// Writes phi_guard = 2 * midpoint - phi_interior and copies that value into
  /// the remaining radial guards at physical Y rows. Neither interior phi nor
  /// solver Vort is modified. Does nothing for periodic X.
  void applyPotentialBoundary();

  /// Compute the current relaxation rates before the auxiliary chi solve.
  ///
  /// Reads the encoded Vort midpoints, the converged phi in the first interior
  /// cell, the model's normalised phi_boundary_velocity, and mesh dx/g11. Uses Z means
  /// and the centred inward face derivative to compute kappa * grad_n(phi).
  /// At an averaged core wall, both the midpoint and physical derivative use
  /// wall-face J dy-weighted Y averages over the connected surface.
  ///
  /// Overwrites wall_rate_inner and wall_rate_outer, indexed by local Y, in
  /// normalised potential/time units. Unused walls and Y guard entries are zero.
  /// These rates supply chi's DC wall values and finally()'s boundary RHS.
  /// Fixed and free-slip rates are zero. Preserves phi/Vort
  /// and recomputes every rate without retaining earlier RHS results.
  void computeWallRates();

  /// Write the precomputed midpoint rates into the time solver's boundary RHS.
  ///
  /// Reads wall_rate_inner and wall_rate_outer from the current transform_impl()
  /// evaluation. Writes each rate into every physical radial ddt(Vort) guard
  /// column and Z cell at its local physical Y row. Preserves interior and Y
  /// boundary derivatives, which finally() must initialise separately. No wall
  /// physics or MPI reduction is performed here; periodic X leaves the RHS alone.
  void writeWallRates();

  /// Shared viscous current, energy flux and heating, using physical vorticity guards.
  void calculateViscosity();

  /// Sum thermal sources and kinetic energy changes from momentum and particle sources
  Field3D sumIonEnergySources(const Options& allspecies) const;

  /// Use the Laplacian as a preconditioner for a finite-volume residual step.
  Field3D correctPotential(const Field3D& rhs, const Field3D& potential) const;
  Field3D correctPotential(const Field3D& rhs, const Field3D& potential,
                           Laplacian& approximate_inverse) const;

  /// Recover the auxiliary current potential and update its solver diagnostics.
  ///
  /// @param rhs Physical interior vorticity rate in normalised units, excluding
  ///            the encoded radial wall-state rates.
  ///
  /// Reads boussinesq_coeff, polarisation_corrections and the configured inverse.
  /// In relaxing mode it also reads the current wall_rate_inner/outer as chi's
  /// DC midpoint values, approximating d(Phi)/dt by d(phi_b)/dt at the wall and
  /// neglecting d(psi_bar_b)/dt. Fixed boundaries use zero DC values.
  /// Overwrites chi, chi_iterations and chi_residual using a fixed solve sequence
  /// from zero, so the result depends on current inputs rather than RHS history.
  void solvePolarisationPotential(const Field3D& rhs);

  /// Calculate ion perpendicular kinetic energy, K⊥ = 1/2 c |∇⊥Φ|^2
  ///
  /// This should be consistent with the vorticity discretisation,
  /// so that the integration by parts of Vort * Phi over the domain
  /// gives the kinetic perpendicular energy plus boundary terms.
  ///
  Field3D calculateKEperp(const Field3D& phi, const Field3D& psi) const;

  /// Store the most recent converged potential and perpendicular kinetic energy
  /// as the starting guess for subsequent right-hand-side evaluations.
  ///
  /// Called from outputVars after an RHS evaluation at the accepted output state.
  /// The guess remains fixed between outputs, including during each step's
  /// nonlinear and linear iterations, so the RHS is deterministic within a step.
  void acceptWarmStart();
};

namespace {
RegisterComponent<EvolveVorticityIonEnergy> registercomponentevolvevorticityionenergy;
}

#endif // EVOLVE_VORTICITY_ION_ENERGY_H
