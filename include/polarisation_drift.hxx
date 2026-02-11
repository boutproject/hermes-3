#pragma once
#ifndef POLARISATION_DRIFT_H
#define POLARISATION_DRIFT_H

#include "component.hxx"

class Laplacian;

/// Calculates polarisation drift terms for all charged species, both
/// ions and electrons.
///
/// Approximates the polarisation drift by a generalised flow potential `phi_pol`
///
///  v_pol = - (A / (Z * B^2)) * Grad_perp(phi_pol)
///
/// phi_pol is approximately the time derivative of the electric potential
/// in the frame of the flow, plus an ion diamagnetic contribution
///
/// phi_pol is calculated using:
///
/// Div(mass_density / B^2 * Grad_perp(phi_pol)) = Div(Jpar) + Div(Jdia) + ...
///
/// Where the divergence of currents on the right is calculated from:
///  - species[...]["momentum"] The parallel momentum of charged species
///  - DivJdia,   diamagnetic current, calculated in vorticity component
///  - DivJcol    collisional current, calculated in vorticity component
///  - DivJextra  Other currents, eg. 2D parallel closures
///
/// The mass_density quantity is the sum of density * atomic mass for all
/// charged species (ions and electrons)
struct PolarisationDrift : public Component {
  //
  PolarisationDrift(std::string name, Options &options, Solver *UNUSED(solver));

  void outputVars(Options &state) override;

  // The following functions are public for unit testing

  /// Calculate divergence of all currents except polarisation
  Field3D calcDivJ(GuardedOptions& state);

  /// Calculate energy exchange term nonlinear in pressure
  /// due to compression of polarisation drift
  ///   (3 / 2) ddt(Pi) += (Pi * m_i / n0 / Z) * DivJ
  ///
  /// Adds energy_source for all species that have charge, mass and pressure
  /// Throws a BoutException if boussinesq=true
  void diamagneticCompression(GuardedOptions& state, Field3D DivJ);

  /// Solve for time derivative of potential
  /// Using Div(mass_density / B^2 Grad_perp(dphi/dt)) = DivJ
  Field3D calcMassDensity(GuardedOptions& state);

  /// Calculate poloidal drift potential-flow approximation
  Field3D calcPolFlowPotential(Field3D mass_density, Field3D DivJ);

  /// Polarisation drift approximated by a potential flow
  ///
  /// v_p = - (m_i / (Z_i * B^2)) * Grad(phi_pol)
  ///
  /// Sets density_source, energy_source and momentum_source
  /// for all species with mass and charge.
  void polarisationAdvection(GuardedOptions& state, Field3D phi_pol);
private:
  std::unique_ptr<Laplacian> phiSolver; // Laplacian solver in X-Z

  Field2D Bsq; // Cached SQ(coord->Bxy)
  
  // Diagnostic outputs
  bool diagnose; ///< Save diagnostic outputs?
  Field3D DivJ; ///< Divergence of all other currents
  Field3D phi_pol; ///< Polarisation drift potential
  Options diagnostics; ///< Other diagnostic outputs

  bool boussinesq; // If true, assume a constant mass density in Jpol
  BoutReal average_atomic_mass; // If boussinesq=true, mass density to use
  BoutReal density_floor; // Minimum mass density if boussinesq=false
  bool advection; // Advect fluids by an approximate polarisation velocity?
  bool diamagnetic_polarisation; // Calculate compression terms?

  /// Inputs
  ///
  /// - species
  ///   - ...  All species with both charge and mass
  ///     - AA
  ///     - charge
  ///     - density
  ///     - momentum (optional)
  ///
  /// - fields
  ///   - DivJextra  (optional)
  ///   - DivJdia    (optional)
  ///   - DivJcol    (optional)
  ///
  /// Sets
  ///
  /// - species
  ///   - ...  All species with both charge and mass
  ///     - density_source
  ///     - energy_source    (if pressure set)
  ///     - momentum_source  (if momentum set)
  ///
  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<PolarisationDrift> registercomponentpolarisationdrift("polarisation_drift");
}

#endif // POLARISATION_DRIFT_H
