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
  void transform(Options &state) override;

  void outputVars(Options &state) override;
private:
  std::unique_ptr<Laplacian> phiSolver; // Laplacian solver in X-Z

  Field3D Bsq; // Cached SQ(coord->Bxy)
  
  // Diagnostic outputs
  bool diagnose; ///< Save diagnostic outputs?
  Field3D DivJ; //< Divergence of all other currents
  Field3D phi_pol; //< Polarisation drift potential

  bool boussinesq; // If true, assume a constant mass density in Jpol
  BoutReal average_atomic_mass; // If boussinesq=true, mass density to use
  BoutReal density_floor; // Minimum mass density if boussinesq=false
  bool advection; // Advect fluids by an approximate polarisation velocity?
  bool diamagnetic_polarisation; // Calculate compression terms?
};

namespace {
RegisterComponent<PolarisationDrift> registercomponentpolarisationdrift("polarisation_drift");
}

#endif // POLARISATION_DRIFT_H
