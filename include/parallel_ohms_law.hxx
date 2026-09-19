#pragma once
#ifndef PARALLEL_OHMS_LAW
#define PARALLEL_OHMS_LAW

#include "component.hxx"

/// Caution! This component is not carfully evaluated for multispecies cases with multiple ions.

/// Use parallel Ohm's law to calculate the parallel current. 
/// This law basically comes from the electron parallel momentum equation 
/// by neglecting inertial (small by m_e / m_i) and viscous terms 
/// (see Phys. Plasmas 10, 4744–4757 (2003) by Simakov & Catto).
/// Use the parallel current to calculate the electron parallel velocity
///
///   Jpar = - 1/eta * ∇phi + 1/(eta e n_e) * ∇P_e + (0.71 k_B)/(eta e) * ∇T_e  
///
/// where phi is the electric potential and eta is the resistivity.
///
/// After normalization the equation takes the form: 
///
///   Jpar = - ∇phi / eta  +  ∇P_e / (eta n_e) + 0.71 ∇T_e / eta 
///
/// where the normalization parameter for the resistivity eta is:
///
///   eta_norm = Tnorm / (rho_s0 * Nnorm * Cs0 * qe)
///
/// Then uses this parallel current to calculate the parallel electron velocity.
///
/// Note: This needs to be called after collisions and other
///       components which impose forces on electrons
///
struct ParallelOhmsLaw : public Component {

  ParallelOhmsLaw(std::string name, Options& alloptions, Solver*);

  /// Required inputs
  /// - species
  ///   - e
  ///     - pressure
  ///     - density
  ///     - momentum_source [optional]
  ///     Asserts that charge = -1
  ///
  /// Sets in the input 
  /// - species
  ///   - <only for e>  if phi is set.  if both density and charge are set
  ///     - electron parallel velocity
  /// 

  void calculateResistivity(Options &options, Field3D &Ne);

  void transform(Options &state) override;

  void finally(const Options &state) override {
    // Get the velocity with boundary condition applied.
    // This is for output only
    Ve = get<Field3D>(state["species"]["e"]["velocity"]);
    NVe = get<Field3D>(state["species"]["e"]["momentum"]);
  }

  /// Save output diagnostics
  void outputVars(Options& state) override;
private:
  std::string name; ///< Name of this species
  bool diagnose; ///< Output additional fields
  BoutReal resistivity_floor;

  BoutReal Tnorm;    // Temperature normalisation [eV]
  BoutReal Nnorm;    // Density normalisation [m^-3]
  BoutReal rho_s0;   // Length normalisation [m]
  BoutReal Omega_ci; // Frequency normalisation [s^-1]
  BoutReal Cs0;      // Velocity normalisation [m/s]

  bool spitzer_resistivity;  // Use Spitzer formula for resistivity. Otherwise, get collision frequency from collisions component. 

  Field3D jpar; ///< Parallel current
  Field3D eta; ///< Parallel plasma resistivity

  Field3D Ve;   ///< Electron parallel velocity
  Field3D NVe;  ///< Electron parallel momentum (normalised)

};

namespace {
RegisterComponent<ParallelOhmsLaw> registercomponentparallelohmslaw("parallel_ohms_law");
}

#endif // ELECTRON_FORCE_BALANCE
