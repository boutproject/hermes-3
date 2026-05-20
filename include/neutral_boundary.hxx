#pragma once
#ifndef NEUTRAL_BOUNDARY_H
#define NEUTRAL_BOUNDARY_H

#include "component.hxx"

/// Per-species boundary condition for neutral particles at
/// sheath (Y) boundaries.
///
/// Sets boundary conditions:
/// - Free boundary conditions on logarithm
///   of density, temperature and pressure
/// - No-flow boundary conditions on velocity
///   and momentum.
///
/// Adds an energy sink corresponding to a flux of heat to the walls.
///
/// Heat flux into the wall is
///   q = gamma_heat * n * T * v_th
///
/// where v_th = sqrt(eT/m) is the thermal speed
///

/// Core ionising: neutrals that reach the core boundary are removed from the domain 
/// and their flux is added to the boundary condition on the flux of ions from the core.
/// 
/// Inputs
///
///   - <name>
///     - species    A comma-separated list of species to ionise
///   - <species>
///     - ionise_as  The species to ionise into
///     - ionise_multiplier   The ionised flux multiplier, between 0 and 1
///     - ionise_energy       The energy of the ionised particles [eV]
///
struct NeutralBoundary : public Component {
  NeutralBoundary(std::string name, Options& options, Solver*);

  void outputVars(Options& state) override;

private:

  std::string name; ///< Short name of species e.g "d"

  BoutReal Tnorm; // Temperature normalisation [eV]

  BoutReal target_energy_refl_factor, sol_energy_refl_factor,
      pfr_energy_refl_factor; ///< Fraction of energy retained after reflection
  BoutReal target_fast_refl_fraction, sol_fast_refl_fraction,
      pfr_fast_refl_fraction; ///< Fraction of neutrals undergoing fast reflection

  Field3D target_energy_source, wall_energy_source; ///< Diagnostic for power loss

  bool diagnose; ///> Save diagnostic variables?

  bool lower_y; ///< Boundary condition at lower y?
  bool upper_y; ///< Boundary condition at upper y?
  bool sol;     ///< Boundary condition at sol?
  bool pfr;     ///< Boundary condition at pfr?

  ///
  /// state
  ///  - species
  ///    - <name>
  ///      - density       Free boundary
  ///      - temperature   Free boundary
  ///      - pressure      Free boundary
  ///      - velocity [if set] Zero boundary
  ///      - momentum [if set] Zero boundary
  ///      - energy_source  Adds wall losses
  ///

  /// Core ionising 
  struct IoniseChannel {
    std::string from; ///< The species name to be ionised (neutral)
    std::string to;   ///< Species to ionising to (ion)

    /// Flux multiplier (ionising fraction).
    BoutReal core_multiplier;
    BoutReal core_energy; ///< Energy of ionisingd particle (normalised to Tnorm)

    // Ionising particle and energy sources for the different sources of ionising
    // These sources are per-channel and added to the `to` species
    /// Ionising density source for core ionising
    Field3D core_ionising_density_source = 0.0;
    /// Ionising energy source for core ionising
    Field3D core_ionising_energy_source = 0.0;
  };
  std::vector<IoniseChannel> channels; // Ionising channels

  bool core_ionising{false}; ///< Flags for enabling ionisinging in different regions

  BoutReal density_floor,
      pressure_floor; ///< minimum values for Nn, Pn to avoid divide by zero

  BoutReal gamma_i; /// Sheath heat transmission coefficient

  Field3D density_source, energy_source; ///< Ionising particle and energy sources for all locations

  Field3D energy_flow_xlow;   ///< Cell edge fluxes used for calculating
  Field3D particle_flow_xlow; ///< Radial core particle fluxes for ionising calc. 
                              ///< Does it need y flow?


  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<NeutralBoundary> registercomponentneutralboundary("neutral_boundary");
}

#endif // NEUTRAL_BOUNDARY_H
