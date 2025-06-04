#pragma once
#ifndef RECYCLING_H
#define RECYCLING_H

#include "component.hxx"

/// Convert fluxes of species at boundaries
///
/// Since this must be calculated after boundary fluxes (e.g. sheath),
/// it is included as a top-level component
/// 
struct Recycling : public Component {
  
  /// Inputs
  /// 
  ///   - <name>
  ///     - species    A comma-separated list of species to recycle
  ///   - <species>
  ///     - recycle_as  The species to recycle into
  ///     - recycle_multiplier   The recycled flux multiplier, between 0 and 1
  ///     - recycle_energy       The energy of the recycled particles [eV]
  ///
  Recycling(std::string name, Options &alloptions, Solver *);

  /// Inputs
  ///
  /// - species
  ///   - <species>
  ///    - density
  ///    - velocity
  ///
  /// Outputs
  ///
  /// - species
  ///  - <species>
  ///   - density_source
  ///
  void transform(Options &state) override;
  void outputVars(Options &state) override;

private:

  struct RecycleChannel {
    std::string from; ///< The species name to recycle
    std::string to;   ///< Species to recycle to

    /// Flux multiplier (recycling fraction). 
    /// Combination of recycling fraction and species change e.g h+ -> h2 results in 0.5 multiplier
    BoutReal target_multiplier, sol_multiplier, pfr_multiplier, pump_multiplier; 
    BoutReal target_energy, sol_energy, pfr_energy; ///< Energy of recycled particle (normalised to Tnorm)
    BoutReal target_fast_recycle_fraction, pfr_fast_recycle_fraction, sol_fast_recycle_fraction;   ///< Fraction of ions undergoing fast reflection
    BoutReal target_fast_recycle_energy_factor, sol_fast_recycle_energy_factor, pfr_fast_recycle_energy_factor;   ///< Fraction of energy retained by fast recycled neutrals

    // Recycling particle and energy sources for the different sources of recycling
    // These sources are per-channel and added to the `to` species
    Field3D target_recycle_density_source, target_recycle_energy_source;
    Field3D wall_recycle_density_source, wall_recycle_energy_source;  ///< Recycling particle and energy sources for pfr + sol recycling
    Field3D pump_density_source, pump_energy_source;  ///< Recycling particle and energy sources for pump recycling
  };

  std::vector<RecycleChannel> channels; // Recycling channels

  bool target_recycle, sol_recycle, pfr_recycle, neutral_pump;  ///< Flags for enabling recycling in different regions
  bool diagnose; ///< Save additional post-processing variables?

  BoutReal density_floor, pressure_floor; ///< minimum values for Nn, Pn to avoid divide by zero

  Field3D density_source, energy_source; ///< Recycling particle and energy sources for all locations
  Field3D energy_flow_ylow, energy_flow_xlow; ///< Cell edge fluxes used for calculating fast recycling energy source
  Field3D particle_flow_xlow; ///< Radial wall particle fluxes for recycling calc. No need to get poloidal from here, it's calculated from sheath velocity

  Field2D is_pump; ///< 1 = pump, 0 = no pump. Works only in SOL/PFR. Provided by user in grid file.
};

namespace {
RegisterComponent<Recycling> registercomponentrecycling("recycling");
}

#endif // RECYCLING_H
