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
    BoutReal target_multiplier, sol_multiplier, pfr_multiplier; 
    BoutReal target_energy, sol_energy, pfr_energy; ///< Energy of recycled particle (normalised to Tnorm)
    BoutReal target_fast_recycle_fraction, pfr_fast_recycle_fraction, sol_fast_recycle_fraction;   ///< Fraction of ions undergoing fast reflection
    BoutReal target_fast_recycle_energy_factor, sol_fast_recycle_energy_factor, pfr_fast_recycle_energy_factor;   ///< Fraction of energy retained by fast recycled neutrals
  };

  std::vector<RecycleChannel> channels; // Recycling channels

  bool target_recycle, sol_recycle, pfr_recycle;  ///< Flags for enabling recycling in different regions
  bool diagnose; ///< Save additional post-processing variables?

  Field3D density_source, energy_source; ///< Recycling particle and energy sources for all locations
  Field3D heatflow_ylow, heatflow_xlow; ///< Cell edge fluxes used for calculating fast recycling energy source

  Field3D target_recycle_density_source, target_recycle_energy_source;  ///< Recycling particle and energy sources for target recycling only
  Field3D wall_recycle_density_source, wall_recycle_energy_source;  ///< Recycling particle and energy sources for pfr + sol recycling

  Field3D radial_particle_outflow, radial_energy_outflow;  ///< Radial fluxes coming from evolve_density and evolve_pressure used in recycling calc
  

};

namespace {
RegisterComponent<Recycling> registercomponentrecycling("recycling");
}

#endif // RECYCLING_H
