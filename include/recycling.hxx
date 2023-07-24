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
    BoutReal multiplier; ///< Flux multiplier. Combination of recycling fraction and species
                         ///< change e.g h+ -> h2 results in 0.5 multiplier
    BoutReal energy; ///< Energy of recycled particle (normalised to Tnorm)
  };

  std::vector<RecycleChannel> channels; // Recycling channels
  bool sol_recycling;  ///< Flag for reycling in the edge
  bool diagnose; ///< Outputting diagnostics?
  Field3D density_source, energy_source; ///< Recycling particle and energy sources for all locations
  Field3D sol_recycling_density_source, sol_recycling_energy_source;  ///< Recycling particle and energy sources for edge recycling only
  Field3D target_recycling_density_source, target_recycling_energy_source;  ///< Recycling particle and energy sources for target recycling only
  Field3D radial_particle_flow, radial_energy_flow;  ///< Radial fluxes coming from evolve_density and evolve_pressure used in recycling calc
  
};

namespace {
RegisterComponent<Recycling> registercomponentrecycling("recycling");
}

#endif // RECYCLING_H
