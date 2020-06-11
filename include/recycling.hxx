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
  ///     - recycle_multiplier   The recycled flux multiplier
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

private:

  struct RecycleChannel {
    std::string from; ///< The species name to recycle
    std::string to;   ///< Species to recycle to
    BoutReal multiplier; ///< Flux multiplier. Combination of recycling fraction and species
                         ///< change e.g h+ -> h2 results in 0.5 multiplier
  };

  std::vector<RecycleChannel> channels; // Recycling channels
  
};

namespace {
RegisterComponent<Recycling> registercomponentrecycling("recycling");
}

#endif // RECYCLING_H
