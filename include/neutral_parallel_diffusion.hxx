#pragma once
#ifndef NEUTRAL_PARALLEL_DIFFUSION_H
#define NEUTRAL_PARALLEL_DIFFUSION_H

#include "component.hxx"

/// Add effective diffusion of neutrals in a 1D system,
/// by projecting cross-field diffusion onto parallel distance.
struct NeutralParallelDiffusion : public Component {
  NeutralParallelDiffusion(std::string name, Options &alloptions, Solver *) : name(name) {
    auto& options = alloptions[name];
    dneut = options["dneut"]
                .doc("cross-field diffusion projection (B  / Bpol)^2")
                .as<BoutReal>();
  }

  ///
  /// Inputs
  ///  - species
  ///    - <name>
  ///      - AA
  ///      - collision_frequency
  ///      - density
  ///      - temperature
  ///      - pressure
  ///      - velocity     [optional]
  ///      - momentum     [if velocity set]
  ///
  /// Sets
  ///  - species
  ///    - <name>
  ///      - density_source
  ///      - energy_source
  ///      - momentum_source  [if velocity set]
  void transform(Options &state) override;

private:
  std::string name; ///< The name of the neutral species
  BoutReal dneut; ///< cross-field diffusion projection (B  / Bpol)^2
};

namespace {
RegisterComponent<NeutralParallelDiffusion>
    register_component_neutral_parallel_diffusion("neutral_parallel_diffusion");
}

#endif // NEUTRAL_PARALLEL_DIFFUSION_H
