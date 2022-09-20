#pragma once
#ifndef NEUTRAL_PARALLEL_DIFFUSION_H
#define NEUTRAL_PARALLEL_DIFFUSION_H

#include "component.hxx"

/// Add effective diffusion of neutrals in a 1D system,
/// by projecting cross-field diffusion onto parallel distance.
///
/// Note: This needs to be calculated after the collision frequency,
/// so is a collective component. This therefore applies diffusion
/// to all neutral species i.e. those with no (or zero) charge
///
/// Diagnostics
/// -----------
///
/// If diagnose = true then the following outputs are saved for
/// each neutral species
///
///  - D<name>_Dpar   Parallel diffusion coefficient e.g. Dhe_Dpar
///  - S<name>_Dpar   Density source due to diffusion
///  - E<name>_Dpar   Energy source due to diffusion
///  - F<name>_Dpar   Momentum source due to diffusion
///
struct NeutralParallelDiffusion : public Component {
  NeutralParallelDiffusion(std::string name, Options &alloptions, Solver *) {
    auto& options = alloptions[name];
    dneut = options["dneut"]
                .doc("cross-field diffusion projection (B  / Bpol)^2")
                .as<BoutReal>();

    diagnose = options["diagnose"]
      .doc("Output additional diagnostics?")
      .withDefault<bool>(false);
  }

  ///
  /// Inputs
  ///  - species
  ///    - <all neutrals>    # Applies to all neutral species
  ///      - AA
  ///      - collision_frequency
  ///      - density
  ///      - temperature
  ///      - pressure     [optional, or density * temperature]
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
  BoutReal dneut; ///< cross-field diffusion projection (B  / Bpol)^2

  bool diagnose; ///< Output diagnostics?

  /// Per-species diagnostics
  struct Diagnostics {
    Field3D Dn; ///< Diffusion coefficient
    Field3D S; ///< Particle source
    Field3D E; ///< Energy source
    Field3D F; ///< Momentum source
  };

  /// Store diagnostics for each species
  std::map<std::string, Diagnostics> diagnostics;
};

namespace {
RegisterComponent<NeutralParallelDiffusion>
    register_component_neutral_parallel_diffusion("neutral_parallel_diffusion");
}

#endif // NEUTRAL_PARALLEL_DIFFUSION_H
