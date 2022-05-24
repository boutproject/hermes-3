#pragma once
#ifndef SOLKIT_NEUTRAL_PARALLEL_DIFFUSION_H
#define SOLKIT_NEUTRAL_PARALLEL_DIFFUSION_H

#include "component.hxx"

/// Add effective diffusion of neutrals in a 1D system
///
/// This version is intended to match the calculation of neutral diffusion
/// in SOL-KiT (ca 2022).
/// 
struct SOLKITNeutralParallelDiffusion : public Component {
  ///
  /// alloptions
  ///   - units
  ///     - eV
  ///     - meters
  ///     - inv_meters_cubed
  ///   - <name>
  ///     - neutral_temperature [eV]
  ///     
  SOLKITNeutralParallelDiffusion(std::string name, Options &alloptions, Solver *) {
    auto Tnorm = get<BoutReal>(alloptions["units"]["eV"]);
    auto& options = alloptions[name];
    neutral_temperature = options["neutral_temperature"]
      .doc("Neutral atom temperature [eV]")
      .withDefault(3.0)
      / Tnorm; // Normalise

    auto Nnorm = get<BoutReal>(alloptions["units"]["inv_meters_cubed"]);
    auto rho_s0 = get<BoutReal>(alloptions["units"]["meters"]);
    area_norm = 1. / (Nnorm * rho_s0);
  }

  ///
  /// Inputs
  ///  - species
  ///    - <all neutrals>    # Applies to all neutral species
  ///      - AA
  ///      - density
  ///
  /// Sets
  ///  - species
  ///    - <name>
  ///      - density_source
  void transform(Options &state) override;

private:
  BoutReal neutral_temperature;  ///< Fixed neutral t
  BoutReal area_norm; ///< Area normalisation [m^2]
};

namespace {
RegisterComponent<SOLKITNeutralParallelDiffusion>
    register_solkit_neutral_parallel_diffusion("solkit_neutral_parallel_diffusion");
}

#endif // SOLKIT_NEUTRAL_PARALLEL_DIFFUSION_H
