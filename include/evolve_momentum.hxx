#pragma once
#ifndef EVOLVE_MOMENTUM_H
#define EVOLVE_MOMENTUM_H

#include "component.hxx"

/// Evolve parallel momentum
struct EvolveMomentum : public Component {
  EvolveMomentum(std::string name, Options &options, Solver *solver);

  /// This sets in the state
  /// - species
  ///   - <name>
  ///     - momentum
  ///     - velocity  if density is defined
  void transform(Options &state) override;
  
  /// Calculate ddt(NV).
  ///
  /// Inputs
  /// - species
  ///   - <name>
  ///     - density
  ///     - velocity
  ///     - pressure (optional)
  ///     - momentum_source (optional)
  ///     - sound_speed (optional, used for numerical dissipation)
  ///     - temperature (only needed if sound_speed not provided)
  /// - fields
  ///   - phi (optional)
  void finally(const Options &state) override;

  void outputVars(Options &state) override;
private:
  std::string name;     ///< Short name of species e.g "e"

  Field3D NV;           ///< Species parallel momentum (normalised, evolving)
  Field3D NV_solver;    ///< Momentum as input from solver
  Field3D V;            ///< Species parallel velocity

  Field3D momentum_source; ///< From other components. Stored for diagnostic output

  bool bndry_flux;      // Allow flows through boundaries?
  bool poloidal_flows;  // Include ExB flow in Y direction?

  BoutReal density_floor;

  BoutReal hyper_z;  ///< Hyper-diffusion

  bool diagnose; ///< Output additional diagnostics?
};

namespace {
RegisterComponent<EvolveMomentum> registercomponentevolvemomentum("evolve_momentum");
}

#endif // EVOLVE_MOMENTUM_H
