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
  void finally(const Options &state) override;
private:
  std::string name;     ///< Short name of species e.g "e"

  Field3D NV;           ///< Species parallel momentum (normalised, evolving)
  
  bool bndry_flux;      // Allow flows through boundaries?
  bool poloidal_flows;  // Include ExB flow in Y direction?
};

namespace {
RegisterComponent<EvolveMomentum> registercomponentevolvemomentum("evolve_momentum");
}

#endif // EVOLVE_MOMENTUM_H
