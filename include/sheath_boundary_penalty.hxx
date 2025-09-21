#pragma once
#ifndef SHEATH_BOUNDARY_PENALTY_H
#define SHEATH_BOUNDARY_PENALTY_H

#include "component.hxx"
#include <bout/region.hxx>

/// Penalty method for imposing sheath boundary conditions over an
/// arbitrary shaped wall that is immersed inside the structured mesh.
///
/// Notes:
///  - This should be applied AFTER sheath boundaries have been applied
///    An error should be raised if this is done in the wrong order.
///
struct SheathBoundaryPenalty : public Component {
  /// # Input options
  /// - <name> e.g. "sheath_boundary_penalty"
  ///   - penalty_timescale    Timescale in seconds
  ///
  /// # Reads from the mesh
  /// - penalty_mask[x,y,z]    A 3D field defining the shape of the boundary
  ///                          Equal 0 inside the plasma, 1 in the wall
  SheathBoundaryPenalty(std::string name, Options& options, Solver*);

  /// # Inputs
  ///   - fields
  ///     - phi  [optional]
  ///       If present then electron sheath current terms are calculated
  ///   - species
  ///     - <name>
  ///       - AA
  ///       - charge
  ///       - density
  ///       - temperature
  ///       - pressure [optional]
  ///       - velocity [optional]
  ///       - momentum [optional]
  ///
  /// # Outputs
  ///   - species
  ///     - <name>
  ///       - density_source     Adds to existing
  ///       - momentum_source    Adds to existing
  ///       - energy_source      Adds to existing
  ///       - density_penalty    Particle source term (negative)
  ///       - momentum_penalty   Momentum source term
  ///       - energy_penalty     Energy source term
  ///
  /// The *_penalty terms can be used in the recycling component
  /// to implement volumetric recycling.
  void transform(Options& state) override;

  /// Save diagnostics
  ///
  /// Always saves `penalty_mask` as a time-independent 3D field
  ///
  /// If `diagnose = true`, also saves:
  ///  - S<species>_penalty    Particle source (negative)
  ///  - F<species>_penalty    Momentum source
  ///  - R<species>_penalty    Energy source
  ///
  void outputVars(Options& state) override;

private:
  /// Mask function that is 0 in the plasma, 1 in the wall
  Field3D penalty_mask;

  /// Cell indices where penalty_mask > 0
  Region<Ind3D> penalty_region;

  BoutReal gamma_e, gamma_i; // Sheath heat transmission factors

  /// Timescale of penalisation [normalised]
  BoutReal penalty_timescale;

  /// Diagnostics for output
  Options diagnostics;
  bool diagnose; ///< Save diagnostics?
};

namespace {
RegisterComponent<SheathBoundaryPenalty>
    registercomponentsheathboundarypenalty("sheath_boundary_penalty");
}

#endif // SHEATH_BOUNDARY_PENALTY_H
