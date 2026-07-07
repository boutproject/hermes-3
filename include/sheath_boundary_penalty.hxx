#pragma once
#ifndef SHEATH_BOUNDARY_PENALTY_H
#define SHEATH_BOUNDARY_PENALTY_H

#include "component.hxx"
#include <bout/mesh.hxx>
#include <bout/region.hxx>

/// Penalty method for imposing sheath boundary conditions over an
/// arbitrary shaped wall that is immersed inside the structured mesh.
///
/// Notes:
///  - This should be applied AFTER sheath boundaries have been applied
///    An error should be raised if this is done in the wrong order.
///
struct SheathBoundaryPenalty : public Component {
  struct PenaltyMaskData {
    Field3D mask;
    Region<Ind3D> region;
  };

  struct PenaltySourceData {
    Field3D density;
    Field3D momentum;
    Field3D energy;
  };

  /// # Input options
  /// - <name> e.g. "sheath_boundary_penalty"
  ///   - penalty_timescale    Timescale in seconds
  ///   - sheath_ion_polytropic       Ion polytropic coefficient in Bohm sound speed
  ///   - sheath_electron_polytropic  Electron polytropic coefficient in Bohm sound speed
  ///
  /// # Reads from the mesh
  /// - penalty_mask[x,y,z]    A 3D field defining the shape of the boundary
  ///                          Equal 0 inside the plasma, 1 in the wall
  SheathBoundaryPenalty(std::string name, Options& options, Solver*);

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

  // The following functions are public for unit testing

  /// Cell indices where mask > threshold in the interior region
  static Region<Ind3D> buildPenaltyRegion(const Field3D& mask, BoutReal threshold = 1e-5);

  /// Apply boundary conditions to the mask and calculate the penalty region
  static PenaltyMaskData preparePenaltyMask(Field3D mask, BoutReal threshold = 1e-5);

  /// Calculate the field-aligned penalty region and extend the mask into Y guard cells
  static PenaltyMaskData prepareFieldAlignedPenaltyMask(Field3D mask_fa, Mesh& mesh,
                                                        BoutReal threshold = 1e-5);

  /// Calculate the volumetric penalty terms for a species
  static PenaltySourceData calculateVolumetricPenalty(
      const PenaltyMaskData& penalty_data, const Field3D& Ni, const Field3D& Ti,
      const Field3D& Vi, BoutReal Mi, BoutReal gamma_i, BoutReal penalty_timescale,
      const Field3D& density_source, const Field3D& momentum_source,
      const Field3D& energy_source, BoutReal density_floor = 1e-5);

  /// Calculate electron surface momentum penalty terms in field-aligned coordinates
  static Field3D calculateElectronSurfaceMomentumPenalty(
      const PenaltyMaskData& penalty_fa_data, const Field3D& Ne_fa, const Field3D& Te_fa,
      const Field3D& Ve_fa, const Field3D& phi_fa, BoutReal Me,
      BoutReal penalty_timescale, BoutReal density_floor = 1e-5,
      BoutReal mask_threshold = 1e-5);

  /// Calculate ion surface momentum penalty terms in field-aligned coordinates
  static Field3D calculateIonSurfaceMomentumPenalty(
      const PenaltyMaskData& penalty_fa_data, const Field3D& Ni_fa, const Field3D& Ti_fa,
      const Field3D& Te_fa, const Field3D& Vi_fa, BoutReal Mi,
      BoutReal sheath_ion_polytropic, BoutReal sheath_electron_polytropic,
      BoutReal penalty_timescale, BoutReal density_floor = 1e-5,
      BoutReal mask_threshold = 1e-5);

private:
  /// Mask function and region in the original mesh coordinates
  PenaltyMaskData penalty_data;

  BoutReal gamma_e, gamma_i;      // Sheath heat transmission factors
  BoutReal sheath_ion_polytropic; ///< Ion polytropic coefficient in Bohm sound speed
  BoutReal
      sheath_electron_polytropic; ///< Electron polytropic coefficient in Bohm sound speed

  /// Timescale of penalisation [normalised]
  BoutReal penalty_timescale;

  /// Diagnostics for output
  Options diagnostics;
  bool diagnose; ///< Save diagnostics?

  // Field-aligned fields for surface terms
  bool surface_terms;
  PenaltyMaskData penalty_fa_data;

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
  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<SheathBoundaryPenalty>
    registercomponentsheathboundarypenalty("sheath_boundary_penalty");
}

#endif // SHEATH_BOUNDARY_PENALTY_H
