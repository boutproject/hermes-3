#pragma once
#ifndef ION_VISCOSITY_H
#define ION_VISCOSITY_H

#include "component.hxx"

/// Ion viscosity terms
///
/// Adds a viscosity to all species which are not electrons
///
/// Uses Braginskii collisional form, combined with a
/// SOLPS-like flux limiter.
///
/// Needs to be calculated after collisions, because collision
/// frequency is used to calculate parallel viscosity
///
/// The ion stress tensor Pi_ci is split into perpendicular and
/// parallel pieces:
///
///    Pi_ci = Pi_ciperp + Pi_cipar
///
/// In the parallel ion momentum equation the Pi_cipar term
/// is solved as a parallel diffusion, so is treated separately
/// All other terms are added to Pi_ciperp, even if they are
/// not really parallel parts
struct IonViscosity : public Component {
  /// Inputs
  /// - <name>
  ///   - eta_limit_alpha: float, default -1
  ///         Flux limiter coefficient. < 0 means off.
  ///   - perpendicular: bool, default false
  ///         Include perpendicular flows?
  ///         Requires curvature vector and phi potential
  ///
  IonViscosity(std::string name, Options& alloptions, Solver*);

  /// Inputs
  /// - species
  ///   - <name>   (skips "e")
  ///     - pressure  (skips if not present)
  ///     - velocity  (skips if not present)
  ///     - collision_frequency
  ///
  /// Sets in the state
  /// - species
  ///   - <name>
  ///     - momentum_source
  ///
  void transform(Options &state) override;

  /// Save variables to the output
  void outputVars(Options &state) override;
private:
  BoutReal eta_limit_alpha; ///< Flux limit coefficient
  bool perpendicular; ///< Include perpendicular flow? (Requires phi)
  std::vector<std::string> collision_names; ///< Collisions used for collisionality
  std::string viscosity_collisions_mode;  ///< Collision selection, either legacy or braginskii
  Field3D nu;   ///< Collision frequency for conduction
  Vector2D Curlb_B; ///< Curvature vector Curl(b/B)

  bool diagnose; ///< Output additional diagnostics?

  /// Per-species diagnostics
  struct Diagnostics {
    Field3D Pi_ciperp; ///< Perpendicular part of Pi scalar
    Field3D Pi_cipar;  ///< Parallel part of Pi scalar
    Field3D DivJ;      ///< Divergence of current in vorticity equation
  };

  /// Store diagnostics for each species
  std::map<std::string, Diagnostics> diagnostics;
};

namespace {
RegisterComponent<IonViscosity> registercomponentionviscosity("ion_viscosity");
}

#endif
