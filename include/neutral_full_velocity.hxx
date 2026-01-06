
#pragma once
#ifndef NEUTRAL_FULL_VELOCITY_H
#define NEUTRAL_FULL_VELOCITY_H

#include <string>

#include "component.hxx"

#include <bout/vector2d.hxx>

/// Neutral gas model, evolving three components of velocity as axisymmetric fields
///
/// Evolves neutral density, pressure and velocity as Field2D quantities
struct NeutralFullVelocity : public Component {
  NeutralFullVelocity(const std::string& name, Options& options, Solver* solver);

  /// Modify the given simulation state
  void transform(Options& state) override;

  /// Use the final simulation state to update internal state
  /// (e.g. time derivatives)
  void finally(const Options& state) override;

  /// Add extra fields for output, or set attributes e.g docstrings
  void outputVars(Options& state) override;

private:
  Coordinates* coord; // Coordinate system

  std::string name; // Name of this species
  BoutReal AA;      // Atomic mass

  BoutReal density_floor; ///< Floor when dividing by density
  BoutReal temperature_floor;
  BoutReal pressure_floor; ///< Minimum Pn used when dividing Pn by Nn to get Tn.

  Field2D Nn2D;                // Neutral gas density (evolving)
  Field2D Pn2D;                // Neutral gas pressure (evolving)
  Vector2D Vn2D;               // Neutral gas velocity
  Vector2D Vn2D_contravariant; ///< Neutral gas velocity v^x, v^y, v^z
  Field2D Tn2D;

  // Transformation to cylindrical coordinates
  Field2D Rxy;
  BoutReal sigma_Bp; // Sign of poloidal field

  // Grad x = Txr * Grad R + Txz * Grad Z
  // Grad y = Tyr * Grad R + Tyz * Grad Z
  Field2D Txr, Txz;
  Field2D Tyr, Tyz;

  // Grad R = Urx * Grad x + Ury * Grad y
  // Grad Z = Uzx * Grad x + Uzy * Grad y
  Field2D Urx, Ury;
  Field2D Uzx, Uzy;

  BoutReal adiabatic_index;    // Ratio of specific heats
  BoutReal neutral_viscosity;  // Neutral gas viscosity
  BoutReal neutral_conduction; // Neutral gas thermal conduction
  BoutReal neutral_gamma;      // Heat transmission for neutrals


  std::vector<std::string> collision_names; ///< Collisions used for collisionality
  std::string diffusion_collisions_mode;  ///< Collision selection, either afn or multispecies
  Field2D nu; ///< Collisionality to use for diffusion

  Field2D Dnn; ///< Diffusion coefficient
  Field2D kappa_n, eta_n; ///< Neutral conduction and viscosity

  BoutReal flux_limit; ///< Diffusive flux limit
  BoutReal neutral_lmax;
  BoutReal diffusion_limit;    ///< Maximum diffusion coefficient

  // Toroidal advection
  bool toroidal_flow;      ///< Evolve toroidal flow?
  bool momentum_advection; ///< Include advection of momentum?
  bool curved_torus;       ///< Include toroidal curvature in momentum advection?
  bool constant_transport_coef; ///< Use constant transport coefficients?

  bool zero_timederivs; ///< Set the time derivatives to zero?
  bool output_ddt; ///< Save time derivatives?  

  bool diagnose; ///< Output additional diagnostics?
  Field2D Vnpar; ///< Parallel flow velocity diagnostic

  Field2D density_source, pressure_source; ///< External input source
  Field2D Sn, Sp, Snv; ///< Particle, pressure and momentum source

};

namespace {
RegisterComponent<NeutralFullVelocity>
    registersolverneutralfullvelocity("neutral_full_velocity");
}

#endif // NEUTRAL_FULL_VELOCITY_H
