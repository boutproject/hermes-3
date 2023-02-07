
#pragma once
#ifndef FULL_VELOCITY_H
#define FULL_VELOCITY_H

#include <string>

#include "component.hxx"

#include "vector2d.hxx"

/// Neutral gas model, evolving three components of velocity as axisymmetric fields
///
/// Evolves neutral density, pressure and velocity
struct NeutralFullVelocity : public Component {
  NeutralFullVelocity(const std::string& name, Options& options, Solver *solver);
  
  /// Modify the given simulation state
  void transform(Options &state) override;
  
  /// Use the final simulation state to update internal state
  /// (e.g. time derivatives)
  void finally(const Options &state) override;

  /// Add extra fields for output, or set attributes e.g docstrings
  void outputVars(Options &state) override;
private:
  Coordinates *coord;   // Coordinate system

  std::string name; // Name of this species
  BoutReal AA; // Atomic mass

  BoutReal Tnorm;
  
  Field2D Nn2D;         // Neutral gas density (evolving)
  Field2D Pn2D;         // Neutral gas pressure (evolving)
  Vector2D Vn2D;        // Neutral gas velocity
  Field2D Tn2D;
  Field2D DivV2D;       // Divergence of gas velocity

  // Transformation to cylindrical coordinates
  
  // Grad x = Txr * Grad R + Txz * Grad Z
  // Grad y = Tyr * Grad R + Tyz * Grad Z
  Field2D Txr, Txz;
  Field2D Tyr, Tyz;
  
  // Grad R = Urx * Grad x + Ury * Grad y
  // Grad Z = Uzx * Grad x + Uzy * Grad y
  Field2D Urx, Ury;
  Field2D Uzx, Uzy;

  BoutReal gamma_ratio;    // Ratio of specific heats
  BoutReal neutral_viscosity; // Neutral gas viscosity
  BoutReal neutral_bulk;   // Neutral gas bulk viscosity
  BoutReal neutral_conduction; // Neutral gas thermal conduction
  BoutReal neutral_gamma; // Heat transmission for neutrals
  
  // Outflowing boundaries for neutrals
  bool outflow_ydown; // Allow neutral outflows?
};

namespace {
RegisterComponent<NeutralFullVelocity> registersolverfullvelocity("full_velocity");
}

#endif // FULL_VELOCITY_H
