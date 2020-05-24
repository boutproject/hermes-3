#pragma once
#ifndef COLLISIONS_H
#define COLLISIONS_H

#include <field3d.hxx>

#include "component.hxx"

/// Calculates the collision rate of each species
/// with all other species
/// 
/// 
struct Collisions : public Component {
  Collisions(std::string name, Options& alloptions, Solver*);

  void transform(Options &state) override;

private:
  BoutReal Tnorm; // Temperature normalisation [eV]
  BoutReal Nnorm; // Density normalisation [m^-3]
  BoutReal rho_s0;  // Length normalisation [m]
  BoutReal Omega_ci; // Frequency normalisation [s^-1]

  /// Which types of collisions to include?
  bool electron_electron, electron_ion, electron_neutral, ion_ion, ion_neutral,
      neutral_neutral;

  /// Update collision frequencies, momentum and energy exchange
  /// nu_12    normalised frequency
  void collide(Options &species1, Options &species2, const Field3D &nu_12);
};



#endif // COLLISIONS_H
