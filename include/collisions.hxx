#pragma once
#ifndef COLLISIONS_H
#define COLLISIONS_H

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
};



#endif // COLLISIONS_H
