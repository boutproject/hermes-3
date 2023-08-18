#pragma once
#ifndef COLLISIONS_H
#define COLLISIONS_H

#include <bout/field3d.hxx>

#include "component.hxx"

/// Calculates the collision rate of each species
/// with all other species
/// 
/// Important: Be careful when including both ion_neutral collisions
///            and reactions such as charge exchange, since that may
///            result in double counting. Similarly for
///            electron_neutral collisions and ionization reactions.
///
struct Collisions : public Component {
  ///
  /// @param alloptions Settings, which should include:
  ///    - units
  ///      - eV
  ///      - inv_meters_cubed
  ///      - meters
  ///      - seconds
  ///
  /// The following boolean options under alloptions[name] control
  /// which collisions are calculated:
  ///
  ///   - electron_electron
  ///   - electron_ion
  ///   - electron_neutral
  ///   - ion_ion
  ///   - ion_neutral
  ///   - neutral_neutral
  ///
  /// There are also switches for other terms:
  ///
  ///   - frictional_heating    Include R dot v heating term as energy source? (includes Ohmic heating)
  ///
  Collisions(std::string name, Options& alloptions, Solver*);

  void transform(Options &state) override;

  /// Add extra fields for output, or set attributes e.g docstrings
  void outputVars(Options &state) override;

private:
  BoutReal Tnorm; // Temperature normalisation [eV]
  BoutReal Nnorm; // Density normalisation [m^-3]
  BoutReal rho_s0;  // Length normalisation [m]
  BoutReal Omega_ci; // Frequency normalisation [s^-1]

  /// Which types of collisions to include?
  bool electron_electron, electron_ion, electron_neutral, ion_ion, ion_neutral,
      neutral_neutral;

  /// Include frictional heating term?
  bool frictional_heating;

  /// Calculated collision rates saved for post-processing and use by other components
  /// Saved in options, the BOUT++ dictionary-like object
  Options collision_rates;

  /// Save more diagnostics?
  bool diagnose;

  /// Update collision frequencies, momentum and energy exchange
  /// nu_12    normalised frequency
  /// momentum_coefficient   Leading coefficient on parallel friction
  ///                        e.g 0.51 for electron-ion with Zi=1
  void collide(Options &species1, Options &species2, const Field3D &nu_12, BoutReal momentum_coefficient);
};

namespace {
RegisterComponent<Collisions> registercomponentcollisions("collisions");
}

#endif // COLLISIONS_H
