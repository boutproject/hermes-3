#pragma once
#ifndef COLLISIONS_H
#define COLLISIONS_H

#include <string>

#include <bout/bout_types.hxx>
#include <bout/field3d.hxx>
#include <bout/options.hxx>

#include "component.hxx"

/// Calculates the collision rate of each species
/// with all other species, using the Braginskii equation.
///
/// Important: Be careful when including both ion_neutral collisions
///            and reactions such as charge exchange, since that may
///            result in double counting. Similarly for
///            electron_neutral collisions and ionization reactions.
///
struct BraginskiiCollisions : public Component {
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
  /// The user can also set
  ///
  ///   - ei_multiplier   arbitrary multiplier on electron-ion collision rate
  ///
  BraginskiiCollisions(const std::string& name, Options& alloptions, Solver*);

  void transform(Options& state) override;

  /// Add extra fields for output, or set attributes e.g docstrings
  void outputVars(Options& state) override;

private:
  BoutReal Tnorm;    // Temperature normalisation [eV]
  BoutReal Nnorm;    // Density normalisation [m^-3]
  BoutReal rho_s0;   // Length normalisation [m]
  BoutReal Omega_ci; // Frequency normalisation [s^-1]

  /// Which types of collisions to include?
  bool electron_electron, electron_ion, electron_neutral, ion_ion, ion_neutral,
      neutral_neutral;

  BoutReal ei_multiplier; // Arbitrary user-set multiplier on electron-ion collisions
  BoutReal density_floor;

  /// Calculated collision rates saved for post-processing and use by other components
  /// Saved in options, the BOUT++ dictionary-like object
  Options collision_rates;

  /// Save more diagnostics?
  bool diagnose;

  /// Update collision frequencies
  /// nu_12    normalised frequency
  void collide(Options& species1, Options& species2, const Field3D& nu_12);
};

namespace {
RegisterComponent<BraginskiiCollisions>
    registercomponentbraginskiicollisions("braginskii_collisions");
}

#endif // COLLISIONS_H
