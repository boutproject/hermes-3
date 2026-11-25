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
struct BraginskiiCollisions : public NamedComponent<BraginskiiCollisions> {
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
  ///   - ion_neutral_sigma          ion-neutral cross-section [m^2], default 5e-19
  ///   - ion_neutral_vrel_reduced   use the mean relative speed of two Maxwellians
  ///                                in the ion-neutral rate, default false
  ///
  BraginskiiCollisions(const std::string& name, Options& alloptions, Solver*);

  /// Add extra fields for output, or set attributes e.g docstrings
  void outputVars(Options& state) override;

  static constexpr auto type = "braginskii_collisions";

private:
  BoutReal Tnorm;    // Temperature normalisation [eV]
  BoutReal Nnorm;    // Density normalisation [m^-3]
  BoutReal rho_s0;   // Length normalisation [m]
  BoutReal Omega_ci; // Frequency normalisation [s^-1]

  /// Which types of collisions to include?
  bool electron_electron, electron_ion, electron_neutral, ion_ion, ion_neutral,
      neutral_neutral;

  BoutReal ei_multiplier; // Arbitrary user-set multiplier on electron-ion collisions

  /// Ion-neutral hard-sphere rate: cross-section [m^2] (default 5e-19, the value
  /// used before) and choice of the relative speed. If ion_neutral_vrel_reduced is
  /// false (default), v_rel = sqrt(T1/m1 + T2/m2) as before; if true, the mean
  /// relative speed of two Maxwellians, v_rel = sqrt(8/pi (T1/m1 + T2/m2))
  /// = sqrt(8 T / (pi mu)) with mu = m1 m2 / (m1 + m2), T = (m2 T1 + m1 T2) / (m1 + m2).
  BoutReal ion_neutral_sigma;
  bool ion_neutral_vrel_reduced;

  /// Calculated collision rates saved for post-processing and use by other components
  /// Saved in options, the BOUT++ dictionary-like object
  Options collision_rates;

  /// Save more diagnostics?
  bool diagnose;

  void transform_impl(GuardedOptions& state) override;

  /// Update collision frequencies, momentum and energy exchange
  /// nu_12    normalised frequency
  void collide(GuardedOptions& species1, GuardedOptions& species2, const Field3D& nu_12);
};

namespace {
RegisterComponent<BraginskiiCollisions> registercomponentbraginskiicollisions;
}

#endif // COLLISIONS_H
