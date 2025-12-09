#pragma once
#ifndef SOUND_SPEED_H
#define SOUND_SPEED_H

#include "component.hxx"

#include <bout/constants.hxx>

/// Calculate the system sound speed
///
/// This uses the sum of all species pressures and mass densities
/// so should run after those have been set.
struct SoundSpeed : public Component {
  SoundSpeed(std::string name, Options& alloptions, Solver*)
      : Component({readOnly("species:{all_species}:pressure", Regions::Interior),
                   writeFinal("sound_speed"), writeFinal("fastest_wave")}) {
    Options &options = alloptions[name];
    electron_dynamics = options["electron_dynamics"]
      .doc("Include electron sound speed?")
      .withDefault<bool>(true);

    alfven_wave = options["alfven_wave"]
      .doc("Include Alfven wave speed?")
      .withDefault<bool>(false);
    if (alfven_wave) {
      // Calculate normalisation factor
      const auto& units = alloptions["units"];
      const auto Bnorm = get<BoutReal>(units["Tesla"]);
      const auto Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
      const auto Cs0 = get<BoutReal>(units["meters"]) / get<BoutReal>(units["seconds"]);
      beta_norm = Bnorm / sqrt(SI::mu0 * Nnorm * SI::Mp) / Cs0;
    }

    temperature_floor = options["temperature_floor"]
      .doc("Minimum temperature when calculating sound speeds [eV]")
      .withDefault(0.0);

    fastest_wave_factor = options["fastest_wave_factor"]
      .doc("Multiply the fastest wave by this factor, affecting lax flux strength")
      .withDefault(1.0);

    if (temperature_floor > 0.0) {
      temperature_floor /= get<BoutReal>(alloptions["units"]["eV"]);
    }

    if (electron_dynamics) {
      setPermissions(readIfSet("species:{all_species}:AA"));
      // FIXME: Only read if AA is set
      setPermissions(readIfSet("species:{all_species}:{opt_inputs}", Regions::Interior));
    } else {
      setPermissions(readIfSet("species:{non_electrons}:AA"));
      // FIXME: Only read if AA is set
      setPermissions(
          readIfSet("species:{all_species}:{non_electrons}", Regions::Interior));
    }
    substitutePermissions("opt_inputs", {"density", "temperature"});
  }

private:
  bool electron_dynamics; ///< Include electron sound speed?
  bool alfven_wave; ///< Include Alfven wave speed?
  BoutReal beta_norm{0.0}; ///< Normalisation factor for Alfven speed
  BoutReal temperature_floor; ///< Minimum temperature when calculating speed
  BoutReal fastest_wave_factor; ///< Multiply the fastest wave by this factor

  /// This sets in the state
  /// - sound_speed     The collective sound speed, based on total pressure and total mass
  /// density
  /// - fastest_wave    The highest species sound speed at each point in the domain
  ///
  /// Optional inputs:
  /// - species
  ///   - ...    // Iterates over all species
  ///     - density
  ///     - AA       // Atomic mass
  ///     - pressure
  ///     - temperature
  ///
  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<SoundSpeed> registercomponentsoundspeed("sound_speed");
}

#endif // SOUND_SPEED_H
