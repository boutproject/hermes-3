#pragma once
#ifndef SOUND_SPEED_H
#define SOUND_SPEED_H

#include "component.hxx"

/// Calculate the system sound speed
///
/// This uses the sum of all species pressures and mass densities
/// so should run after those have been set.
struct SoundSpeed : public Component {
  SoundSpeed(std::string, Options&, Solver*) {}
  
  /// This sets in the state
  /// - sound_speed     The collective sound speed, based on total pressure and total mass density
  /// - fastest_wave    The highest species sound speed at each point in the domain
  ///
  /// Optional inputs:
  /// - species
  ///   - ...    // Iterates over all species
  ///     - density
  ///     - AA       // Atomic mass
  ///     - pressure
  ///
  void transform(Options &state) override;
};

namespace {
RegisterComponent<SoundSpeed> registercomponentsoundspeed("sound_speed");
}

#endif // SOUND_SPEED_H
