
#include "../include/sound_speed.hxx"
#include "../include/component.hxx"
#include "../include/guarded_options.hxx"
#include "../include/hermes_utils.hxx"

#include <bout/bout_types.hxx>
#include <bout/field.hxx>
#include <bout/field3d.hxx>
#include <bout/mesh.hxx>
#include <bout/utils.hxx>

#include <cmath>

void SoundSpeed::transform_impl(GuardedOptions& state) {
  Field3D total_pressure = 0.0;
  Field3D total_density = 0.0;

  Field3D fastest_wave = 0.0;
  for (auto& kv : state["species"].getChildren()) {
    const GuardedOptions species = kv.second;

    if (species.isSet("pressure")) {
      total_pressure += GET_NOBOUNDARY(Field3D, species["pressure"]);
    }

    if ((kv.first == "e") and !electron_dynamics) {
      // Exclude electron sound speed, but include electron pressure in
      // collective sound speed calculation (total_pressure).
      continue;
    }

    if (species.isSet("AA")) {
      auto AA = get<BoutReal>(species["AA"]); // Atomic mass number

      if (species.isSet("density")) {
        total_density +=
            GET_NOBOUNDARY(Field3D, species["density"]) * get<BoutReal>(species["AA"]);
      }

      if (species.isSet("temperature")) {
        auto T = GET_NOBOUNDARY(Field3D, species["temperature"]);
        for (const auto& i : fastest_wave.getRegion("RGN_NOBNDRY")) {
          const BoutReal sound_speed = sqrt(softFloor(T[i], temperature_floor) / AA);
          fastest_wave[i] = BOUTMAX(fastest_wave[i], sound_speed);
        }
      }
    }
  }

  total_density = softFloor(total_density, 1e-10);
  Field3D sound_speed = sqrt(total_pressure / total_density);
  for (const auto& i : fastest_wave.getRegion("RGN_NOBNDRY")) {
    fastest_wave[i] = BOUTMAX(fastest_wave[i], sound_speed[i]);
  }

  if (alfven_wave) {
    auto* coord = fastest_wave.getCoordinates();
    for (const auto& i : fastest_wave.getRegion("RGN_NOBNDRY")) {
      const BoutReal alfven_speed = beta_norm * coord->Bxy[i] / sqrt(total_density[i]);
      fastest_wave[i] = BOUTMAX(fastest_wave[i], alfven_speed);
    }
  }

  // Communicate guard cells
  fastest_wave.getMesh()->communicate(fastest_wave, sound_speed);
  fastest_wave.resetRegion();
  sound_speed.resetRegion();

  set(state["sound_speed"], sound_speed);
  set(state["fastest_wave"], fastest_wave * fastest_wave_factor);
}
