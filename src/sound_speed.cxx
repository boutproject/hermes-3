
#include "../include/sound_speed.hxx"

namespace {
BoutReal floor(BoutReal value, BoutReal min) {
  if (value < min)
    return min;
  return value;
}
} // namespace

void SoundSpeed::transform(Options &state) {
  Field3D total_pressure = 0.0;
  Field3D total_density = 0.0;

  Field3D fastest_wave = 0.0;
  for (auto& kv : state["species"].getChildren()) {
    const Options& species = kv.second;

    if (species.isSet("pressure")) {
      total_pressure += GET_NOBOUNDARY(Field3D, species["pressure"]);
    }

    if ((kv.first == "e") and !electron_dynamics) {
      // Exclude electron sound speed, but include electron pressure in
      // collective sound speed calculation (total_pressure).
      continue;
    }

    if (species.isSet("density") and species.isSet("AA")) {
      total_density += GET_NOBOUNDARY(Field3D, species["density"]) * get<BoutReal>(species["AA"]);

      if (species.isSet("pressure")) {
        // Both pressure and mass density available
        auto P = GET_NOBOUNDARY(Field3D, species["pressure"]);
        auto N = GET_NOBOUNDARY(Field3D, species["density"]);
        auto AA = get<BoutReal>(species["AA"]);
        for (auto& i : fastest_wave.getRegion("RGN_NOBNDRY")) {
          BoutReal sound_speed = sqrt(P[i] / (floor(N[i], 1e-5) * AA));
          fastest_wave[i] = BOUTMAX(fastest_wave[i], sound_speed);
        }
      }
    }
  }

  Field3D sound_speed = sqrt(total_pressure / floor(total_density, 1e-10));
  for (auto& i : fastest_wave.getRegion("RGN_NOBNDRY")) {
    fastest_wave[i] = BOUTMAX(fastest_wave[i], sound_speed[i]);
  }
  set(state["sound_speed"], sound_speed);
  set(state["fastest_wave"], fastest_wave);
}
