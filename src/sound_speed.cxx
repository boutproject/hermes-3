
#include "../include/sound_speed.hxx"

void SoundSpeed::transform(Options &state) {
  Field3D total_pressure = 0.0;
  Field3D total_density = 0.0;

  Field3D fastest_wave = 0.0;
  for (auto& kv : state["species"].getChildren()) {
    const Options& species = kv.second;
    
    if (species.isSet("pressure")) {
      total_pressure += GET_NOBOUNDARY(Field3D, species["pressure"]);
    }

    if (species.isSet("density") and species.isSet("AA")) {
      total_density += GET_NOBOUNDARY(Field3D, species["density"]) * get<BoutReal>(species["AA"]);

      if (species.isSet("pressure")) {
        // Both pressure and mass density available
        auto P = GET_NOBOUNDARY(Field3D, species["pressure"]);
        auto N = GET_NOBOUNDARY(Field3D, species["density"]);
        auto AA = get<BoutReal>(species["AA"]);
        for (auto& i : fastest_wave.getRegion("RGN_NOBNDRY")) {
          BoutReal sound_speed = sqrt(P[i] / (N[i] * AA));
          fastest_wave[i] = BOUTMAX(fastest_wave[i], sound_speed);
        }
      }
    }
  }

  Field3D sound_speed = sqrt(total_pressure / floor(total_density, 1e-10));
  set(state["sound_speed"], sound_speed);
  set(state["fastest_wave"], fastest_wave);
}
