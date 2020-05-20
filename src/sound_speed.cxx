
#include "../include/sound_speed.hxx"

void SoundSpeed::transform(Options &state) {
  Field3D total_pressure = 0.0;
  Field3D total_density = 0.0;
  
  for (auto& kv : state["species"].getChildren()) {
    const Options& species = kv.second;
    
    if (species.isSet("pressure")) {
      total_pressure += get<Field3D>(species["pressure"]);
    }

    if (species.isSet("density") and species.isSet("AA")) {
      total_density += get<Field3D>(species["density"]) * get<BoutReal>(species["AA"]);
    }
  }

  set(state["sound_speed"], sqrt(total_pressure / total_density));
}
