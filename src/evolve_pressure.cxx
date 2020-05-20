
#include <bout/fv_ops.hxx>
#include <difops.hxx>

#include "../include/evolve_pressure.hxx"
#include "../include/div_ops.hxx"

EvolvePressure::EvolvePressure(std::string name, Options &alloptions, Solver *solver) : name(name) {
  AUTO_TRACE();
  
  solver->add(P, std::string("P") + name);

  auto& options = alloptions[name];
  
  bndry_flux = options["bndry_flux"]
                   .doc("Allow flows through radial boundaries")
                   .withDefault<bool>(true);

  poloidal_flows = options["poloidal_flows"]
                       .doc("Include poloidal ExB flow")
                       .withDefault<bool>(true);

  thermal_conduction = options["thermal_conduction"]
                           .doc("Include parallel heat conduction?")
                           .withDefault<bool>(true);
}

void EvolvePressure::transform(Options &state) {
  AUTO_TRACE();
  
  auto& species = state["species"][name];

  // Check for low pressures, ensure that Pi >= 0
  P = floor(P, 0.0);
  
  set(species["pressure"], P);

  // Calculate temperature
  N = get<Field3D>(species["density"]);
  T = P / N;
  set(species["temperature"], T);
}

void EvolvePressure::finally(const Options &state) {
  AUTO_TRACE();

  /// Get the section containing this species
  const auto& species = state["species"][name];
  
  if (state.isSection("fields") and state["fields"].isSet("phi")) {
    // Electrostatic potential set -> include ExB flow

    Field3D phi = get<Field3D>(state["fields"]["phi"]);

    ddt(P) = -Div_n_bxGrad_f_B_XPPM(P, phi, bndry_flux, poloidal_flows, true);
  } else {
    ddt(P) = 0.0;
  }

  if (species.isSet("velocity")) {
    Field3D V = get<Field3D>(species["velocity"]);
    
    // Typical wave speed used for numerical diffusion
    Field3D sound_speed;
    if (state.isSet("sound_speed")) {
      Field3D sound_speed = get<Field3D>(state["sound_speed"]);
    } else {
      Field3D T = get<Field3D>(species["temperature"]);
      sound_speed = sqrt(T);
    }
    
    ddt(P) -= FV::Div_par(P, V, sound_speed);

    // Work done. This balances energetically a term in the momentum equation
    ddt(P) -= (2. / 3) * P * Div_par(V);
  }

  // Parallel heat conduction
  if (thermal_conduction) {
    
    // Calculate ion collision times
    Field3D tau = 1. / get<Field3D>(species["collision_rate"]);
    
    // Parallel heat conduction
    Field3D kappa_par = 3.9 * P * tau;
    
    ddt(P) += (2. / 3) * FV::Div_par_K_Grad_par(kappa_par, T);
  }
  
  //////////////////////
  // Other sources

  if (species.isSet("energy_source")) {
    ddt(P) += (2./3) * get<Field3D>(species["energy_source"]);
  }
}
