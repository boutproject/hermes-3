
#include <bout/fv_ops.hxx>
#include <bout/difops.hxx>
#include <bout/output_bout_types.hxx>
#include <bout/initialprofiles.hxx>
#include "../include/div_ops.hxx"
#include "../include/binormal_stpm.hxx"
#include "../include/hermes_utils.hxx"
#include "../include/hermes_build_config.hxx"

using bout::globals::mesh;

BinormalSTPM::BinormalSTPM(std::string name, Options& alloptions, Solver* solver)
  : name(name) {
  AUTO_TRACE();
  auto& options = alloptions[name];
  const Options& units = alloptions["units"];
  const BoutReal rho_s0 = units["meters"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();

  const BoutReal diffusion_norm = rho_s0 * rho_s0 * Omega_ci; // m^2/s
  
  Theta = options["Theta"]
    .doc("Field-line Pitch defined by Feng et al.")
    .withDefault(1e-3);

  chi = options["chi"]
    .doc("Anomalous heat diffusion.")
    .withDefault(3.)
    /diffusion_norm;

  D = options["D"]
    .doc("Anomalous density diffusion.")
    .withDefault(1.)
    /diffusion_norm;
  
  nu = options["nu"]
    .doc("Anomalous momentum diffusion.")
    .withDefault(1.)
    /diffusion_norm;

  Vbn = options["Vbn"]
    .doc("Binormal velocity.")
    .withDefault(0.);

  chi_Theta = chi/Theta;
  D_Theta = D/Theta;
  nu_Theta = nu/Theta;
  
}

void BinormalSTPM::transform(Options& state) {
  AUTO_TRACE();
  Options& allspecies = state["species"];
  // Loop through all species
  for (auto& kv : allspecies.getChildren()) {
    const auto& species_name = kv.first;

    Options& species = allspecies[species_name];

    const Field3D T = get<Field3D>(species["temperature"]);
    const Field3D NV = get<Field3D>(species["momentum"]);
    const Field3D N = get<Field3D>(species["density"]);
    
    add(species["pressure_source"], (1/Theta) * FV::Div_par_K_Grad_par(chi_Theta, T, false));
    add(species["pressure_source"], (3/(2*Theta)) * Grad_par(N*T*Vbn/Theta));    

    add(species["momentum_source"], FV::Div_par_K_Grad_par(nu_Theta, NV, false));
    
    add(species["density_source"], (1/Theta) * FV::Div_par_K_Grad_par(D_Theta, N, false));
    add(species["density_source"], (1/Theta) * Grad_par(N*Vbn/Theta));
  }
}



