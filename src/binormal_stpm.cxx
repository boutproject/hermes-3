
#include <bout/fv_ops.hxx>
#include <bout/output_bout_types.hxx>
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
    auto AA = get<BoutReal>(species["AA"]);

    const Field3D N = species.isSet("density")
      ? GET_NOBOUNDARY(Field3D, species["density"])
      : 0.0;
    const Field3D T = species.isSet("temperature")
      ? GET_NOBOUNDARY(Field3D, species["temperature"])
      : 0.0;
    const Field3D NV = species.isSet("momentum")
      ? GET_NOBOUNDARY(Field3D, species["momentum"])
      : 0.0;
    
    add(species["pressure_source"],
	(2. / 3) * (1/Theta) * FV::Div_par_K_Grad_par(chi_Theta*N, T, false));

    add(species["momentum_source"],
	(1/Theta) * FV::Div_par_K_Grad_par(AA*nu_Theta, NV, false));
    
    add(species["density_source"],
	(1/Theta) * FV::Div_par_K_Grad_par(D_Theta, N, false));

  }
}

