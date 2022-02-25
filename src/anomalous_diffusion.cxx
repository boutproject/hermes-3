#include "../include/anomalous_diffusion.hxx"

#include <bout/fv_ops.hxx>
#include "../include/div_ops.hxx"

AnomalousDiffusion::AnomalousDiffusion(std::string name, Options &alloptions, Solver *) : name(name) {
  // Normalisations
  const Options& units = alloptions["units"];
  const BoutReal rho_s0 = units["meters"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();

  const BoutReal diffusion_norm = rho_s0 * rho_s0 * Omega_ci; // m^2/s
  
  Options& options = alloptions[name];
  
  include_D = options.isSet("anomalous_D");
  if (include_D) {
    anomalous_D = options["anomalous_D"]
                      .doc("Anomalous particle diffusion coefficient [m^2/s]")
                      .as<Field2D>()
                  / diffusion_norm;
  }

  include_chi = options.isSet("anomalous_chi");
  if (include_chi) {
    anomalous_chi = options["anomalous_chi"]
                        .doc("Anomalous thermal diffusion coefficient [m^2/s]")
                        .as<Field2D>()
                    / diffusion_norm;
  }

  include_nu = options.isSet("anomalous_nu");
  if (include_nu) {
    anomalous_nu = options["anomalous_nu"]
                       .doc("Anomalous momentum diffusion coefficient [m^2/s]")
                       .as<Field2D>()
                   / diffusion_norm;
  }
}

void AnomalousDiffusion::transform(Options &state) {

  Options& species = state["species"][name];

  // Diffusion operates on 2D (axisymmetric) profiles
  const Field3D N = GET_NOBOUNDARY(Field3D, species["density"]);
  const Field2D N2D = DC(N);

  const Field3D T = species.isSet("temperature")
                        ? GET_NOBOUNDARY(Field3D, species["temperature"])
                        : 0.0;
  const Field2D T2D = DC(T);

  const Field3D V =
      species.isSet("velocity") ? GET_NOBOUNDARY(Field3D, species["velocity"]) : 0.0;
  const Field2D V2D = DC(V);
  
  if (include_D) {
    // Particle diffusion. Gradients of density drive flows of particles,
    // momentum and energy
    add(species["density_source"],
        Div_a_Laplace_perp_upwind(anomalous_D, N2D));

    // Note: Upwind operators used, or unphysical increases
    // in temperature and flow can be produced
    add(species["momentum_source"],
        Div_a_Laplace_perp_upwind(V2D * anomalous_D, N2D));

    add(species["energy_source"],
        Div_a_Laplace_perp_upwind((3./2)*T2D * anomalous_D, N2D));
  }

  if (include_chi) {
    // Gradients in temperature which drive energy flows
    add(species["energy_source"],
        Div_a_Laplace_perp_upwind(anomalous_chi * N2D, T2D));
  }

  if (include_nu) {
    // Gradients in slow speed which drive momentum flows
    add(species["momentum_source"],
        Div_a_Laplace_perp_upwind(anomalous_nu * N2D, V2D));
  }
}
