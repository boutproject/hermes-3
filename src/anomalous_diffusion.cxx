#include "../include/anomalous_diffusion.hxx"

#include "../include/div_ops.hxx"
#include <bout/fv_ops.hxx>
#include <bout/output_bout_types.hxx>

using bout::globals::mesh;

AnomalousDiffusion::AnomalousDiffusion(std::string name, Options& alloptions, Solver*)
    : name(name) {
  // Normalisations
  const Options& units = alloptions["units"];
  const BoutReal rho_s0 = units["meters"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();

  const BoutReal diffusion_norm = rho_s0 * rho_s0 * Omega_ci; // m^2/s

  Options& options = alloptions[name];

  // Set in the mesh or options (or both)
  anomalous_D = 0.0;
  include_D = (mesh->get(anomalous_D, std::string("D_") + name) == 0)
              || options.isSet("anomalous_D");
  // Option overrides mesh value
  anomalous_D = options["anomalous_D"]
                    .doc("Anomalous particle diffusion coefficient [m^2/s]")
                    .withDefault(anomalous_D)
                / diffusion_norm;

  anomalous_chi = 0.0;
  include_chi = (mesh->get(anomalous_chi, std::string("chi_") + name) == 0)
                || options.isSet("anomalous_chi");
  anomalous_chi = options["anomalous_chi"]
                      .doc("Anomalous thermal diffusion coefficient [m^2/s]")
                      .withDefault(anomalous_chi)
                  / diffusion_norm;

  anomalous_nu = 0.0;
  include_nu = (mesh->get(anomalous_nu, std::string("nu_") + name) == 0)
               || options.isSet("anomalous_nu");
  anomalous_nu = options["anomalous_nu"]
                     .doc("Anomalous momentum diffusion coefficient [m^2/s]")
                     .withDefault(anomalous_nu)
                 / diffusion_norm;

  anomalous_sheath_flux = options["anomalous_sheath_flux"]
                              .doc("Allow anomalous diffusion into sheath?")
                              .withDefault<bool>(false);

  diagnose = alloptions[name]["diagnose"]
                   .doc("Output additional diagnostics?")
                   .withDefault<bool>(false);
}

void AnomalousDiffusion::transform(Options& state) {
  AUTO_TRACE();

  Options& species = state["species"][name];

  // Diffusion operates on 2D (axisymmetric) profiles
  // Note: Includes diffusion in Y, so set boundary fluxes
  // to zero by imposing neumann boundary conditions.
  const Field3D N = GET_NOBOUNDARY(Field3D, species["density"]);
  Field2D N2D = DC(N);

  const Field3D T = species.isSet("temperature")
                        ? GET_NOBOUNDARY(Field3D, species["temperature"])
                        : 0.0;
  Field2D T2D = DC(T);

  const Field3D V =
      species.isSet("velocity") ? GET_NOBOUNDARY(Field3D, species["velocity"]) : 0.0;
  Field2D V2D = DC(V);

  if (!anomalous_sheath_flux) {
    // Apply Neumann Y boundary condition, so no additional flux into boundary
    // Note: Not setting radial (X) boundaries since those set radial fluxes
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        N2D(r.ind, mesh->ystart - 1, jz) = N2D(r.ind, mesh->ystart, jz);
        T2D(r.ind, mesh->ystart - 1, jz) = T2D(r.ind, mesh->ystart, jz);
        V2D(r.ind, mesh->ystart - 1, jz) = V2D(r.ind, mesh->ystart, jz);
      }
    }
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        N2D(r.ind, mesh->yend + 1, jz) = N2D(r.ind, mesh->yend, jz);
        T2D(r.ind, mesh->yend + 1, jz) = T2D(r.ind, mesh->yend, jz);
        V2D(r.ind, mesh->yend + 1, jz) = V2D(r.ind, mesh->yend, jz);
      }
    }
  }

  if (include_D) {
    // Particle diffusion. Gradients of density drive flows of particles,
    // momentum and energy
    add(species["density_source"], Div_a_Grad_perp_upwind(anomalous_D, N2D));

    // Note: Upwind operators used, or unphysical increases
    // in temperature and flow can be produced
    add(species["momentum_source"], Div_a_Grad_perp_upwind(V2D * anomalous_D, N2D));

    add(species["energy_source"],
        Div_a_Grad_perp_upwind((3. / 2) * T2D * anomalous_D, N2D));
  }

  if (include_chi) {
    // Gradients in temperature which drive energy flows
    add(species["energy_source"], Div_a_Grad_perp_upwind(anomalous_chi * N2D, T2D));
  }

  if (include_nu) {
    // Gradients in slow speed which drive momentum flows
    add(species["momentum_source"], Div_a_Grad_perp_upwind(anomalous_nu * N2D, V2D));
  }

  if (diagnose) {

    // void outputVars(Options& state) override {
      AUTO_TRACE();
      // Normalisations
      // auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
      // auto rho_s0 = get<BoutReal>(state["rho_s0"]);
      
      // Save particle, momentum and energy channels

      set_with_attrs(state[{std::string("anomalous_D_") + name}], anomalous_D,
                      {{"time_dimension", "t"},
                      {"units", "m^2 s^-1"},
                      {"conversion", rho_s0 * rho_s0 * Omega_ci},
                      {"standard_name", "anomalous density diffusion"},
                      {"long_name", std::string("Anomalous density diffusion of ") + name},
                      {"source", "anomalous_diffusion"}});
    // }
  }
}
