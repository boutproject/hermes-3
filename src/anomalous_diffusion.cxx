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
      N2D(r.ind, mesh->ystart - 1) = N2D(r.ind, mesh->ystart);
      T2D(r.ind, mesh->ystart - 1) = T2D(r.ind, mesh->ystart);
      V2D(r.ind, mesh->ystart - 1) = V2D(r.ind, mesh->ystart);
    }
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      N2D(r.ind, mesh->yend + 1) = N2D(r.ind, mesh->yend);
      T2D(r.ind, mesh->yend + 1) = T2D(r.ind, mesh->yend);
      V2D(r.ind, mesh->yend + 1) = V2D(r.ind, mesh->yend);
    }
  }

  Field3D flow_xlow, flow_ylow; // Flows through cell faces

  if (include_D) {
    // Particle diffusion. Gradients of density drive flows of particles,
    // momentum and energy. The implementation here is equivalent to an
    // advection velocity
    //
    //  v_D = - D Grad_perp(N) / N

    add(species["density_source"], Div_a_Grad_perp_upwind_flows(anomalous_D, N2D,
                                                                flow_xlow, flow_ylow));
    add(species["particle_flow_xlow"], flow_xlow);
    add(species["particle_flow_ylow"], flow_ylow);

    // Note: Upwind operators used, or unphysical increases
    // in temperature and flow can be produced
    auto AA = get<BoutReal>(species["AA"]);
    add(species["momentum_source"], Div_a_Grad_perp_upwind_flows(AA * V2D * anomalous_D, N2D,
                                                                 flow_xlow, flow_ylow));
    add(species["momentum_flow_xlow"], flow_xlow);
    add(species["momentum_flow_ylow"], flow_ylow);

    add(species["energy_source"],
        Div_a_Grad_perp_upwind_flows((3. / 2) * T2D * anomalous_D, N2D,
                                     flow_xlow, flow_ylow));
    add(species["energy_flow_xlow"], flow_xlow);
    add(species["energy_flow_ylow"], flow_ylow);
  }

  if (include_chi) {
    // Gradients in temperature that drive energy flows
    add(species["energy_source"], Div_a_Grad_perp_upwind_flows(anomalous_chi * N2D, T2D, flow_xlow, flow_ylow));
    add(species["energy_flow_xlow"], flow_xlow);
    add(species["energy_flow_ylow"], flow_ylow);
  }

  if (include_nu) {
    // Gradients in flow speed that drive momentum flows
    auto AA = get<BoutReal>(species["AA"]);
    add(species["momentum_source"], Div_a_Grad_perp_upwind_flows(anomalous_nu * AA * N2D, V2D, flow_xlow, flow_ylow));
    add(species["momentum_flow_xlow"], flow_xlow);
    add(species["momentum_flow_ylow"], flow_ylow);
  }
}

void AnomalousDiffusion::outputVars(Options& state) {
  AUTO_TRACE();
  // Normalisations
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto rho_s0 = get<BoutReal>(state["rho_s0"]);

  if (diagnose) {

      AUTO_TRACE();
      // Save particle, momentum and energy channels

      set_with_attrs(state[{std::string("anomalous_D_") + name}], anomalous_D,
                      {{"time_dimension", "t"},
                      {"units", "m^2 s^-1"},
                      {"conversion", rho_s0 * rho_s0 * Omega_ci},
                      {"standard_name", "anomalous density diffusion"},
                      {"long_name", std::string("Anomalous density diffusion of ") + name},
                      {"source", "anomalous_diffusion"}});

      set_with_attrs(state[{std::string("anomalous_Chi_") + name}], anomalous_chi,
                      {{"time_dimension", "t"},
                      {"units", "m^2 s^-1"},
                      {"conversion", rho_s0 * rho_s0 * Omega_ci},
                      {"standard_name", "anomalous thermal diffusion"},
                      {"long_name", std::string("Anomalous thermal diffusion of ") + name},
                      {"source", "anomalous_diffusion"}});

      set_with_attrs(state[{std::string("anomalous_nu_") + name}], anomalous_nu,
                      {{"time_dimension", "t"},
                      {"units", "m^2 s^-1"},
                      {"conversion", rho_s0 * rho_s0 * Omega_ci},
                      {"standard_name", "anomalous momentum diffusion"},
                      {"long_name", std::string("Anomalous momentum diffusion of ") + name},
                      {"source", "anomalous_diffusion"}});
  }
}

