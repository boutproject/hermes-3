#include "../include/anomalous_diffusion_3d.hxx"

#include "../include/div_ops.hxx"
#include <bout/fv_ops.hxx>
#include <bout/output_bout_types.hxx>

using bout::globals::mesh;


AnomalousDiffusion3D::AnomalousDiffusion3D(std::string name, Options& alloptions, Solver*)
  : name(name), dagp(FCI::getDagp_fv(alloptions, mesh)){
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

  if (include_D) {
    output_info.write("\tUsing mean(anomalous_D) = {}\n", mean(anomalous_D));
  }
  if (include_chi) {
    output_info.write("\tUsing mean(anomalous_chi) = {}\n", mean(anomalous_chi));
  }
  if (include_nu) {
    output_info.write("\tUsing mean(anomalous_nu) = {}\n", mean(anomalous_nu));
  }
  
  diagnose = options["diagnose"]
                   .doc("Output additional diagnostics?")
                   .withDefault<bool>(false);

  
  
}

void AnomalousDiffusion3D::transform(Options& state) {
  AUTO_TRACE();

  Options& species = state["species"][name];

  // Diffusion operates on 2D (axisymmetric) profiles
  // Note: Includes diffusion in Y, so set boundary fluxes
  // to zero by imposing neumann boundary conditions.
  Field3D N = GET_NOBOUNDARY(Field3D, species["density"]);

  Field3D T = species.isSet("temperature")
                        ? GET_NOBOUNDARY(Field3D, species["temperature"])
                        : 0.0;

  Field3D V =
      species.isSet("velocity") ? GET_NOBOUNDARY(Field3D, species["velocity"]) : 0.0;
  // if (!anomalous_sheath_flux) {
  //   // Apply Neumann Y boundary condition, so no additional flux into boundary
  //   // Note: Not setting radial (X) boundaries since those set radial fluxes
  //   // yboundary.iter([&](auto& region) {
  //   //   for (auto& pnt : region) {
  //   // 	const auto& i = pnt.ind();
  //   yboundary.iter([&](auto& region) {                                                                                                                                                                                                    
  //     for (auto& pnt : region) {
  // 	pnt.neumann_o2(N);
  // 	pnt.neumann_o2(T);
  // 	pnt.neumann_o2(V);
  //     }
  //   });

  Field3D flow_xlow, flow_zlow; // Flows through cell faces

  if (include_D) {
    // Particle diffusion. Gradients of density drive flows of particles,
    // momentum and energy. The implementation here is equivalent to an
    // advection velocity
    //
    //  v_D = - D Grad_perp(N) / N

    add(species["density_source"], (*dagp)(anomalous_D, N,
					flow_xlow, flow_zlow));
    add(species["particle_flow_xlow"], flow_xlow);
    add(species["particle_flow_zlow"], flow_zlow);

    // Note: Upwind operators used, or unphysical increases
    // in temperature and flow can be produced
    auto AA = get<BoutReal>(species["AA"]);
    add(species["momentum_source"], (*dagp)(AA * V * anomalous_D, N,
					 flow_xlow, flow_zlow));
    add(species["momentum_flow_xlow"], flow_xlow);
    add(species["momentum_flow_zlow"], flow_zlow);

    add(species["energy_source"],
        (*dagp)((3. / 2) * T * anomalous_D, N,
	     flow_xlow, flow_zlow));
    add(species["energy_flow_xlow"], flow_xlow);
    add(species["energy_flow_zlow"], flow_zlow);
  }

  if (include_chi) {
    // Gradients in temperature that drive energy flows
    add(species["energy_source"], (*dagp)(anomalous_chi * N, T, flow_xlow, flow_zlow));
    add(species["energy_flow_xlow"], flow_xlow);
    add(species["energy_flow_zlow"], flow_zlow);
  }

  if (include_nu) {
    // Gradients in flow speed that drive momentum flows
    auto AA = get<BoutReal>(species["AA"]);
    add(species["momentum_source"], (*dagp)(anomalous_nu * AA * N, V, flow_xlow, flow_zlow));
    add(species["momentum_flow_xlow"], flow_xlow);
    add(species["momentum_flow_zlow"], flow_zlow);
  }
}

void AnomalousDiffusion3D::outputVars(Options& state) {
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

