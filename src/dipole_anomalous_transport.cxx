#include "../include/dipole_anomalous_transport.hxx"

#include "../include/div_ops.hxx"
#include <bout/fv_ops.hxx>
#include <bout/output_bout_types.hxx>

using bout::globals::mesh;

DipoleAnomalousDiffusion::DipoleAnomalousDiffusion(std::string name, Options& alloptions, Solver*)
    : Component({readOnly("species:{name}:density"),
                 readIfSet("species:{name}:{optional}"),
                 readWrite("species:{name}:{output}")}),
      name(name) {
  // Normalisations
  const Options& units = alloptions["units"];
  const BoutReal rho_s0 = units["meters"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();

  const BoutReal diffusion_norm = rho_s0 * rho_s0 * Omega_ci; // m^2/s
  Options& options = alloptions[name];

  // Set in the mesh or options (or both)
  dipole_anomalous_D = 0.0;
 
  //U3D = compute_U3D();
  include_D = (mesh->get(dipole_anomalous_D, std::string("dipole_anomalous_D_") + name) == 0)
              || options.isSet("dipole_anomalous_D");
  // Option overrides mesh value
  dipole_anomalous_D = options["dipole_anomalous_D"]
                    .doc("Dipole anomalous particle diffusion coefficient [m^2/s]")
                    .withDefault(dipole_anomalous_D)
                / diffusion_norm;
  dipole_gamma = 5.0/3.0;
  dipole_gamma = options["dipole_gamma"]
                    .doc("Dipole gamma")
                    .withDefault(dipole_gamma);

  dipole_anomalous_chi = 0.0;
  include_chi = (mesh->get(dipole_anomalous_chi, std::string("dipole_chi_") + name) == 0)
                || options.isSet("dipole_anomalous_chi");
  dipole_anomalous_chi = options["dipole_anomalous_chi"]
                      .doc("Dipole anomalous thermal diffusion coefficient [m^2/s]")
                      .withDefault(dipole_anomalous_chi)
                  / diffusion_norm;
  density_floor = options["density_floor"].doc("Minimum density floor").withDefault(1e-7);
  local_U = options["local_U"]
                   .doc("Compute local U?")
                   .withDefault<bool>(false);
   U2D = compute_U2D(local_U);
  mesh->communicate(U2D);
  // dipole_anomalous_nu = 0.0;
  // include_nu = (mesh->get(dipole_anomalous_nu, std::string("dipole_nu_") + name) == 0)
  //              || options.isSet("dipole_anomalous_nu");
  // dipole_anomalous_nu = options["dipole_anomalous_nu"]
  //                    .doc("Dipole anomalous momentum diffusion coefficient [m^2/s]")
  //                    .withDefault(dipole_anomalous_nu)
  //                / diffusion_norm;

  // dipole_anomalous_sheath_flux = options["dipole_anomalous_sheath_flux"]
  //                             .doc("Allow dipole anomalous diffusion into sheath?")
  //                             .withDefault<bool>(false);

  diagnose = alloptions[name]["diagnose"]
                   .doc("Output additional diagnostics?")
                   .withDefault<bool>(false);

  substitutePermissions("name", {name});
  substitutePermissions("optional", {"temperature", "velocity", "pressure"});
  std::vector<std::string> output_vars;
  if (include_D) {
    output_vars.push_back("density_source");
    output_vars.push_back("particle_flow_xlow");
    output_vars.push_back("particle_flow_ylow");
  }
  if (include_D or include_chi) {
    output_vars.push_back("energy_source");
    output_vars.push_back("energy_flow_xlow");
    output_vars.push_back("energy_flow_ylow");
  }
  if (include_D or include_nu) {
    setPermissions(readOnly(fmt::format("species:{}:AA", name)));
    output_vars.push_back("momentum_source");
    output_vars.push_back("momentum_flow_xlow");
    output_vars.push_back("momentum_flow_ylow");
  }
  substitutePermissions("output", output_vars);
}

void DipoleAnomalousDiffusion::transform_impl(GuardedOptions& state) {

  GuardedOptions species = state["species"][name];

  // Diffusion operates on 2D (axisymmetric) profiles
  // Note: Includes diffusion in Y, so set boundary fluxes
  // to zero by imposing neumann boundary conditions.
  const Field3D N = GET_VALUE(Field3D, species["density"]);
  Field2D N2D = DC(N);



  const Field3D T = species.isSet("temperature")
                        ? GET_VALUE(Field3D, species["temperature"])
                        : 0.0;
  Field2D T2D = DC(T);

  const Field3D P = species.isSet("pressure")
                        ? GET_VALUE(Field3D, species["pressure"])
                        : N * T;
  Field2D P2D = DC(P);

  Field2D Npowg = pow(softFloor(N2D, density_floor), dipole_gamma-1.0);
  Field2D S2D = DC(P)/Npowg;

  // const Field3D V =
  //   species.isSet("velocity") ? GET_NOBOUNDARY(Field3D, species["velocity"]) : 0.0;
  // Field2D V2D = DC(V);

  // if (!dipole_anomalous_sheath_flux) {
  //   // Apply Neumann Y boundary condition, so no additional flux into boundary
  //   // Note: Not setting radial (X) boundaries since those set radial fluxes
  //   for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
  //     N2D(r.ind, mesh->ystart - 1) = N2D(r.ind, mesh->ystart);
  //     T2D(r.ind, mesh->ystart - 1) = T2D(r.ind, mesh->ystart);
  //     2D(r.ind, mesh->ystart - 1) = P2D(r.ind, mesh->ystart);
  //     V2D(r.ind, mesh->ystart - 1) = V2D(r.ind, mesh->ystart);
  //   }
  //   for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
  //     N2D(r.ind, mesh->yend + 1) = N2D(r.ind, mesh->yend);
  //     T2D(r.ind, mesh->yend + 1) = T2D(r.ind, mesh->yend);
  //     V2D(r.ind, mesh->yend + 1) = V2D(r.ind, mesh->yend);
  //   }
  // }
  transport_on = isnegative_grad_perp(P2D);
  mesh->communicate(transport_on);
  Field3D flow_xlow, flow_ylow; // Flows through cell faces

  if (include_chi && include_D) {
    // Particle diffusion. Gradients of density drive flows of particles,
    // momentum and energy. The implementation here is equivalent to an
    // advection velocity
    //
    //  v_D = - D_dp Grad_perp(N * dV) / N
    //transport_on = isnegative_grad_perp(P2D);
    Field3D Gamma_N = 1/U2D *transport_on *Div_a_Grad_perp_flows(dipole_anomalous_D, N2D*U2D,
                                                               flow_xlow, flow_ylow);
    // Field3D Gamma_N =Div_a_Grad_perp_upwind_flows(dipole_anomalous_D, N2D*U2D,flow_xlow, flow_ylow);
    add(species["density_source"], Gamma_N);
    add(species["particle_flow_xlow"], flow_xlow * transport_on);
    add(species["particle_flow_ylow"], flow_ylow * transport_on);

    // Note: Upwind operators used, or unphysical increases
    // in temperature and flow can be produced
    // auto AA = get<BoutReal>(species["AA"]);
    // add(species["momentum_source"], 1/U3D *Div_a_Grad_perp_upwind_flows(AA * V2D * dipole_anomalous_D, N2D*U2D,
    //                                                              flow_xlow, flow_ylow));
    // add(species["momentum_flow_xlow"], flow_xlow);
    // add(species["momentum_flow_ylow"], flow_ylow);

    // add(species["energy_source"],
    //     Div_a_Grad_perp_upwind_flows((3. / 2) * T2D * dipole_anomalous_D, N2D,
    //                                  flow_xlow, flow_ylow));
    // add(species["energy_flow_xlow"], flow_xlow);
    // add(species["energy_flow_ylow"], flow_ylow);
  }

  if (include_chi && include_D) {
    // Gradients in temperature that drive energy flows
    // Field3D Gamma_E = 1/U3D * pow(N2D,dipole_gamma-1.0) * isnegative_grad_perp(P2D) * Div_a_Grad_perp_upwind_flows(dipole_anomalous_chi, S2D * U2D, flow_xlow, flow_ylow);
    // Gamma_E += 1/U3D * T * (dipole_gamma-1) *isnegative_grad_perp(P2D) * Div_a_Grad_perp_upwind_flows(dipole_anomalous_D, N2D * U2D, flow_xlow, flow_ylow);
    Field3D q_inward = (3. / 2.) * 1.0/U2D * Npowg *transport_on *Div_a_Grad_perp_flows(dipole_anomalous_chi, S2D * U2D, flow_xlow, flow_ylow);
    add(species["energy_source"], q_inward);
    add(species["energy_flow_xlow"], flow_xlow* transport_on);
    add(species["energy_flow_ylow"], flow_ylow* transport_on);
    q_inward = (3. / 2.) * 1.0/U2D *transport_on * T * (dipole_gamma-1.0) * Div_a_Grad_perp_flows(dipole_anomalous_D, N2D * U2D, flow_xlow, flow_ylow);
    add(species["energy_source"], q_inward);
    add(species["energy_flow_xlow"], flow_xlow* transport_on);
    add(species["energy_flow_ylow"], flow_ylow* transport_on);
  }

  // if (include_nu) {
  //   // Gradients in flow speed that drive momentum flows
  //   auto AA = get<BoutReal>(species["AA"]);
  //   add(species["momentum_source"], Div_a_Grad_perp_upwind_flows(dipole_anomalous_nu * AA * N2D, V2D, flow_xlow, flow_ylow));
  //   add(species["momentum_flow_xlow"], flow_xlow);
  //   add(species["momentum_flow_ylow"], flow_ylow);
  // }
}

void DipoleAnomalousDiffusion::outputVars(Options& state) {
  // Normalisations
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto rho_s0 = get<BoutReal>(state["rho_s0"]);

  if (diagnose) {

      // Save particle, momentum and energy channels

      set_with_attrs(state[{std::string("dipole_anomalous_D_") + name}], dipole_anomalous_D,
                      {{"time_dimension", "t"},
                      {"units", "m^2 s^-1"},
                      {"conversion", rho_s0 * rho_s0 * Omega_ci},
                      {"standard_name", "dipole anomalous density diffusion"},
                      {"long_name", std::string("Dipole anomalous density diffusion of ") + name},
                      {"source", "dipole_anomalous_transport"}});

      // set_with_attrs(state[{std::string("dipole_anomalous_Chi_") + name}], U3D,
      //                 {{"time_dimension", "t"},
      //                 {"units", "m^2 s^-1"},
      //                 {"conversion", rho_s0},
      //                 {"standard_name", "Magnetic flux surface volume"},
      //                 {"long_name", "Magnetic flux surface volume"",
      //                 {"source", "dipole_anomalous_transport"}});

      set_with_attrs(state[{std::string("dipole_anomalous_Chi_") + name}], dipole_anomalous_chi,
                      {{"time_dimension", "t"},
                      {"units", "m^2 s^-1"},
                      {"conversion", rho_s0 * rho_s0 * Omega_ci},
                      {"standard_name", "dipole anomalous thermal diffusion"},
                      {"long_name", std::string("Dipole anomalous thermal diffusion of ") + name},
                      {"source", "dipole_anomalous_transport"}});

      set_with_attrs(state["U"], U2D,
                      {{"source", "dipole_anomalous_transport"}});
      set_with_attrs(state["transport_on"], transport_on,
                      {{"source", "dipole_anomalous_transport"}});
  }
}

const Field2D compute_U2D(bool local_U) {
  Field2D U;
  U.allocate();

  Coordinates* coord = mesh->getCoordinates();
  BoutReal s;
  for (int i = mesh->xstart; i <= mesh->xend; i++) {
    s = 0.0;
    if (!local_U) {
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      s += coord->dy(i, j) / coord->Bxy(i, j);
    }
    
      for (int j = mesh->ystart; j <= mesh->yend; j++) {
        U(i, j) = s;
      }
    }
    else
    {
      for (int j = mesh->ystart; j <= mesh->yend; j++) {
      U(i, j) = 1.0 / coord->Bxy(i, j);
    }
    }
  }
  for (int j = mesh->ystart; j <= mesh->yend; j++) {
      U(mesh->xstart-1, j) = U(mesh->xstart, j);
          U(mesh->xstart-2, j) = U(mesh->xstart, j);
  }
  for (int j = mesh->ystart; j <= mesh->yend; j++) {
      U(mesh->xend+1, j) = U(mesh->xend, j);
      U(mesh->xend+2, j) = U(mesh->xend, j);
  }




  // for (int i = mesh->xstart; i <= mesh->xend; i++) {
  //   for (int j = mesh->ystart; j <= mesh->yend; j++) {
  //     U(i, j) = coord->dy(i, j) / coord->Bxy(i, j);
  //   }
  // }
  
  return U;

}

// const Field3D compute_U3D() {
//   Field3D U;
//   U.allocate();

//   Coordinates* coord = mesh->getCoordinates();

//   for (int i = mesh->xstart; i <= mesh->xend; i++) {
//     for (int j = mesh->ystart; j <= mesh->yend; j++) {
//       for (int k = 0; k < mesh->LocalNz; k++) {
//         U(i, j, k) = coord->dy(i, j,k) / coord->Bxy(i, j,k);
//       }
//     }
//   }
//   return U;
  
// }

const Field2D isnegative_grad_perp(const Field2D& P) {


  Field2D result{zeroFrom(P)};

  Coordinates* coord = mesh->getCoordinates();

  
  // Flux in x

  int xs = mesh->xstart - 1;
  int xe = mesh->xend;

  for (int i = xs; i <= xe; i++){
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
        if ((P(i + 1, j) - P(i, j)) < 0)
        {result(i,j) = 1.0;}
        else
        {result(i,j) = 0.0;}
    }
  }
  return result;
  }


// const Field3D isnegative_grad_perp(const Field3D& P) {
//   Field3D result{zeroFrom(P)};

//   Coordinates* coord = mesh->getCoordinates();

 
//   // Flux in x

//   int xs = mesh->xstart - 1;
//   int xe = mesh->xend;

//   for (int i = xs; i <= xe; i++){
//     for (int j = mesh->ystart; j <= mesh->yend; j++) {
//       for (int k = 0; k < mesh->LocalNz; k++) {
//         // Calculate flux from i to i+1
//         if ((P(i + 1, j, k) - P(i, j, k)) * (coord->Bxy(i + 1, j, k) - coord->Bxy(i, j, k)) < 0)
//         {result(i,j,k) = 1.0;}
//         else
//         {result(i,j,k) = 0.0;}
//     }
//     }
//   }
//   return result;
//   }
  