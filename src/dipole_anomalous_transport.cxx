#include "../include/dipole_anomalous_transport.hxx"

#include "../include/div_ops.hxx"
#include <bout/constants.hxx>
#include <bout/fv_ops.hxx>
#include <bout/output_bout_types.hxx>
using bout::globals::mesh;

DipoleAnomalousDiffusion::DipoleAnomalousDiffusion(std::string name, Options& alloptions,
                                                  
                                                   Solver*)
    : Component({readOnly("species:{name}:density"),
                 readIfSet("species:{name}:{optional}"),
                 readWrite("species:{name}:{output}")}),
      name(name) {
  // Normalisations
  const Options& units = alloptions["units"];
  const BoutReal rho_s0 = units["meters"];
  const BoutReal Bnorm = units["Tesla"];
  const BoutReal Tnorm = units["eV"];
  const BoutReal Nnorm = units["inv_meters_cubed"];

  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();

  const BoutReal diffusion_norm = rho_s0 * rho_s0 * Omega_ci; // m^2/s
  const BoutReal velocity_norm = rho_s0 * Omega_ci; // m/s
  Options& options = alloptions[name];

  // Set in the mesh or options (or both)
  dipole_anomalous_D = 0.0;
  dipole_quasilinear_D = 0.0;
  mu0_normalized = SI::mu0 * SI::qe * Tnorm * Nnorm / SQ(Bnorm);

  include_D = true; //(mesh->get(dipole_anomalous_D, std::string("dipole_anomalous_D_") +
                    //name) == 0)

  transport_off_factor = options["transport_off_factor"].doc("transport_off_factor").withDefault(0.0);

  dipole_anomalous_D = options["dipole_anomalous_D"]
                           .doc("Dipole anomalous particle diffusion coefficient [m^2/s]")
                           .withDefault(dipole_anomalous_D)
                       / diffusion_norm / (rho_s0 * rho_s0 * Bnorm * Bnorm);

  zero_gradient = options["zero_gradient"].doc("Apply zero gradient to anomalous diffusion coefficient at outer boundary?").withDefault<bool>(true);              
  BoutReal psi_threshold = options["psi_threshold"].doc("Threshold value for psi").withDefault(0.01);
  // factor_B = options["factor_B"].doc("Factor multiplying B in anomalous diffusion coefficient").withDefault(2.0);
  BoutReal dipole_model = options["dipole_model"].doc("Which dipole model to use? 0: no dipole transport, 1: simple diffusion, 2: quasilinear").withDefault(2);
  BoutReal qs_norm;
  Field3D B = mesh->getCoordinates()->Bxy;
  B2D = DC(B);
  Kn = 0;
  Field2D rescale_qs_D = zeroFrom(dipole_quasilinear_D) +1.0;
  Coordinates* coord = dipole_quasilinear_D.getCoordinates();
  if (dipole_model == 1) {
    qs_norm = diffusion_norm
                  * (rho_s0 * rho_s0 * Bnorm * Bnorm);
    BOUT_FOR(i, dipole_quasilinear_D.getRegion("RGN_ALL")) {
      rescale_qs_D[i] = coord->g11[i];
    }
    }
    else if (dipole_model == 2) {
      qs_norm = diffusion_norm * (rho_s0 * rho_s0);
      BOUT_FOR(i, dipole_quasilinear_D.getRegion("RGN_ALL")) {
        rescale_qs_D[i] = coord->g11[i] / B2D[i] / B2D[i];
      }
    } else if (dipole_model == 3) {
      qs_norm = diffusion_norm * (rho_s0 * rho_s0 * Bnorm * Bnorm * Bnorm * Bnorm);
      BOUT_FOR(i, dipole_quasilinear_D.getRegion("RGN_ALL")) {
        rescale_qs_D[i] = coord->g11[i] * B2D[i] * B2D[i];
      }
    } else {
      qs_norm = diffusion_norm;
      BOUT_FOR(i, dipole_quasilinear_D.getRegion("RGN_ALL")) { rescale_qs_D[i] = 1.0; }
    }

  dipole_quasilinear_D =
        options["dipole_quasilinear_D"]
            .doc("Dipole quasilinear particle diffusion coefficient [Wb^2/s]")
            .withDefault(dipole_quasilinear_D)
        / qs_norm;

  BOUT_FOR(i, dipole_quasilinear_D.getRegion("RGN_ALL")) {
    dipole_quasilinear_D[i] = dipole_quasilinear_D[i] / rescale_qs_D[i];;
  }
    
    mesh->communicate(coord->g11);
    // BOUT_FOR(i, R2D.getRegion("RGN_ALL")) { R2D[i] = R2D[i] / (coord->g11[i]); }
    Field2D Psixy;
    Psixy.allocate();
    mesh->get(Psixy, "psixy"); // get Psi
    BoutReal psi_axis, psi_bndry;
    mesh->get(psi_axis, "psi_axis");   // axis flux
    mesh->get(psi_bndry, "psi_bndry"); // edge flux
    output.write("psi_axis: {:e}, psi_bndry: {:e}\n", psi_axis, psi_bndry);

    BoutReal psi_core =
        options["psi_core"].doc("Core flux value for normalisation").withDefault(-0.48);
    BoutReal psi_edge =
        options["psi_edge"].doc("Edge flux value for normalisation").withDefault(-0.11);
    BoutReal pow_psi = options["pow_psi"]
                           .doc("Power of psi to use in diffusion coefficient")
                           .withDefault(1.0);
    output.write("psi_core: {:e}, psi_edge: {:e}", psi_core, psi_edge);
    Psixy = (Psixy - psi_core) / (psi_edge - psi_core);
    mesh->communicate(Psixy);
    output.write("psixy: {:e}, {:e}", Psixy(2, 2), Psixy(10, 30));
    BOUT_FOR(i, dipole_quasilinear_D.getRegion("RGN_ALL")) {
      dipole_quasilinear_D[i] = dipole_quasilinear_D[i] 
                                * pow(abs(Psixy[i]) + psi_threshold, pow_psi)
                                ;
    }

    dipole_eq_D = options["dipole_eq_D"]
                      .doc("Dipole equilibrium particle diffusion coefficient [Wb^2/s]")
                      .withDefault(0.0)
                  / diffusion_norm / (rho_s0 * rho_s0 * Bnorm * Bnorm);
    dipole_quasilinear_v = 0.0;
    dipole_quasilinear_v = options["dipole_quasilinear_v"]
                               .doc("Dipole quasilinear particle inward velocity [Wb/s]")
                               .withDefault(dipole_quasilinear_v)
                           / velocity_norm;
    mesh->communicate(dipole_quasilinear_v);
    // mesh->get(R2D, "R", 0.0, true);
   
    // if (mesh->firstX()) {
    //   for (int j = mesh->ystart; j <= mesh->yend; j++) {
        // dipole_quasilinear_D(mesh->xstart - 2, j) = 0.0;
        // dipole_quasilinear_D(mesh->xstart - 1, j) = 0.0;
        // dipole_quasilinear_D(mesh->xstart, j) = 0.0;
    //   }
    // }
    
    Coordinates* coord5 = dipole_eq_D.getCoordinates();
    mesh->communicate(coord5->g11);
    BOUT_FOR(i, dipole_eq_D.getRegion("RGN_ALL")) {
      dipole_eq_D[i] = dipole_eq_D[i] / (coord5->g11[i]);
    }
  Coordinates* coord2 = dipole_anomalous_D.getCoordinates();
    mesh->communicate(coord2->g11);
    BOUT_FOR(i, dipole_anomalous_D.getRegion("RGN_ALL")) {
      dipole_anomalous_D[i] = dipole_anomalous_D[i] / (coord->g11[i]) ;
    }
    // Coordinates* coord3 = dipole_quasilinear_v.getCoordinates();
    // mesh->communicate(coord3->g11);
    // // BOUT_FOR(i, R2D.getRegion("RGN_ALL")) { R2D[i] = R2D[i] / (coord->g11[i]); }
    // // BOUT_FOR(i, dipole_quasilinear_v.getRegion("RGN_ALL")) {
    // //   dipole_quasilinear_v[i] = dipole_quasilinear_v[i] / sqrt(coord3->g11[i]);
    // // }
    // output.write("dipole_quasilinear_v: {:e}", dipole_quasilinear_v(2,2));
    //     // BOUT_FOR(i, dipole_quasilinear_D.getRegion("RGN_GUARDS")) {
    //     //   dipole_quasilinear_D[i] = 0.0;
    //     // }
    // // loghtheta = log(sqrt(mesh->getCoordinates()->g22));
    //     dipole_gamma = ;
    dipole_gamma = options["dipole_gamma"].doc("Dipole gamma").withDefault(5.0 / 3.0);
    dipole_upwind = options["dipole_upwind"].doc("Dipole gamma").withDefault(false);
    dipole_quasilinear_chi = 0.0;
    dipole_anomalous_chi = 0.0;
    include_chi = true; //(mesh->get(dipole_anomalous_chi, std::string("dipole_chi_")
                        //+ name) == 0)
    //|| options.isSet("dipole_anomalous_chi");
    dipole_anomalous_chi =
        options["dipole_anomalous_chi"]
            .doc("Dipole anomalous thermal diffusion coefficient [m^2/s]")
            .withDefault(dipole_anomalous_chi)
        / diffusion_norm / (rho_s0 * rho_s0 * Bnorm * Bnorm);

    dipole_quasilinear_chi =
        options["dipole_quasilinear_chi"]
            .doc("Dipole quasilinear thermal diffusion coefficient [m^2/s]")
            .withDefault(dipole_quasilinear_chi)
        / diffusion_norm / (rho_s0 * rho_s0 * Bnorm * Bnorm);
    // model_U = options["model_U"].doc("Use U instead of B in diffusion coefficient?").withDefault(false);
    // Coordinates* coord2 = dipole_quasilinear_chi.getCoordinates();
    // mesh->communicate(coord2->g11);
    // BOUT_FOR(i, dipole_quasilinear_chi.getRegion("RGN_ALL")) {
    //   dipole_quasilinear_chi[i] =
    //       dipole_quasilinear_chi[i] / (coord2->g11[i]);
    // }
    // BOUT_FOR(i, dipole_anomalous_chi.getRegion("RGN_ALL")) {
    //   dipole_anomalous_chi[i] = dipole_anomalous_chi[i] / (coord2->g11[i]);
    // }

  
    // output.write("mu0_normalized: {:f}", mu0_normalized);
    // // if (mesh->firstX()) {
    // //   for (int j = mesh->ystart; j <= mesh->yend; j++) {
    // //     B2D(mesh->xstart - 2, j) = B2D(mesh->xstart + 3, j);
    // //     B2D(mesh->xstart - 1, j) = B2D(mesh->xstart + 3, j);
    // //     B2D(mesh->xstart, j) = B2D(mesh->xstart + 3, j);
    // //   }
    // // }
    // // if (mesh->lastX()) {
    // //   for (int j = mesh->ystart; j <= mesh->yend; j++) {
    // //     B2D(mesh->xend, j) = B2D(mesh->xend - 3, j);
    // //     B2D(mesh->xend - 1, j) = B2D(mesh->xend - 3, j);
    // //     B2D(mesh->xend - 2, j) = B2D(mesh->xend - 3, j);
    // //   }
    // // }
    // mesh->communicate(B2D);
    // mesh->communicate(loghtheta);
    // dipole_quasilinear_D =
    //     dipole_quasilinear_D/U2D;


    // if (zero_gradient && mesh->firsX()) {
    //   for (int j = mesh->ystart; j <= mesh->yend; j++) {
    //     dipole_quasilinear_D(mesh->xend + 2, j) = 0.0;
    //     dipole_quasilinear_D(mesh->xend + 1, j) = 0.0;
    //     dipole_quasilinear_D(mesh->xend, j) = 0.0;
    //   }
    // }


    mesh->communicate(dipole_anomalous_D);
    mesh->communicate(dipole_quasilinear_chi);
    mesh->communicate(dipole_quasilinear_D);
 
    // dipole_div_form = options["dipole_div_form"]
    //                       .doc("Include dipole anomalous diffusion?")
    //                       .withDefault(false);
    // dipole_naive =
    //     options["dipole_naive"].doc("Include dipole naive diffusion?").withDefault(true);

    // dipole_anomalous = options["dipole_anomalous"]
    //                        .doc("Include dipole anomalous diffusion?")
    //                        .withDefault(false);
    // density_floor =
    //     options["density_floor"].doc("Minimum density floor").withDefault(1e-7);

    // dipole_anomalous_nu = 0.0;
    // include_nu = (mesh->get(dipole_anomalous_nu, std::string("dipole_nu_") + name)
    // == 0)
    //              || options.isSet("dipole_anomalous_nu");
    // dipole_anomalous_nu = options["dipole_anomalous_nu"]
    //                    .doc("Dipole anomalous momentum diffusion coefficient
    //                    [m^2/s]") .withDefault(dipole_anomalous_nu)
    //                / diffusion_norm;

    // dipole_anomalous_sheath_flux = options["dipole_anomalous_sheath_flux"]
    //                             .doc("Allow dipole anomalous diffusion into
    //                             sheath?") .withDefault<bool>(false);
   
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

      void DipoleAnomalousDiffusion::transform_impl(GuardedOptions & state) {

        GuardedOptions species = state["species"][name];

        // Diffusion operates on 2D (axisymmetric) profiles
        // Note: Includes diffusion in Y, so set boundary fluxes
        // to zero by imposing neumann boundary conditions.
        const Field3D N = GET_VALUE(
            Field3D, species["density"]);
        Field2D N2D = DC(N);
        //N2D.applyBoundary(); this apply BC but prevent convergence to steady state
        mesh->communicate(N2D);
        const Field3D T = species.isSet("temperature")
                              ? GET_VALUE(Field3D, species["temperature"])
                              : 0.0;
        Field2D T2D = DC(T);
        mesh->communicate(T2D);

        // mesh->communicate(N2D);
        // mesh->communicate(B2D);
        // Field2D v_B = factor_B *dipole_quasilinear_D * N2D;
        // Field2D v_htheta = dipole_quasilinear_D * N2D ;
        // mesh->communicate(v_B);
        // mesh->communicate(v_htheta);
        // Field2D logB2D = log(B2D);
        // dipole_quasilinear_D_p = dipole_eq_D * B2D * B2D / mu0_normalized / T2D;
        // mesh -> communicate(dipole_quasilinear_D_p);
        mesh->communicate(dipole_quasilinear_v);
        const Field3D P = species.isSet("pressure")
                              ? GET_VALUE(Field3D, species["pressure"])
                              : N * T;
        Field2D P2D = DC(P);
        mesh->communicate(P2D);

        // const Field3D V =
        //   species.isSet("velocity") ? GET_NOBOUNDARY(Field3D, species["velocity"]) :
        //   0.0;
        // Field2D V2D = DC(V);

    
        transport_on = 1.0;
        mesh->communicate(transport_on);
        
            // -- quasilinear anomalous transport --
            Field3D S_N = Div_a_Grad_perp_upwind_flows(dipole_anomalous_D, N2D, flow_xlow,
                                                       flow_ylow);

        add(species["density_source"], S_N);
        add(species["particle_flow_xlow"], flow_xlow);
        add(species["particle_flow_ylow"], flow_ylow);

        // -- quasilinear dipole transport --
        Field2D htheta2D = 1/sqrt(mesh->getCoordinates()->g22);
        mesh->communicate(htheta2D);
        Kn = compute_Kn(N2D, B2D, htheta2D);
        mesh->communicate(Kn);
        transport_on = set_quasilinear_transport_on(Kn, transport_off_factor);
        mesh->communicate(transport_on);
        Field2D nthetaB = N2D * htheta2D / B2D;
        mesh->communicate(nthetaB);
        if (zero_gradient && mesh->firstX()) {
          for (int j = mesh->ystart; j <= mesh->yend; j++) {
            nthetaB(mesh->xstart - 2, j) = nthetaB(mesh->xstart, j);
            nthetaB(mesh->xstart - 1, j) = nthetaB(mesh->xstart, j);
          }
        }
        Field2D _dipole_quasilinear_D =
            DC(dipole_quasilinear_D) * B2D / htheta2D * transport_on;
        mesh->communicate(_dipole_quasilinear_D);
        
        S_N = 
              Div_a_Grad_perp_upwind_flows(_dipole_quasilinear_D, nthetaB, flow_xlow,
                                             flow_ylow);

        add(species["density_source"], S_N);
        add(species["particle_flow_xlow"], flow_xlow );
        add(species["particle_flow_ylow"], flow_ylow );



        // Field3D Gamma_N, q_inward, q_anomalous, Gamma_anomalous;
        // if (include_chi && include_D) {
        // Particle diffusion. Gradients of density drive flows of particles,
        // momentum and energy. The implementation here is equivalent to an
        // advection velocity
        //
        //  v_D = - D_dp Grad_perp(N * dV) / N
        // transport_on = isnegative_grad_perp(P2D);
        // if (dipole_naive) {

        // // -- quasilinear dipole transport --
        // if (mesh->firstX()) {
        //   for (int j = mesh->ystart; j <= mesh->yend; j++) {
        //     // transport_on(mesh->xstart - 2, j) = 0.0;
        //     // transport_on(mesh->xstart - 1, j) = 0.0;
        //     // transport_on(mesh->xstart, j) = 0.0;
        //     dipole_eq_D(mesh->xstart - 2, j) = 0.0;
        //     dipole_eq_D(mesh->xstart - 1, j) = 0.0;
        //     dipole_eq_D(mesh->xstart, j) = 0.0;
        //   }
        // }
        // if (mesh->lastX()) {
        //   for (int j = mesh->ystart; j <= mesh->yend; j++) {
        //     // transport_on(mesh->xend + 2, j) = 0.0;
        //     // transport_on(mesh->xend + 1, j) = 0.0;
        //     // transport_on(mesh->xend, j) = 0.0;
        //     dipole_eq_D(mesh->xend + 2, j) = 0.0;
        //     dipole_eq_D(mesh->xend+ 1, j) = 0.0;
        //     dipole_eq_D(mesh->xend, j) = 0.0;
        //   }
        // }
        // if (zero_gradient && mesh->firstX()) {
        //   for (int j = mesh->ystart; j <= mesh->yend; j++) {
        //     v_B(mesh->xstart - 2, j) = 0.0;
        //     v_B(mesh->xstart - 1, j) = 0.0;
        //     v_B(mesh->xstart, j) = 0.0;
        //     v_htheta(mesh->xstart - 2, j) = 0.0;
        //     v_htheta(mesh->xstart - 1, j) = 0.0;
        //     v_htheta(mesh->xstart, j) = 0.0;
        //     // dipole_eq_D(mesh->xstart - 2, j) = 0.0;
        //     // dipole_eq_D(mesh->xstart - 1, j) = 0.0;
        //     // dipole_eq_D(mesh->xstart, j) = 0.0;
        //   }
        // }
        // mesh->communicate(transport_on);
        // mesh->communicate(dipole_eq_D);
        // if (model_U) {
        //   Gamma_N = -transport_on
        //             * Div_a_Grad_perp_upwind_flows(v_B, logU2D, flow_xlow, flow_ylow);
        // }
        // else
        // {
        // Gamma_N = -transport_on
        //           * Div_a_Grad_perp_upwind_flows(v_B, logB2D, flow_xlow, flow_ylow);
        // }
        // add(species["density_source"], Gamma_N);
        // add(species["particle_flow_xlow"], flow_xlow * transport_on);
        // add(species["particle_flow_ylow"], flow_ylow * transport_on);

        
        //   Gamma_N = -transport_on
        //             * Div_a_Grad_perp_upwind_flows(v_htheta, loghtheta, flow_xlow, flow_ylow);
        
        // add(species["density_source"], Gamma_N);
        // add(species["particle_flow_xlow"], flow_xlow * transport_on);
        // add(species["particle_flow_ylow"], flow_ylow * transport_on);

        // Gamma_N = Div_a_Grad_perp_upwind_flows(dipole_quasilinear_D_p, logB2D,
        //                                          flow_xlow, flow_ylow);
        // add(species["density_source"], Gamma_N);
        // add(species["particle_flow_xlow"], flow_xlow );
        // add(species["particle_flow_ylow"], flow_ylow);

        // Gamma_N = Div_a_Grad_perp_upwind_flows(dipole_eq_D, N2D,
        //                                          flow_xlow, flow_ylow);
        // add(species["density_source"], Gamma_N);
        // add(species["particle_flow_xlow"], flow_xlow);
        // add(species["particle_flow_ylow"], flow_ylow);

        // add(species["density_source"], -Conv_flows(dipole_quasilinear_v, N2D, flow_xlow));
        // add(species["particle_flow_xlow"], flow_xlow);
    

        // +Div_a_Grad_perp_upwind_flows(dipole_quasilinear_D, N2D, flow_xlow, flow_ylow)
        //     - Div_a_Grad_perp_upwind_flows(dipole_quasilinear_D_p, P2D, flow_xlow,
        //                                    flow_ylow);
        //- 
        //+ Div_a_Grad_perp_upwind_flows(dipole_anomalous_D, N2D, flow_xlow,
        //                 flow_ylow);
        // }
        // else {
        //   if (dipole_upwind) {
        //     Gamma_N = 1 / U2D * transport_on
        //               * Div_a_Grad_perp_upwind_flows(dipole_quasilinear_D,
        //               N2D * U2D,
        //                                              flow_xlow, flow_ylow);
        //   } else {
        //     if (!dipole_div_form) {
        //       Gamma_N = 1 / U2D * transport_on
        //                 * Div_a_Grad_perp_flows(dipole_quasilinear_D, N2D *
        //                 U2D, flow_xlow,
        //                                         flow_ylow);
        //     } else {
        //       Gamma_N = transport_on
        //                 * Div_a_Grad_perp_flows(1 / U2D *
        //                 dipole_quasilinear_D, N2D * U2D,
        //                                         flow_xlow, flow_ylow);
        //     }
        //   }
        // }
        // Gamma_N = 1/U2D *transport_on *Div_a_Grad_per

        // Field3D Gamma_N =Div_a_Grad_perp_upwind_flows(dipole_anomalous_D,
        // N2D*U2D,flow_xlow, flow_ylow);
        // add(species["density_source"], Gamma_N);
        // add(species["particle_flow_xlow"], flow_xlow * transport_on);
        // add(species["particle_flow_ylow"], flow_ylow * transport_on);

        // Note: Upwind operators used, or unphysical increases
        // in temperature and flow can be produced
        // auto AA = get<BoutReal>(species["AA"]);
        // add(species["momentum_source"], 1/U3D *Div_a_Grad_perp_upwind_flows(AA * V2D
        // * dipole_anomalous_D, N2D*U2D,
        //                                                              flow_xlow,
        //                                                              flow_ylow));
        // add(species["momentum_flow_xlow"], flow_xlow);
        // add(species["momentum_flow_ylow"], flow_ylow);

        // add(species["energy_source"],
        //     Div_a_Grad_perp_upwind_flows((3. / 2) * T2D * dipole_anomalous_D, N2D,
        //                                  flow_xlow, flow_ylow));
        // add(species["energy_flow_xlow"], flow_xlow);
        // add(species["energy_flow_ylow"], flow_ylow);

        // Gradients in temperature that drive energy flows
        // Field3D Gamma_E = 1/U3D * pow(N2D,dipole_gamma-1.0) *
        // isnegative_grad_perp(P2D) *
        // Div_a_Grad_perp_upwind_flows(dipole_anomalous_chi, S2D * U2D, flow_xlow,
        // flow_ylow); Gamma_E += 1/U3D * T * (dipole_gamma-1)
        // *isnegative_grad_perp(P2D) * Div_a_Grad_perp_upwind_flows(dipole_anomalous_D,
        // N2D * U2D, flow_xlow, flow_ylow);
        //   Field2D Npowg = pow(softFloor(N2D, density_floor), dipole_gamma - 1.0);
        //   Field2D U2Dpowg2 = pow(U2D, dipole_gamma - 1.0);
        //   mesh->communicate(Npowg);
        //   mesh->communicate(U2Dpowg2);
        //   Field2D TU2D = T2D * U2Dpowg2;
        //   Field2D S2D = DC(P) / Npowg;
        //   mesh->communicate(S2D);
        //   mesh->communicate(TU2D);
        //   if (dipole_naive) {
        //     q_inward = transport_on
        //                * Div_a_Grad_perp_flows(dipole_quasilinear_chi * N2D, TU2D,
        //                                        flow_xlow, flow_ylow);
        //   } else {
        //     if (dipole_upwind) {
        //       q_inward = (3. / 2.) * 1.0 / U2D * Npowg * transport_on
        //                  * Div_a_Grad_perp_upwind_flows(dipole_quasilinear_chi, S2D * U2D,
        //                                                 flow_xlow, flow_ylow);

        //       // q_inward = transport_on
        //       //            * Div_a_Grad_perp_upwind_flows((3. / 2.) * 1.0 / U2D * Npowg
        //       //                                               *  dipole_quasilinear_chi,
        //       //                                           S2D * U2D, flow_xlow,
        //       //                                           flow_ylow);

        //     } else {
        //       if (!dipole_div_form) {
        //         q_inward = (3. / 2.) * 1.0 / U2D * Npowg * transport_on
        //                    * Div_a_Grad_perp_flows(dipole_quasilinear_chi, S2D * U2D,
        //                                            flow_xlow, flow_ylow);

        //       } else {
        //         q_inward = transport_on
        //                    * Div_a_Grad_perp_flows((3. / 2.) * 1.0 / U2D * Npowg
        //                                                * dipole_quasilinear_chi,
        //                                            S2D * U2D, flow_xlow, flow_ylow);
        //       }
        //     }
        //   }
        //   // ef_dipole_quasilinear_xlow += flow_xlow * transport_on;
        //   // ef_dipole_quasilinear_ylow += flow_ylow * transport_on;
        //   add(species["energy_source"], q_inward);
        //   add(species["energy_flow_xlow"], flow_xlow * transport_on);
        //   add(species["energy_flow_ylow"], flow_ylow * transport_on);
        //   if (dipole_naive) {
        //     q_inward =
        //         transport_on
        //         * Div_a_Grad_perp_upwind_flows((3. / 2) * T2D * dipole_quasilinear_D,
        //                                        N2D * U2D, flow_xlow, flow_ylow);
        //   } else {
        //     if (dipole_upwind) {
        //       q_inward = (3. / 2.) * 1.0 / U2D * T * transport_on
        //                  * Div_a_Grad_perp_upwind_flows((dipole_gamma - 1.0)
        //                                                     * dipole_quasilinear_D,
        //                                                 N2D * U2D, flow_xlow, flow_ylow);
        //       // q_inward = transport_on
        //       //            * Div_a_Grad_perp_upwind_flows((3. / 2.) * 1.0 / U2D * T *
        //       //            (dipole_gamma - 1.0)
        //       //                                        * dipole_quasilinear_D,
        //       //                                    N2D * U2D, flow_xlow, flow_ylow);
        //     } else {
        //       // q_inward = (3. / 2.) * 1.0 / U2D * transport_on * T * (dipole_gamma
        //       // - 1.0)
        //       //            * Div_a_Grad_perp_flows(dipole_quasilinear_D, N2D * U2D,
        //       //              flow_xlow, flow_ylow);
        //       if (!dipole_div_form) {
        //         q_inward = (dipole_gamma - 1.0) * (3. / 2.) * 1.0 / U2D * T2D
        //                    * transport_on
        //                    * Div_a_Grad_perp_flows(dipole_quasilinear_D, N2D * U2D,
        //                                            flow_xlow, flow_ylow);
        //       } else {
        //         q_inward = transport_on
        //                    * Div_a_Grad_perp_flows((3. / 2.) * 1.0 / U2D * T2D
        //                                                * (dipole_gamma - 1.0)
        //                                                * dipole_quasilinear_D,
        //                                            N2D * U2D, flow_xlow, flow_ylow);
        //       }
        //     }
        //   }
        //   // ef_dipole_quasilinear_xlow += flow_xlow * transport_on;
        //   // ef_dipole_quasilinear_ylow += flow_ylow * transport_on;
        //   add(species["energy_source"], q_inward);
        //   add(species["energy_flow_xlow"], flow_xlow * transport_on);
        //   add(species["energy_flow_ylow"], flow_ylow * transport_on);
        //   if (dipole_anomalous) {
        //     // if (dipole_upwind) {
        //     //   Gamma_anomalous = Div_a_Grad_perp_upwind_flows(dipole_anomalous_D, N2D,
        //     //                                                  flow_xlow, flow_ylow);
        //     // } else {
        //     //   Gamma_anomalous =
        //     //       Div_a_Grad_perp_flows(dipole_anomalous_D, N2D, flow_xlow, flow_ylow);
        //     //   ;
        //     // }
        //     // add(species["density_source"], Gamma_anomalous);
        //     // // pf_dipole_anomalous_xlow += flow_xlow ;
        //     // // pf_dipole_anomalous_ylow += flow_ylow ;
        //     // add(species["particle_flow_xlow"], flow_xlow);
        //     // add(species["particle_flow_ylow"], flow_ylow);

        //     if (dipole_upwind) {
        //       q_anomalous = Div_a_Grad_perp_upwind_flows(dipole_anomalous_chi * N2D, T2D,
        //                                                  flow_xlow, flow_ylow);
        //     } else {
        //       q_anomalous = Div_a_Grad_perp_flows(dipole_anomalous_chi * N2D, T2D,
        //                                           flow_xlow, flow_ylow);
        //     }
        //     add(species["energy_source"], q_anomalous);
        //     add(species["energy_flow_xlow"], flow_xlow);
        //     add(species["energy_flow_ylow"], flow_ylow);
        //     if (dipole_upwind) {
        //       q_anomalous = Div_a_Grad_perp_upwind_flows(
        //           (3. / 2) * T2D * dipole_anomalous_D, N2D, flow_xlow, flow_ylow);
        //     } else {
        //       q_anomalous = Div_a_Grad_perp_flows((3. / 2) * T2D * dipole_anomalous_D,
        //                                           N2D, flow_xlow, flow_ylow);
        //     }
        //     add(species["energy_source"], q_anomalous);
        //     add(species["energy_flow_xlow"], flow_xlow);
        //     add(species["energy_flow_ylow"], flow_ylow);
        //   }
        // }
      }

      void DipoleAnomalousDiffusion::outputVars(Options & state) {
        // Normalisations
        auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
        auto rho_s0 = get<BoutReal>(state["rho_s0"]);
        auto Bnorm = get<BoutReal>(state["Bnorm"]);
        auto Nnorm = get<BoutReal>(state["Nnorm"]);

        if (diagnose) {

      // Save particle, momentum and energy channels

          set_with_attrs(
              state[{std::string("dipole_quasilinear_D_") + name}], dipole_quasilinear_D,
              {{"time_dimension", "t"},
               {"units", "m^2 s^-1"},
               {"conversion",
                rho_s0 * rho_s0 * Omega_ci * rho_s0 * rho_s0 * Bnorm * Bnorm},
               {"standard_name", "dipole quasilinear density diffusion"},
               {"long_name",
                std::string("Dipole quasilinear density diffusion of ") + name},
               {"source", "dipole_anomalous_transport"}});

          // set_with_attrs(state[{std::string("dipole_anomalous_Chi_") + name}], U3D,
          //                 {{"time_dimension", "t"},
          //                 {"units", "m^2 s^-1"},
          //                 {"conversion", rho_s0},
          //                 {"standard_name", "Magnetic flux surface volume"},
          //                 {"long_name", "Magnetic flux surface volume"",
          //                 {"source", "dipole_anomalous_transport"}});

          set_with_attrs(
              state[{std::string("dipole_quasilinear_Chi_") + name}],
              dipole_quasilinear_chi,
              {{"time_dimension", "t"},
               {"units", "m^2 s^-1"},
               {"conversion", rho_s0 * rho_s0 * Omega_ci},
               {"standard_name", "dipole quasilinear thermal diffusion"},
               {"long_name",
                std::string("Dipole quasilinear thermal diffusion of ") + name},
               {"source", "dipole_anomalous_transport"}});

          // set_with_attrs(
          //     state["U"], U2D,
          //     {{"source", "dipole_anomalous_transport"}, {"units", "m^3 Wb^-1"}});
          //,
          //{"conversion", rho_s0 / Bnorm} });

          set_with_attrs(state[fmt::format("transport_on_{}", name)], transport_on,
                         {
                             {"source", "dipole_anomalous_transport"},
                             {"time_dimension", "t"},
                             {"units", "1"},
                             {"long_name", "Dipole transport on/off switch"},
                         });

          if (Kn.isAllocated()) {
            set_with_attrs(
                state[fmt::format("Kn_{}", name)], Kn,
                {{"time_dimension", "t"},
                 {"units", ""},
                 {"standard_name", "Kn"},
                 {"long_name", name + " dipole quasilinear Kn"},
                 {"species", name},
                 {"source", "dipole_anomalous_transport"}});
          }

          if (flow_xlow.isAllocated()) {
            set_with_attrs(
                state[fmt::format("pf{}_dipole_xlow", name)], flow_xlow,
                {{"time_dimension", "t"},
                 {"units", "s^-1"},
                 {"conversion", rho_s0 * SQ(rho_s0) * Nnorm * Omega_ci},
                 {"standard_name", "dipole particle flow"},
                 {"long_name", name + " dipole transport particle flow in X. "},
                 {"species", name},
                 {"source", "dipole_anomalous_transport"}});
          }
          if (flow_ylow.isAllocated()) {
            set_with_attrs(
                state[fmt::format("pf{}_dipole_ylow", name)], flow_ylow,
                {{"time_dimension", "t"},
                 {"units", "s^-1"},
                 {"conversion", rho_s0 * SQ(rho_s0) * Nnorm * Omega_ci},
                 {"standard_name", "particle flow"},
                 {"long_name", name + " dipole transport particle flow in Y. "},
                 {"species", name},
                 {"source", "dipole_anomalous_transport"}});
          }
        }
      }

      const void compute_U2D(Field2D & U, bool zero_inner_gradient_U,
                             bool zero_outer_gradient_U) {
        // Load and save coordinate variables
        Coordinates* coord = mesh->getCoordinates();
        if (!(mesh->sourceHasVar("dVdpsi"))) {
          throw BoutException(
              "Grid input does not contain dVdpsi needed for dipole transport.\n");
        }
        mesh->get(U, "dVdpsi", 0.0, true);
        mesh->communicate(U);
        for (int i = mesh->xstart - 2; i <= mesh->xend + 2; i++) {
          for (int j = mesh->ystart; j <= mesh->yend; j++) {
            U(i, j) = U(i, j) / 230.0; // Normalise to value at inner boundary
          }
        }
        if (mesh->firstX() and zero_inner_gradient_U) {
          for (int j = mesh->ystart; j <= mesh->yend; j++) {

            U(mesh->xstart - 2, j) = U(mesh->xstart + 1, j);
            U(mesh->xstart - 1, j) = U(mesh->xstart + 1, j);
            U(mesh->xstart, j) = U(mesh->xstart + 1, j);
          }
        }
        if (mesh->lastX() and zero_outer_gradient_U) {
          for (int j = mesh->ystart; j <= mesh->yend; j++) {

            U(mesh->xend + 2, j) = U(mesh->xend - 1, j);
            U(mesh->xend + 1, j) = U(mesh->xend - 1, j);
            U(mesh->xend, j) = U(mesh->xend - 1, j);
          }
        }
      }
      const Field2D isnegative_dnBthetadpsi(const Field2D& N2D, const Field2D& B2D, const Field2D& htheta2D,
                                            BoutReal transport_off_factor) {

        Field2D result{zeroFrom(N2D)};


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
  