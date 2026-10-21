#include "../include/dipole_anomalous_transport.hxx"

#include "../include/div_ops.hxx"
#include <bout/constants.hxx>
#include <bout/fv_ops.hxx>
#include <bout/output_bout_types.hxx>
using bout::globals::mesh;

DipoleAnomalousDiffusion::DipoleAnomalousDiffusion(std::string name, Options& alloptions,
                                                   Solver*)
    : NamedComponent(name, {readOnly("species:{name}:density"),
                            readIfSet("species:{name}:{optional}"),
                            readWrite("species:{name}:{output}")}),
      name(name) {
  // Normalisations
  const Options& units = alloptions["units"];
  const BoutReal rho_s0 = units["meters"];
  const BoutReal Bnorm = units["Tesla"];

  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();

  const BoutReal diffusion_norm = rho_s0 * rho_s0 * Omega_ci; // m^2/s
  Options& options = alloptions[name];

  Field3D B = mesh->getCoordinates()->Bxy();
  B2D = DC(B);
  Coordinates* coord = mesh->getCoordinates();
  auto g11 = coord->g11();
  mesh->communicate(g11);

  // Converts an input diffusion coefficient to the normalised D used in the
  // operators, depending on the model:
  // 1: D = value/(R^2 B^2), value in [Wb^2 rad^-2 s^-1] (OSTR-Mod.1 eq. 38)
  // 2: D = value/R^2, value in [m^4 s^-1]
  // 3: D = value/(R^2 B^4), value in [T^4 m^4 s^-1]
  // other: D = value in [m^2 s^-1]
  // R B is taken from g11 = (R Bp)^2, so B = Bp is assumed (no toroidal field)
  const std::string model_doc =
      "1: D = value/(R^2 B^2) with value in [Wb^2 rad^-2 s^-1] (OSTR-Mod.1 eq. 38), "
      "2: D = value/R^2 with value in [m^4 s^-1], 3: D = value/(R^2 B^4) with value "
      "in [T^4 m^4 s^-1], other: D = value in [m^2 s^-1]";
  auto model_scaling = [&](int model) {
    Field2D scale{zeroFrom(B2D)};
    BOUT_FOR(i, scale.getRegion("RGN_ALL")) {
      if (model == 1) {
        scale[i] = 1. / (diffusion_norm * SQ(rho_s0 * Bnorm) * g11[i]);
      } else if (model == 2) {
        scale[i] = SQ(B2D[i]) / (diffusion_norm * SQ(rho_s0) * g11[i]);
      } else if (model == 3) {
        scale[i] = 1. / (diffusion_norm * SQ(SQ(Bnorm) * rho_s0) * g11[i] * SQ(B2D[i]));
      } else {
        scale[i] = 1. / diffusion_norm;
      }
    }
    return scale;
  };

  const int dipole_anomalous_model =
      options["dipole_anomalous_model"]
          .doc("Anomalous diffusivity model for dipole_anomalous_D and "
               "dipole_anomalous_chi. "
               + model_doc)
          .withDefault(1);
  const Field2D anomalous_scaling = model_scaling(dipole_anomalous_model);

  dipole_anomalous_D = 0.0;
  dipole_anomalous_D = options["dipole_anomalous_D"]
                           .doc("Dipole anomalous particle diffusion; units set by "
                                "dipole_anomalous_model")
                           .withDefault(dipole_anomalous_D)
                       * anomalous_scaling;

  dipole_anomalous_chi = 0.0;
  dipole_anomalous_chi = options["dipole_anomalous_chi"]
                             .doc("Dipole anomalous thermal diffusion; units set by "
                                  "dipole_anomalous_model")
                             .withDefault(dipole_anomalous_chi)
                         * anomalous_scaling;

  const bool dipole_inner_noflux =
      options["dipole_inner_noflux"]
          .doc("Default for dipole_inner_noflux_particle and dipole_inner_noflux_heat")
          .withDefault<bool>(false);
  dipole_inner_noflux_particle =
      options["dipole_inner_noflux_particle"]
          .doc("No quasilinear particle (and momentum) flux through the inner boundary?")
          .withDefault<bool>(dipole_inner_noflux);
  // Energy through the wall is 3/2 [n^(gamma-1) Gamma_s + (gamma-1) T Gamma_n]: with
  // this, only the energy carried by escaping particles remains
  dipole_inner_noflux_heat =
      options["dipole_inner_noflux_heat"]
          .doc("No quasilinear entropy (heat) flux through the inner boundary?")
          .withDefault<bool>(dipole_inner_noflux);

  const int dipole_quasilinear_model =
      options["dipole_quasilinear_model"]
          .doc("Quasilinear diffusivity model for dipole_quasilinear_D. " + model_doc)
          .withDefault(1);

  dipole_quasilinear_D = 0.0;
  dipole_quasilinear_D = options["dipole_quasilinear_D"]
                             .doc("Dipole quasilinear diffusion flux function Dpsi; units "
                                  "set by dipole_quasilinear_model")
                             .withDefault(dipole_quasilinear_D)
                         * model_scaling(dipole_quasilinear_model);

  dipole_gamma = options["dipole_gamma"].doc("Adiabatic index").withDefault(5.0 / 3.0);
  dipole_upwind =
      options["dipole_upwind"]
          .doc("Use upwind instead of central Div_a_Grad_perp for dipole fluxes?")
          .withDefault(true);

  // Dipole transport only on closed field lines (periodic in y). On the face between
  // the last closed and first open cell the operators still see half (central) or
  // all-or-nothing (upwind) of the core coefficient
  Field2D closed{zeroFrom(B2D)};
  for (int x = 0; x < mesh->LocalNx; x++) {
    for (int y = 0; y < mesh->LocalNy; y++) {
      closed(x, y) = mesh->periodicY(x) ? 1.0 : 0.0;
    }
  }
  dipole_anomalous_D *= closed;
  dipole_anomalous_chi *= closed;
  dipole_quasilinear_D *= closed;

  mesh->communicate(dipole_anomalous_D);
  mesh->communicate(dipole_anomalous_chi);
  mesh->communicate(dipole_quasilinear_D);

  diagnose = alloptions[name]["diagnose"]
                 .doc("Output additional diagnostics?")
                 .withDefault<bool>(false);

  substitutePermissions("name", {name});
  substitutePermissions("optional", {"temperature", "velocity", "pressure"});
  setPermissions(readOnly(fmt::format("species:{}:AA", name)));
  substitutePermissions("output", {"density_source", "particle_flow_xlow",
                                   "particle_flow_ylow", "energy_source",
                                   "energy_flow_xlow", "energy_flow_ylow",
                                   "momentum_source", "momentum_flow_xlow",
                                   "momentum_flow_ylow"});
}

void DipoleAnomalousDiffusion::transform_impl(GuardedOptions& state) {
  GuardedOptions species = state["species"][name];
  const auto AA = get<BoutReal>(species["AA"]);

  // Diffusion operates on 2D (axisymmetric) profiles
  const Field3D N = GET_VALUE(Field3D, species["density"]);
  Field2D N2D = DC(N);
  // N2D.applyBoundary() would set the BC but prevents convergence to steady state
  mesh->communicate(N2D);
  const Field3D T =
      species.isSet("temperature") ? GET_VALUE(Field3D, species["temperature"]) : 0.0;
  Field2D T2D = DC(T);
  mesh->communicate(T2D);
  const Field3D P =
      species.isSet("pressure") ? GET_VALUE(Field3D, species["pressure"]) : N * T;
  Field2D P2D = DC(P);
  mesh->communicate(P2D);
  const Field3D V =
      species.isSet("velocity") ? GET_NOBOUNDARY(Field3D, species["velocity"]) : 0.0;
  Field2D V2D = DC(V);
  mesh->communicate(V2D);

  auto div_flows = [&](const Field2D& a, const Field2D& f, Field3D& fx, Field3D& fy) {
    return dipole_upwind ? Div_a_Grad_perp_upwind_flows(a, f, fx, fy)
                         : Div_a_Grad_perp_flows(a, f, fx, fy);
  };
  Field3D fx, fy; // Flows through cell faces

  // -- anomalous transport: particle diffusion carrying momentum and energy,
  //    as in anomalous_diffusion --
  add(species["density_source"], div_flows(dipole_anomalous_D, N2D, fx, fy));
  add(species["particle_flow_xlow"], fx);
  add(species["particle_flow_ylow"], fy);

  add(species["momentum_source"], div_flows(AA * V2D * dipole_anomalous_D, N2D, fx, fy));
  add(species["momentum_flow_xlow"], fx);
  add(species["momentum_flow_ylow"], fy);

  add(species["energy_source"],
      div_flows((3. / 2) * T2D * dipole_anomalous_D, N2D, fx, fy));
  add(species["energy_flow_xlow"], fx);
  add(species["energy_flow_ylow"], fy);

  // Gradients in temperature that drive energy flows, as in anomalous_diffusion
  add(species["energy_source"], div_flows(dipole_anomalous_chi * N2D, T2D, fx, fy));
  add(species["energy_flow_xlow"], fx);
  add(species["energy_flow_ylow"], fy);

  // -- quasilinear dipole transport, OSTR-Mod.1 section 4.4 --
  Field2D htheta2D = 1 / sqrt(mesh->getCoordinates()->g22());
  mesh->communicate(htheta2D);
  auto inner_noflux = [&](Field2D& f, bool noflux) {
    if (noflux && mesh->firstX()) {
      for (int j = mesh->ystart; j <= mesh->yend; j++) {
        f(mesh->xstart - 2, j) = f(mesh->xstart, j);
        f(mesh->xstart - 1, j) = f(mesh->xstart, j);
      }
    }
  };
  Field2D nthetaB = N2D * htheta2D / B2D;
  mesh->communicate(nthetaB);
  inner_noflux(nthetaB, dipole_inner_noflux_particle);
  Field2D a_QL = DC(dipole_quasilinear_D) * B2D / htheta2D;
  mesh->communicate(a_QL);

  // F_n, eqs. (52), (54)
  F_n = div_flows(a_QL, nthetaB, flow_xlow, flow_ylow);
  add(species["density_source"], F_n);
  add(species["particle_flow_xlow"], flow_xlow);
  add(species["particle_flow_ylow"], flow_ylow);

  // Parallel momentum carried by the quasilinear particle flux
  add(species["momentum_source"], div_flows(AA * V2D * a_QL, nthetaB, fx, fy));
  add(species["momentum_flow_xlow"], fx);
  add(species["momentum_flow_ylow"], fy);

  // F_s, eq. (53): s = p n^(1-gamma) is conserved by the interchange displacement
  const Field2D Nlim = floor(N2D, 1e-8);
  Field2D sthetaB = P2D * pow(Nlim, 1. - dipole_gamma) * htheta2D / B2D;
  inner_noflux(sthetaB, dipole_inner_noflux_heat);
  F_s = div_flows(a_QL, sthetaB, fx, fy);
  // dp/dt = n^(gamma-1) F_s + (gamma-1) T F_n, energy = 3/2 p, eqs. (55), (57)
  add(species["energy_source"],
      1.5 * (pow(Nlim, dipole_gamma - 1.) * F_s + (dipole_gamma - 1.) * T2D * F_n));

  // Energy flows through cell faces: 3/2 q_p, eq. (58). Only the conservative
  // part; the compression work W_a is in the source but not in the flows
  const Field2D Ngm1 = pow(Nlim, dipole_gamma - 1.);
  Field3D eflow_xlow{zeroFrom(fx)}, eflow_ylow{zeroFrom(fy)};
  for (int i = mesh->xstart; i <= mesh->xend + 1; i++) {
    for (int j = mesh->ystart; j <= mesh->yend + 1; j++) {
      for (int k = 0; k < mesh->LocalNz; k++) {
        eflow_xlow(i, j, k) =
            1.5 * (0.5 * (Ngm1(i - 1, j) + Ngm1(i, j)) * fx(i, j, k)
                   + (dipole_gamma - 1.) * 0.5 * (T2D(i - 1, j) + T2D(i, j))
                         * flow_xlow(i, j, k));
        eflow_ylow(i, j, k) =
            1.5 * (0.5 * (Ngm1(i, j - 1) + Ngm1(i, j)) * fy(i, j, k)
                   + (dipole_gamma - 1.) * 0.5 * (T2D(i, j - 1) + T2D(i, j))
                         * flow_ylow(i, j, k));
      }
    }
  }
  add(species["energy_flow_xlow"], eflow_xlow);
  add(species["energy_flow_ylow"], eflow_ylow);
}

void DipoleAnomalousDiffusion::outputVars(Options& state) {
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto rho_s0 = get<BoutReal>(state["rho_s0"]);
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);

  if (diagnose) {
    set_with_attrs(
        state[{std::string("dipole_quasilinear_D_") + name}], dipole_quasilinear_D,
        {{"time_dimension", "t"},
         {"units", "m^2 s^-1"},
         {"conversion", rho_s0 * rho_s0 * Omega_ci},
         {"standard_name", "dipole quasilinear density diffusion"},
         {"long_name", std::string("Dipole quasilinear density diffusion of ") + name},
         {"source", "dipole_anomalous_transport"}});

    if (F_n.isAllocated()) {
      set_with_attrs(state[fmt::format("F_n_{}_dipole", name)], F_n,
                     {{"time_dimension", "t"},
                      {"units", "m^-3 s^-1"},
                      {"conversion", Nnorm * Omega_ci},
                      {"standard_name", "dipole quasilinear density source"},
                      {"long_name", name + " dipole quasilinear density source F_n"},
                      {"species", name},
                      {"source", "dipole_anomalous_transport"}});
    }
    if (F_s.isAllocated()) {
      // s = p n^(1-gamma)
      set_with_attrs(
          state[fmt::format("F_s_{}_dipole", name)], F_s,
          {{"time_dimension", "t"},
           {"units", "J m^-3 (m^-3)^(1-gamma) s^-1"},
           {"conversion", SI::qe * Tnorm * pow(Nnorm, 2. - dipole_gamma) * Omega_ci},
           {"standard_name", "dipole quasilinear entropy density source"},
           {"long_name", name + " dipole quasilinear entropy density source F_s"},
           {"species", name},
           {"source", "dipole_anomalous_transport"}});
    }

    if (flow_xlow.isAllocated()) {
      set_with_attrs(state[fmt::format("pf{}_dipole_xlow", name)], flow_xlow,
                     {{"time_dimension", "t"},
                      {"units", "s^-1"},
                      {"conversion", rho_s0 * SQ(rho_s0) * Nnorm * Omega_ci},
                      {"standard_name", "dipole particle flow"},
                      {"long_name", name + " dipole transport particle flow in X. "},
                      {"species", name},
                      {"source", "dipole_anomalous_transport"}});
    }
    if (flow_ylow.isAllocated()) {
      set_with_attrs(state[fmt::format("pf{}_dipole_ylow", name)], flow_ylow,
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
