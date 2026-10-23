#include "../include/dipole_anomalous_transport.hxx"

#include "../include/div_ops.hxx"
#include <bout/constants.hxx>
#include <bout/fv_ops.hxx>
#include <bout/output_bout_types.hxx>
using bout::globals::mesh;

namespace {
/// Radial (1/J) d/dx (C df/dx), with the coefficient C = a J g11 given per cell, and
/// its flows through the x faces (y flows are zero: B = Bp, orthogonal grid).
/// Unlike Div_a_Grad_perp(_upwind)_flows, the face value of the product C is taken
/// (central average or upwind cell). For the quasilinear terms C = Dpsi(psi) is a flux
/// function, so the face coefficient is the same at every y: the flux summed over a
/// flux surface is exactly C d(n U)/dpsi and the zero net flux state is n U = const,
/// with or without upwinding
Field3D Div_x_C_Grad_flows(const Field2D& C, const Field2D& f, Field3D& flow_xlow,
                           Field3D& flow_ylow, bool upwind) {
  const Coordinates* coord = mesh->getCoordinates();
  Field3D result{zeroFrom(Field3D{f})};
  flow_xlow = zeroFrom(result);
  flow_ylow = zeroFrom(result);
  for (int i = mesh->xstart - 1; i <= mesh->xend; i++) {
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      const BoutReal grad =
          2. * (f(i + 1, j) - f(i, j)) / (coord->dx()(i, j) + coord->dx()(i + 1, j));
      const BoutReal C_face =
          upwind ? ((grad > 0) ? C(i + 1, j) : C(i, j)) : 0.5 * (C(i, j) + C(i + 1, j));
      const BoutReal fout = C_face * grad; // Flux from i to i+1 is -fout
      for (int k = 0; k < mesh->LocalNz; k++) {
        result(i, j, k) += fout / (coord->dx()(i, j) * coord->J()(i, j));
        result(i + 1, j, k) -= fout / (coord->dx()(i + 1, j) * coord->J()(i + 1, j));
        // Flow is positive in the positive x direction
        flow_xlow(i + 1, j, k) = -fout * coord->dy()(i, j) * coord->dz()(i, j);
      }
    }
  }
  return result;
}
} // namespace

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

  // No quasilinear flux at all through the inner face, whatever the n, p boundary
  // conditions
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

  // The quasilinear flux -C d(n htheta/B) = -C [(htheta/B) dn + n d(htheta/B)] has a
  // geometric (convective, inward pinch) part. Without it at the inner face only
  // -C (htheta/B) dn remains: neumann on n then gives zero flux at the magnet.
  // Same for the entropy s = p n^(1-gamma) (neumann on n and p: zero heat flux)
  const bool dipole_inner_noconvection =
      options["dipole_inner_noconvection"]
          .doc("Default for dipole_inner_noconvection_particle and _heat")
          .withDefault<bool>(false);
  dipole_inner_noconvection_particle =
      options["dipole_inner_noconvection_particle"]
          .doc("Drop the geometric part of the quasilinear particle (and momentum) flux "
               "at the inner boundary, so the density boundary condition sets it?")
          .withDefault<bool>(dipole_inner_noconvection);
  dipole_inner_noconvection_heat =
      options["dipole_inner_noconvection_heat"]
          .doc("Drop the geometric part of the quasilinear entropy (heat) flux at the "
               "inner boundary, so the n and p boundary conditions set it?")
          .withDefault<bool>(dipole_inner_noconvection);

  if ((dipole_inner_noflux_particle and dipole_inner_noconvection_particle)
      or (dipole_inner_noflux_heat and dipole_inner_noconvection_heat)) {
    throw BoutException("{}: dipole_inner_noflux and dipole_inner_noconvection are two "
                        "different inner boundary treatments; choose one per channel "
                        "(particle, heat)",
                        name);
  }

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

  // Optional relaxation in time from D_initial to dipole_quasilinear_D, e.g. to start
  // with a large D (fast, robust relaxation) and reach the target D smoothly
  D_relax_time = options["dipole_quasilinear_D_relax_time"]
                     .doc("Relaxation time [s] of the quasilinear D from "
                          "dipole_quasilinear_D_initial to dipole_quasilinear_D: D(t) = "
                          "D + (D_initial - D) exp(-t / relax_time). <= 0: off")
                     .withDefault(-1.0)
                 / alloptions["units"]["seconds"].as<BoutReal>();
  if (D_relax_time > 0) {
    dipole_quasilinear_D_initial =
        options["dipole_quasilinear_D_initial"]
            .doc("Starting value of the quasilinear D relaxation (units as "
                 "dipole_quasilinear_D)")
            .as<Field2D>()
        * model_scaling(dipole_quasilinear_model);
    setPermissions(readOnly("time"));
  }

  dipole_quasilinear_closure =
      options["dipole_quasilinear_closure"]
          .doc("Scale dipole_quasilinear_D by d/dipole_gamma, d = -dln p/dln U from the "
               "total charged pressure p(psi) and the flux tube volume U = int dl/B")
          .withDefault<bool>(false);
  if (dipole_quasilinear_closure) {
    setPermissions(readIfSet("species:{all_species}:charge"));
    setPermissions(readIfSet("species:{charged}:pressure", Regions::Interior));
  }

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
  if (D_relax_time > 0) {
    dipole_quasilinear_D_initial *= closed;
    mesh->communicate(dipole_quasilinear_D_initial);
  }

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
  // Inner guard cells of f = q htheta/B. noflux: f = f(xstart), zero face gradient.
  // noconvection: q htheta/B with the first cell's htheta/B, so the face gradient is
  // (htheta/B) dq and q's boundary condition (guard values) sets the flux
  auto inner_boundary = [&](Field2D& f, const Field2D& q, bool noflux,
                            bool noconvection) {
    if ((noflux or noconvection) and mesh->firstX()) {
      for (int j = mesh->ystart; j <= mesh->yend; j++) {
        const BoutReal ratio = htheta2D(mesh->xstart, j) / B2D(mesh->xstart, j);
        for (int i = mesh->xstart - 2; i < mesh->xstart; i++) {
          f(i, j) = noflux ? f(mesh->xstart, j) : q(i, j) * ratio;
        }
      }
    }
  };
  Field2D nthetaB = N2D * htheta2D / B2D;
  mesh->communicate(nthetaB);
  inner_boundary(nthetaB, N2D, dipole_inner_noflux_particle,
                 dipole_inner_noconvection_particle);
  const Coordinates* coord = mesh->getCoordinates();
  Field2D D_QL = DC(dipole_quasilinear_D);
  if (D_relax_time > 0) {
    // Absolute simulation time: on restart the relaxation continues where it was
    const BoutReal time = get<BoutReal>(state["time"]);
    D_QL += (DC(dipole_quasilinear_D_initial) - D_QL) * exp(-time / D_relax_time);
  }
  if (dipole_quasilinear_closure) {
    // Interchange closure D(psi) = D0 d/gamma, d = -dln p/dln U (marginal: p U^gamma
    // = const). Total charged pressure: the displacement is common to all species
    GuardedOptions allspecies = state["species"];
    Field2D Ptot{zeroFrom(B2D)};
    for (auto& kv : allspecies.getChildren()) {
      GuardedOptions sp = allspecies[kv.first];
      if (!sp.isSet("charge") || fabs(get<BoutReal>(sp["charge"])) < 1e-5
          || !IS_SET_NOBOUNDARY(sp["pressure"])) {
        continue;
      }
      Ptot += DC(GET_NOBOUNDARY(Field3D, sp["pressure"]));
    }
    // Flux surface sums: U = sum J dy (flux tube volume per unit flux), volume
    // averaged pressure. Not averageY, which assumes no branch cuts
    Field2D U{zeroFrom(B2D)}, pbar{zeroFrom(B2D)};
    for (int x = mesh->xstart; x <= mesh->xend; x++) {
      BoutReal local[2] = {0., 0.}, sum[2];
      for (int y = mesh->ystart; y <= mesh->yend; y++) {
        const BoutReal dV = coord->J()(x, y) * coord->dy()(x, y);
        local[0] += dV;
        local[1] += Ptot(x, y) * dV;
      }
      MPI_Allreduce(local, sum, 2, MPI_DOUBLE, MPI_SUM, mesh->getYcomm(x));
      for (int y = 0; y < mesh->LocalNy; y++) {
        U(x, y) = sum[0];
        pbar(x, y) = sum[1] / sum[0];
      }
    }
    mesh->communicate(U, pbar);

    // Central difference over closed surfaces, one-sided at the domain edges and
    // at the last closed surface
    closure_d = 0.0;
    for (int x = mesh->xstart; x <= mesh->xend; x++) {
      if (!mesh->periodicY(x)) {
        continue;
      }
      const int xl =
          (mesh->periodicY(x - 1) && !(mesh->firstX() && x == mesh->xstart)) ? x - 1 : x;
      const int xr =
          (mesh->periodicY(x + 1) && !(mesh->lastX() && x == mesh->xend)) ? x + 1 : x;
      const int y0 = mesh->ystart;
      const BoutReal dlnU = log(U(xr, y0) / U(xl, y0));
      if (fabs(dlnU) < 1e-12) {
        continue;
      }
      const BoutReal d =
          -log(std::max(pbar(xr, y0), 1e-8) / std::max(pbar(xl, y0), 1e-8)) / dlnU;
      for (int y = 0; y < mesh->LocalNy; y++) {
        closure_d(x, y) = d;
      }
    }
    // Boundary guard cells take the edge value, as the face coefficient uses them
    for (int y = 0; y < mesh->LocalNy; y++) {
      for (int x = 0; x < mesh->xstart; x++) {
        if (mesh->firstX()) {
          closure_d(x, y) = closure_d(mesh->xstart, y);
        }
      }
      for (int x = mesh->xend + 1; x < mesh->LocalNx; x++) {
        if (mesh->lastX()) {
          closure_d(x, y) = closure_d(mesh->xend, y);
        }
      }
    }
    // d < 0 (hollow pressure) would give anti-diffusion
    D_QL *= floor(closure_d, 0.0) / dipole_gamma;
  }
  Field2D a_QL = D_QL * B2D / htheta2D;
  // C = a J g11, a flux function (Dpsi) for model 1, so that upwinding keeps the
  // marginal state n U = const
  const Field2D JG11{coord->J() * coord->g11()};
  auto div_QL = [&](const Field2D& a, const Field2D& f, Field3D& fx, Field3D& fy) {
    return Div_x_C_Grad_flows(a * JG11, f, fx, fy, dipole_upwind);
  };
  mesh->communicate(a_QL);

  // F_n, eqs. (52), (54)
  F_n = div_QL(a_QL, nthetaB, flow_xlow, flow_ylow);
  add(species["density_source"], F_n);
  add(species["particle_flow_xlow"], flow_xlow);
  add(species["particle_flow_ylow"], flow_ylow);

  // Parallel momentum carried by the quasilinear particle flux
  add(species["momentum_source"], div_QL(AA * V2D * a_QL, nthetaB, fx, fy));
  add(species["momentum_flow_xlow"], fx);
  add(species["momentum_flow_ylow"], fy);

  // F_s, eq. (53): s = p n^(1-gamma) is conserved by the interchange displacement
  const Field2D Nlim = floor(N2D, 1e-8);
  const Field2D S2D = P2D * pow(Nlim, 1. - dipole_gamma);
  Field2D sthetaB = S2D * htheta2D / B2D;
  inner_boundary(sthetaB, S2D, dipole_inner_noflux_heat, dipole_inner_noconvection_heat);
  F_s = div_QL(a_QL, sthetaB, fx, fy);
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
    // Fixed at initialisation, so saved once without a time dimension
    set_with_attrs(state[fmt::format("dipole_anomalous_D_{}", name)], dipole_anomalous_D,
                   {{"units", "m^2 s^-1"},
                    {"conversion", rho_s0 * rho_s0 * Omega_ci},
                    {"standard_name", "dipole anomalous density diffusion"},
                    {"long_name", "Dipole anomalous density diffusion of " + name},
                    {"source", "dipole_anomalous_transport"}});
    set_with_attrs(state[fmt::format("dipole_anomalous_chi_{}", name)], dipole_anomalous_chi,
                   {{"units", "m^2 s^-1"},
                    {"conversion", rho_s0 * rho_s0 * Omega_ci},
                    {"standard_name", "dipole anomalous thermal diffusion"},
                    {"long_name", "Dipole anomalous thermal diffusion of " + name},
                    {"source", "dipole_anomalous_transport"}});

    if (dipole_quasilinear_closure) {
      set_with_attrs(state[fmt::format("dipole_d_{}", name)], closure_d,
                     {{"time_dimension", "t"},
                      {"units", ""},
                      {"conversion", 1.0},
                      {"standard_name", "dipole interchange drive"},
                      {"long_name", "d = -dln p/dln U, quasilinear D = D0 d/gamma"},
                      {"species", name},
                      {"source", "dipole_anomalous_transport"}});
    }
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
