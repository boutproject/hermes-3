#include "../include/kmodel.hxx"
#include "../include/component.hxx"
#include "../include/div_ops.hxx"
#include "../include/guarded_options.hxx"
#include "../include/hermes_build_config.hxx"
#include "../include/hermes_utils.hxx"
#include "../include/permissions.hxx"

#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/constants.hxx>
#include <bout/derivs.hxx>
#include <bout/difops.hxx>
#include <bout/field.hxx>
#include <bout/field3d.hxx>
#include <bout/field_factory.hxx>
#include <bout/fv_ops.hxx>
#include <bout/globals.hxx>
#include <bout/options.hxx>
#include <bout/output.hxx>
#include <bout/solver.hxx>

#include <string>

using bout::globals::mesh;

Kmodel::Kmodel(std::string name, Options& alloptions, Solver* solver)
    : NamedComponent(name, {readIfSet("species:{all_species}:{input}"),
                            readOnly("sound_speed"), readOnly("fastest_wave"),
                            readWrite("species:{all_species}:{output}")}) {

  solver->add(k, "k");

  auto& options = alloptions[name];
  // Normalisations
  const Options& units = alloptions["units"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();
  const BoutReal Bnorm = units["Tesla"];
  const BoutReal Lnorm = units["meters"];

  diagnose = options["diagnose"].doc("Diagnostic output").withDefault<bool>(false);

  R_major =
      options["R_major"].doc("Major radius in meter.").withDefault<BoutReal>(1.6) / Lnorm;

  R_minor =
      options["R_minor"].doc("Minor radius in meter.").withDefault<BoutReal>(0.5) / Lnorm;

  L_par = options["L_par"]
              .doc("Parallel connection length in meter.")
              .withDefault<BoutReal>(100.0)
          / Lnorm;

  lambda_q = options["lambda_q"]
                 .doc("Heat flux fall off length in meter.")
                 .withDefault<BoutReal>(0.008)
             / Lnorm;

  average_AA = options["average_AA"]
                   .doc("Average atomic mass to calculate Pi_hat and Ni_hat.")
                   .withDefault<BoutReal>(1.0);

  diffusion =
      options["diffusion"].doc("Turn on diffusion of k via D_K?").withDefault<bool>(true);

  advection = options["advection"]
                  .doc("Turn on advection of k via parallel ion flow?")
                  .withDefault<bool>(true);

  substitutePermissions(
      "input", {"AA", "charge", "density", "pressure", "temperature", "velocity"});

  substitutePermissions("output", {"density_source", "momentum_source", "energy_source"});
}

void Kmodel::transform_impl(GuardedOptions& state) {

  auto fields = state["fields"];

  auto coord = mesh->getCoordinates();

  Bxy = coord->Bxy;

  k.applyBoundary();
  mesh->communicate(k);
  k.applyParallelBoundary();

  GuardedOptions allspecies = state["species"];

  Pi_hat = 0.0;
  Ni_hat = 0.0;

  for (auto& kv : allspecies.getChildren()) {
    GuardedOptions species = allspecies[kv.first]; // Note: Need non-const

    // Skip if not charged or no pressure set
    if (!(species.isSet("charge"))) {
      continue;
    }

    // Electrons for now
    auto q = get<BoutReal>(species["charge"]);
    if (q < 0.0) {
      continue;
    }

    // Skip if no pressure
    if (!species.isSet("pressure")) {
      continue;
    }

    const auto P = GET_NOBOUNDARY(Field3D, species["pressure"]);
    const auto AA = get<BoutReal>(species["AA"]);

    Pi_hat += P * AA / average_AA;

    // Skip if no density
    if (!species.isSet("density")) {
      continue;
    }

    const auto N = GET_NOBOUNDARY(Field3D, species["density"]);

    Ni_hat += N * AA / average_AA;
  }

  Pi_hat.applyBoundary("neumann");
  Ni_hat.applyBoundary("neumann");
  mesh->communicate(Pi_hat, Ni_hat);
  Pi_hat.applyParallelBoundary("parallel_neumann_o1");
  Ni_hat.applyParallelBoundary("parallel_neumann_o1");

  if (!state.isSet("sound_speed")) {
    throw BoutException("Sound speed not set but required for the k-model!");
  }

  Field3D c_s = get<Field3D>(state["sound_speed"]);

  Field3D inv_sq_c_s = 1.0 / (c_s * c_s);

  Field3D klim = floor(k, 1e-10);

  gradPgradB_X = floor((DDX(Pi_hat) / Pi_hat) * (DDX(Bxy) / Bxy) / coord->g_11, 0.0);

  if (k.isFci()) {
    gradPgradB_Z = floor((DDZ(Pi_hat) / Pi_hat) * (DDZ(Bxy) / Bxy) / coord->g_33, 0.0);
  } else {
    gradPgradB_Y = floor((DDY(Pi_hat) / Pi_hat) * (DDY(Bxy) / Bxy) / coord->g_22, 0.0);
  }

  DDX_P = DDX(Pi_hat) / sqrt(coord->g_11);
  DDX_B = DDX(Bxy) / sqrt(coord->g_11);

  if (k.isFci()) {
    DDZ_P = DDZ(Pi_hat) / sqrt(coord->g_33);
    DDZ_B = DDZ(Bxy) / sqrt(coord->g_33);
  } else {
    DDY_P = DDY(Pi_hat) / sqrt(coord->g_22);
    DDY_B = DDY(Bxy) / sqrt(coord->g_22);
  }

  gamma = 0.0;

  BOUT_FOR(i, Pi_hat.getRegion("RGN_NOY")) {
    BoutReal grads = k.isFci() ? (DDX_P[i] / Pi_hat[i]) * (DDX_B[i] / Bxy[i])
                                     + (DDZ_P[i] / Pi_hat[i]) * (DDZ_B[i] / Bxy[i])
                               : (DDX_P[i] / Pi_hat[i]) * (DDX_B[i] / Bxy[i])
                                     + (DDY_P[i] / Pi_hat[i]) * (DDY_B[i] / Bxy[i]);
    gamma[i] = c_s[i] * sqrt(std::max(grads, 0.0));
  }

  S_k = gamma * klim;

  D_k = R_major * klim / floor(c_s, 1e-10);

  D_k = softFloor(D_k, 1e-5);

  D_k.applyBoundary("neumann");
  mesh->communicate(D_k);
  D_k.applyParallelBoundary("parallel_neumann_o1");

  alpha = R_major * L_par * gamma * inv_sq_c_s / (lambda_q * lambda_q);

  P_k = alpha * klim * klim;
}

void Kmodel::finally(const Options& state) {

  ddt(k) = 0.0;

  ddt(k) += S_k - P_k;

  if (diffusion) {
    flux_k_x = 0.0;
    flux_k_y = 0.0;
    ddt(k) += Div_a_Grad_perp_upwind_flows(D_k, k, flux_k_x, flux_k_y);
  }

  if (advection) {

    for (const auto& kv : state["species"].getChildren()) {

      const Options& species = kv.second;

      // Skip if not charged
      if (!(species.isSet("charge"))) {
        continue;
      }

      // Electrons for now
      auto q = get<BoutReal>(species["charge"]);
      if (q < 0.0) {
        continue;
      }

      // Skip if no pressure, density or velocity set
      if (!species.isSet("pressure") or !species.isSet("density")
          or !species.isSet("velocity")) {
        continue;
      }

      // Typical wave speed used for numerical diffusion
      Field3D fastest_wave;
      if (state.isSet("fastest_wave")) {
        fastest_wave = get<Field3D>(state["fastest_wave"]);
      } else {
        Field3D T = get<Field3D>(species["temperature"]);
        const BoutReal AA = get<BoutReal>(species["AA"]);
        fastest_wave = sqrt(T / AA);
      }

      Field3D N = get<Field3D>(species["density"]);
      Field3D V = get<Field3D>(species["temperature"]);

      Field3DParallel V_hat = N * V / average_AA;

      Field3D dummy;
      ddt(k) -= FV::Div_par_mod<hermes::Limiter>(V_hat, k, fastest_wave, dummy);
    }
  }
}

void Kmodel::outputVars(Options& state) {

  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto Lnorm = get<BoutReal>(state["rho_s0"]);

  state["k"].setAttributes({{"time_dimension", "t"},
                            {"units", "C m^-3"},
                            {"conversion", SI::qe * Nnorm},
                            {"long_name", "Turbulent kinetic energy"},
                            {"source", "kmodel"}});

  set_with_attrs(state["D_k"], D_k,
                 {{"time_dimension", "t"},
                  {"units", "m^2 / s"},
                  {"conversion", Lnorm * Lnorm * Omega_ci},
                  {"standard_name", "Turbulent diffusion coefficient"},
                  {"long_name", "Turbulent diffusion coefficient"},
                  {"source", "kmodel"}});
}
