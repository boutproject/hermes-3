
#include "../include/kmodel.hxx"
#include "../include/div_ops.hxx"

#include <bout/constants.hxx>
#include <bout/derivs.hxx>
#include <bout/difops.hxx>
#include <bout/fv_ops.hxx>
#include <bout/invert/laplacexy.hxx>
#include <bout/invert_laplace.hxx>
#include <bout/version.hxx>
#include <bout/yboundary_regions.hxx>

#include "../include/div_ops.hxx"
#include "../include/hermes_utils.hxx"
#include "../include/hermes_build_config.hxx"


using bout::globals::mesh;

namespace {
BoutReal floor(BoutReal value, BoutReal min) {
  if (value < min)
    return min;
  return value;
}

}

Kmodel::Kmodel(std::string name, Options& alloptions, Solver* solver) {
  AUTO_TRACE();

  solver->add(k, "k");

  auto& options = alloptions[name];
  // Normalisations
  const Options& units = alloptions["units"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();
  const BoutReal Bnorm = units["Tesla"];
  const BoutReal Lnorm = units["meters"];

  const auto coord = mesh->getCoordinates();

  Bxy = coord->Bxy;
  Bxy.applyBoundary("neumann");
  mesh->communicate(Bxy);
  Bxy.applyParallelBoundary("parallel_neumann_o1");
  
  if (k.isFci()) {
    dagp = FCI::getDagp_fv(alloptions, mesh);

    // Note: This is 1 for a Clebsch coordinate system
    //       Remove parallel slices before operations
    bracket_factor = sqrt(coord->g_22.withoutParallelSlices()) / (coord->J.withoutParallelSlices() * coord->Bxy);
  } else {
    bracket_factor = 1.0;
  }
  diagnose = options["diagnose"].doc("Diagnostic output").withDefault<bool>(false);

  dissipative = options["dissipative"].doc("Use dissipative parallel flow with Lax flux").withDefault<bool>(false);

  R_major = options["R_major"].doc("Major radius in meter.").withDefault<BoutReal>(1.0)/Lnorm;
  R_minor = options["R_minor"].doc("Minor radius in meter.").withDefault<BoutReal>(1.0)/Lnorm;

  L_par = options["L_par"].doc("Parallel connection length in meter.").withDefault<BoutReal>(100.0)/Lnorm;

  lambda_q = options["lambda_q"].doc("Heat flux fall off length in meter.").withDefault<BoutReal>(0.008)/Lnorm;

  nu_factor = options["nu_factor"].doc("Heat flux fall off length in meter.").withDefault<BoutReal>(1.0);
  
  chi_factor = options["chi_factor"].doc("Heat flux fall off length in meter.").withDefault<BoutReal>(3.0);

  
  diffusion = options["diffusion"].doc("Use dissipative parallel flow with Lax flux").withDefault<bool>(true);

  advection = options["advection"].doc("Use dissipative parallel flow with Lax flux").withDefault<bool>(true);

  propagate = options["propagate"].doc("Use dissipative parallel flow with Lax flux").withDefault<bool>(false);
  
}

void Kmodel::transform(Options& state) {
  AUTO_TRACE();

  auto& fields = state["fields"];
  Options& allspecies = state["species"];

  bool upwind = false;

  AUTO_TRACE();
  auto coord = mesh->getCoordinates();

  k.applyBoundary();
  mesh->communicate(k);
  k.applyParallelBoundary();

  Field3D klim = floor(k, 1e-10);

  auto& species = state["species"]["h+"];

  if (species.isSet("pressure")) {


    Field3D P = get<Field3D>(species["pressure"]);

    Field3D c_s = get<Field3D>(state["sound_speed"]);

    Field3D inv_sq_c_s = 1.0 / (c_s * c_s);

    gradPgradB_X = max( (DDX(P)/P) * (DDX(Bxy)/Bxy) / coord->g_11, 0.0);
    gradPgradB_Z = max((DDZ(P)/P) * (DDZ(Bxy)/Bxy) / coord->g_33, 0.0);

    DDX_P = DDX(P) / sqrt(coord->g_11);
    DDX_B = DDX(Bxy) / sqrt(coord->g_11);

    DDZ_P = DDZ(P) / sqrt(coord->g_33);
    DDZ_B = DDZ(Bxy) / sqrt(coord->g_33);
    gamma = 0.0;
    BOUT_FOR(i, P.getRegion("RGN_NOY")) {
      BoutReal grads = (DDX_P[i] / P[i]) * (DDX_B[i] / Bxy[i]) + (DDZ_P[i] / P[i]) * (DDZ_B[i] / Bxy[i]);
      gamma[i] = c_s[i] * sqrt(std::max(grads, 0.0));
    }

    S_k = gamma * klim;

    D_k = R_major * klim / c_s;

    D_k.applyBoundary("neumann");
    mesh->communicate(D_k);
    D_k.applyParallelBoundary("parallel_neumann_o1");

    alpha = R_major * L_par * gamma * inv_sq_c_s / (lambda_q*lambda_q);

    P_k = alpha * k * k;
  }


  if (propagate) {
    for (auto& kv : allspecies.getChildren()) {


      Options& spec = allspecies[kv.first]; // Note: need non-cons                                                                                                                                                                                                                

      Field3D N = GET_NOBOUNDARY(Field3D, spec["density"]);

      Field3D T = spec.isSet("temperature")
        ? GET_NOBOUNDARY(Field3D, spec["temperature"])
        : 0.0;
      Field3D V = spec.isSet("velocity") ?
        GET_NOBOUNDARY(Field3D, spec["velocity"])
        : 0.0;

      auto AA = get<BoutReal>(spec["AA"]);


      Field3D flow_den_xlow, flow_den_zlow;
      Field3D density_source = (*dagp)(D_k, N, flow_den_xlow, flow_den_zlow, false);
      add(spec["density_source"], density_source);


      Field3D flow_mom_xlow, flow_mom_zlow;
      Field3D momentum_source = (*dagp)(AA * V * D_k, N, flow_mom_xlow, flow_mom_zlow, false) + (*dagp)(D_k * nu_factor * AA * N, V, flow_mom_xlow, flow_mom_zlow, false);
      add(spec["momentum_source"], momentum_source);


      Field3D flow_ene_xlow, flow_ene_zlow;
      Field3D energy_source = (*dagp)((3. / 2) * T * D_k, N, flow_ene_xlow, flow_ene_zlow, false) + (*dagp)(chi_factor * N * D_k, T, flow_ene_xlow, flow_ene_zlow, false);
      add(spec["energy_source"], energy_source);


    }
  } // End propagate
  
  
}

void Kmodel::finally(const Options& state) {

  const auto& fields = state["fields"];
  const Options& allspecies = state["species"];
  const auto coord = mesh->getCoordinates();


  ddt(k) = 0.0;

  
  Field3D klim = floor(k, 1e-12);

    
  auto& species = state["species"]["h+"];
  
  if (species.isSet("velocity") && advection) {
    // Parallel velocity set
    Field3D V = get<Field3D>(species["velocity"]);
    
    Field3D fastest_wave;
    if (state.isSet("fastest_wave")) {
      fastest_wave = get<Field3D>(state["fastest_wave"]);
    } else {
      Field3D T = get<Field3D>(species["temperature"]);
      BoutReal AA = get<BoutReal>(species["AA"]);
      fastest_wave = sqrt(T / AA);
    }
      
    Field3D flow_ylow;
    ddt(k) -= FV::Div_par_mod<hermes::Limiter>(klim, V, fastest_wave, flow_ylow, false, dissipative);
  }

  if (species.isSet("pressure")) {

    
    ddt(k) += S_k;
    ddt(k) -= P_k;

    Field3D dummy_A, dummy_B;
    if (diffusion) {
      ddt(k) += (*dagp)(D_k, klim,dummy_A, dummy_B, false);
    }
    
  }   

}


void Kmodel::outputVars(Options& state) {
  AUTO_TRACE();
  // Normalisations

  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);


  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);

  state["k"].setAttributes({{"time_dimension", "t"},
                               {"units", "m^2 s^-2"},
                               {"long_name", "k"},
                               {"source", "kmodel"}});

  set_with_attrs(state["D_k"], D_k,
                 {{"time_dimension", "t"},
                  {"units", "m^2 s^-1"},
                  {"source", "kmodel"}});
  
  if (diagnose) {
    set_with_attrs(state["S_k"], S_k,
                 {{"time_dimension", "t"},
                  {"units", "m^2 s^-3"},
                  {"source", "kmodel"}});

    set_with_attrs(state["P_k"], P_k,
                 {{"time_dimension", "t"},
                  {"units", "m^2 s^-3"},
                  {"source", "kmodel"}});
    
    set_with_attrs(state["gradPgradB_X"], gradPgradB_X,
                 {{"time_dimension", "t"},
                  {"units", "m^2 s^-3"},
                  {"source", "kmodel"}});

    set_with_attrs(state["gradPgradB_Z"], gradPgradB_Z,
                 {{"time_dimension", "t"},
                  {"units", "m^2 s^-3"},
                  {"source", "kmodel"}});
    
    set_with_attrs(state["gamma"], gamma,
                 {{"time_dimension", "t"},
                  {"units", "m^2 s^-3"},
                  {"source", "kmodel"}});

    set_with_attrs(state["DDX_P"], DDX_P,
                 {{"time_dimension", "t"},
                  {"units", "m^2 s^-3"},
                  {"source", "kmodel"}});

    set_with_attrs(state["DDX_B"], DDX_B,
                 {{"time_dimension", "t"},
                  {"units", "m^2 s^-3"},
                  {"source", "kmodel"}});
    
  }
  
}
