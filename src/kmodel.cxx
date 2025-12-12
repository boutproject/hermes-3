
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

  auto coord = mesh->getCoordinates();

  Bxy = coord->Bxy;
  
  if (k.isFci()) {
    dagp = FCI::getDagp_fv(alloptions, mesh);

    const auto coord = mesh->getCoordinates();
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

  Field3D klim = floor(k, 1e-8);

    
  const Options& species = state["species"]["e"];
    
  if (species.isSet("pressure")) {

    Field3D P = GET_NOBOUNDARY(Field3D, species["pressure"]);

    Field3D c_s = GET_NOBOUNDARY(Field3D, state["sound_speed"]);

    Field3D inv_sq_c_s = 1.0 / (c_s * c_s);

    Field3D gradPgradB = (DDX(P)/P) * (DDX(Bxy)/Bxy) / coord->g_11 + (DDZ(P)/P) * (DDZ(Bxy)/Bxy) / coord->g_33 ;

    gamma = c_s * sqrt(max(gradPgradB, 0.0));

    S_k = gamma * klim;

    D_k = R_major * klim / c_s;

    D_k.applyBoundary("neumann");
    mesh->communicate(D_k);
    D_k.applyParallelBoundary("parallel_neumann_o1");
    
    alpha = R_major * L_par * gamma * inv_sq_c_s / (lambda_q*lambda_q);

    P_k = alpha * k * k;

    set(fields["k"], k);
    set(fields["P_k"], P_k);
    set(fields["S_k"], S_k);
    set(fields["D_k"], D_k);
    
  }
    
  
  
}

void Kmodel::finally(const Options& state) {


  ddt(k) = 0.0;

  // Loop through all species

  const Options& allspecies = state["species"];
    
  
  Field3D klim = floor(k, 1e-8);

    
  auto& species = state["species"]["h+"];

  if (species.isSet("velocity")) {
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
    ddt(k) += (*dagp)(D_k, klim,dummy_A, dummy_B, false);
      
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
    
  }
  
}
