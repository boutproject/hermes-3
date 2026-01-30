#include <bout/constants.hxx>
#include <bout/fv_ops.hxx>
#include <bout/field_factory.hxx>
#include <bout/output_bout_types.hxx>
#include <bout/derivs.hxx>
#include <bout/difops.hxx>
#include <bout/initialprofiles.hxx>

#include "../include/div_ops.hxx"
#include "../include/evolve_tracedensity.hxx"
#include "../include/hermes_utils.hxx"
#include "../include/hermes_build_config.hxx"

using bout::globals::mesh;

EvolveTraceDensity::EvolveTraceDensity(std::string name, Options& alloptions, Solver* solver)
    : name(name) {
  AUTO_TRACE();

  auto& options = alloptions[name];

  
  
  bndry_flux = options["bndry_flux"]
                   .doc("Allow flows through radial boundaries")
                   .withDefault<bool>(true);


  scale_ExB = options["scale_ExB"]
                   .doc("Scale ExB flow?")
                   .withDefault<BoutReal>(1.0);

  output_ddt = options["output_ddt"]
                   .doc("Include ExB advection?")
                   .withDefault<bool>(false);
  
  poloidal_flows =
      options["poloidal_flows"].doc("Include poloidal ExB flow").withDefault<bool>(true);

  density_floor = options["density_floor"].doc("Minimum density floor").withDefault(1e-5);

  low_n_diffuse = options["low_n_diffuse"]
                      .doc("Parallel diffusion at low density")
                      .withDefault<bool>(false);

  low_n_diffuse_perp = options["low_n_diffuse_perp"]
                           .doc("Perpendicular diffusion at low density")
                           .withDefault<bool>(false);

  pressure_floor = density_floor * (1./get<BoutReal>(alloptions["units"]["eV"]));

  low_p_diffuse_perp = options["low_p_diffuse_perp"]
                           .doc("Perpendicular diffusion at low pressure")
                           .withDefault<bool>(false);

  hyper_z = options["hyper_z"].doc("Hyper-diffusion in Z").withDefault(-1.0);

  evolve_log = options["evolve_log"]
                   .doc("Evolve the logarithm of density?")
                   .withDefault<bool>(false);

  isMMS = options["mms"].withDefault<bool>(false);

  dissipative = options["dissipative"].doc("Use dissipative parallel flow with Lax flux").withDefault<bool>(false);
  
  // Evolve the density in time
  solver->add(tN, std::string("tN") + name);
  

  // Charge and mass
  charge = options["charge"].doc("Particle charge. electrons = -1");
  AA = options["AA"].doc("Particle atomic mass. Proton = 1");

  diagnose =
      options["diagnose"].doc("Output additional diagnostics?").withDefault<bool>(false);

  const Options& units = alloptions["units"];
  const BoutReal Nnorm = units["inv_meters_cubed"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();
  const BoutReal Lnorm = units["meters"];
  
  auto& tn_options = alloptions[std::string("tN") + name];

  hyper_n = options["hyper_n"].doc("Hyper-viscosity. < 0 -> off").withDefault(-1.0) / (Lnorm * Lnorm * Lnorm * Lnorm * Omega_ci);
  
  source_normalisation = Nnorm * Omega_ci;
  time_normalisation = 1./Omega_ci;

  
  exb_advection = tn_options["exb_advection"]
                   .doc("Include ExB advection?")
                   .withDefault<bool>(options["exb_advection"].withDefault<bool>(true));
  
  sourceterm = tn_options["sourceterm"].withDefault<bool>(true);
  
  // Try to read the density source from the mesh
  // Units of particles per cubic meter per second
  source = 0.0;
  mesh->get(source, std::string("tN") + name + "_src");
  // Allow the user to override the source from input file
  source = tn_options["source"]
    .doc("Source term in ddt(N" + name + std::string("). Units [m^-3/s]"))
    .withDefault(source)
    / source_normalisation;

  disable_ddt = tn_options["disable_ddt"]
    .withDefault<bool>(false);
  

  if (mesh->isFci()) {
    const auto coord = mesh->getCoordinates();
    // Note: This is 1 for a Clebsch coordinate system
    //       Remove parallel slices before operations
    bracket_factor = sqrt(coord->g_22.withoutParallelSlices()) / (coord->J.withoutParallelSlices() * coord->Bxy);
  } else {
    // Clebsch coordinate system
    bracket_factor = 1.0;
  }

  
}

void EvolveTraceDensity::transform(Options& state) {
  AUTO_TRACE();

  tN.applyBoundary();
  mesh->communicate(tN);
  tN.applyParallelBoundary();

  auto& species = state["species"][name];
  set(species["tracedensity"], floor(tN, 0.0)); // Density in state always >= 0
                 // Atomic mass

  // The particle source needs to be known in other components
  // (e.g when electromagnetic terms are enabled)
  // So evaluate them here rather than in finally()
  final_source = source;
  

  
  
  final_source.allocate(); // Ensure unique memory storage.
  add(species["tracedensity_source"], final_source);
}

void EvolveTraceDensity::finally(const Options& state) {
  AUTO_TRACE();

  auto& species = state["species"][name];

  // Get density boundary conditions
  // but retain densities which fall below zero
  tN.setBoundaryTo(get<Field3D>(species["tracedensity"]));

  if (exb_advection and (fabs(charge) > 1e-5) and
      state.isSection("fields") and state["fields"].isSet("phi")) {
    // Electrostatic potential set and species is charged -> include ExB flow

    Field3D phi = get<Field3D>(state["fields"]["phi"]);

    ddt(tN) = -scale_ExB * Div_n_bxGrad_f_B_XPPM(tN, phi, bndry_flux, poloidal_flows,
                                    true) * bracket_factor; // ExB drift
  } else {
    ddt(tN) = 0.0;
  }

  if (species.isSet("velocity")) {
    // Parallel velocity set
    Field3D V = get<Field3D>(species["velocity"]);

    // Wave speed used for numerical diffusion
    // Note: For simulations where ion density is evolved rather than electron density,
    // the fast electron dynamics still determine the stability.
    Field3D fastest_wave;
    if (state.isSet("fastest_wave")) {
      fastest_wave = get<Field3D>(state["fastest_wave"]);
    } else {
      Field3D T = get<Field3D>(species["temperature"]);
      BoutReal AA = get<BoutReal>(species["AA"]);
      fastest_wave = sqrt(T / AA);
    }
    
    ddt(tN) -= FV::Div_par_mod<hermes::Limiter>(tN, V, fastest_wave, flow_ylow, false, dissipative);
    
    if (state.isSection("fields") and state["fields"].isSet("Apar_flutter")) {
      // Magnetic flutter term
      const Field3D Apar_flutter = get<Field3D>(state["fields"]["Apar_flutter"]);
      // Note: Using -Apar_flutter rather than reversing sign in front,
      //       so that upwinding is handled correctly
      ddt(tN) -= Div_n_g_bxGrad_f_B_XZ(tN, V, -Apar_flutter);
    }
  }

  
  // Collect the external source from above with all the sources from
  // elsewhere (collisions, reactions, etc) for diagnostics
  Sn = get<Field3D>(species["tracedensity_source"]);
  if (sourceterm) {
    ddt(tN) += Sn;
  }
  
  BOUT_FOR(i, tN.getRegion("RGN_NOY")) {
    if ((tN[i] < density_floor * 1e-2) && (ddt(tN)[i] < 0.0)) {
      ddt(tN)[i] = 0.0;
    }
  }
  
  // Scale time derivatives
  if (state.isSet("scale_timederivs")) {
    ddt(tN) *= get<Field3D>(state["scale_timederivs"]);
  }


  
  

#if CHECKLEVEL >= 1
  for (auto& i : tN.getRegion("RGN_NOBNDRY")) {
    if (!std::isfinite(ddt(tN)[i])) {
      throw BoutException("ddt(tN{}) non-finite at {}. Sn={}\n", name, i, Sn[i]);
    }
  }
#endif

  

  if (disable_ddt){
    ddt(tN) = 0.0;
  }
  
}

void EvolveTraceDensity::outputVars(Options& state) {
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);

  state[std::string("tN") + name].setAttributes({{"time_dimension", "t"},
                                                {"units", "m^-3"},
                                                {"conversion", Nnorm},
                                                {"standard_name", "trace density"},
                                                {"long_name", name + "trace number density"},
                                                {"species", name},
                                                {"source", "evolve_tracedensity"}});

  if (output_ddt || diagnose) {
    set_with_attrs(
        state[std::string("ddt(tN") + name + std::string(")")], ddt(tN),
        {{"time_dimension", "t"},
         {"units", "m^-3 s^-1"},
         {"conversion", Nnorm * Omega_ci},
         {"long_name", std::string("Rate of change of ") + name + " number density"},
         {"species", name},
         {"source", "evolve_density"}});
  }
}
