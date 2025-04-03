
#include <bout/constants.hxx>
#include <bout/fv_ops.hxx>
#include <bout/field_factory.hxx>
#include <bout/output_bout_types.hxx>
#include <bout/derivs.hxx>
#include <bout/difops.hxx>
#include <bout/initialprofiles.hxx>

#include "../include/div_ops.hxx"
#include "../include/evolve_density.hxx"
#include "../include/hermes_utils.hxx"
#include "../include/hermes_build_config.hxx"

using bout::globals::mesh;

EvolveDensity::EvolveDensity(std::string name, Options& alloptions, Solver* solver)
    : name(name) {
  AUTO_TRACE();

  auto& options = alloptions[name];

  bndry_flux = options["bndry_flux"]
                   .doc("Allow flows through radial boundaries")
                   .withDefault<bool>(true);

  poloidal_flows =
      options["poloidal_flows"].doc("Include poloidal ExB flow").withDefault<bool>(true);

  density_floor = options["density_floor"].doc("Minimum density floor").withDefault(1e-7);

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

  if (evolve_log) {
    // Evolve logarithm of density
    solver->add(logN, std::string("logN") + name);
    // Save the density to the restart file
    // so the simulation can be restarted evolving density
    // get_restart_datafile()->addOnce(N, std::string("N") + name);

    if (!alloptions["hermes"]["restarting"]) {
      // Set logN from N input options
      initial_profile(std::string("N") + name, N);
      logN = log(N);
    } else {
      // Ignore these settings
      Options::root()[std::string("N") + name].setConditionallyUsed();
    }
  } else {
    // Evolve the density in time
    solver->add(N, std::string("N") + name);
  }

  // Charge and mass
  charge = options["charge"].doc("Particle charge. electrons = -1");
  AA = options["AA"].doc("Particle atomic mass. Proton = 1");

  diagnose =
      options["diagnose"].doc("Output additional diagnostics?").withDefault<bool>(false);

  const Options& units = alloptions["units"];
  const BoutReal Nnorm = units["inv_meters_cubed"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();

  auto& n_options = alloptions[std::string("N") + name];
  source_time_dependent = n_options["source_time_dependent"]
    .doc("Use a time-dependent source?")
    .withDefault<bool>(false);

  source_only_in_core = n_options["source_only_in_core"]
    .doc("Zero the source outside the closed field-line region?")
    .withDefault<bool>(false);

  source_normalisation = Nnorm * Omega_ci;
  time_normalisation = 1./Omega_ci;

  // Try to read the density source from the mesh
  // Units of particles per cubic meter per second
  source = 0.0;
  mesh->get(source, std::string("N") + name + "_src");
  // Allow the user to override the source from input file
  source = n_options["source"]
    .doc("Source term in ddt(N" + name + std::string("). Units [m^-3/s]"))
    .withDefault(source)
    / source_normalisation;

  // If time dependent, parse the function with respect to time from the input file
  if (source_time_dependent) {
    auto str = n_options["source_prefactor"]
      .doc("Time-dependent function of multiplier on ddt(N" + name + std::string(") source."))
      .as<std::string>();
      source_prefactor_function = FieldFactory::get()->parse(str, &n_options);
  }

  // Putting source at first X index would put it in both PFR in core, this ensures only core
  if (source_only_in_core) {
    for (int x = mesh->xstart; x <= mesh->xend; x++) {
      if (!mesh->periodicY(x)) {
        // Not periodic, so not in core
        for (int y = mesh->ystart; y <= mesh->yend; y++) {
          for (int z = mesh->zstart; z <= mesh->zend; z++) {
            source(x, y, z) = 0.0;
          }
        }
      }
    }
  }

  neumann_boundary_average_z = alloptions[std::string("N") + name]["neumann_boundary_average_z"]
    .doc("Apply neumann boundary with Z average?")
    .withDefault<bool>(false);
}

void EvolveDensity::transform(Options& state) {
  AUTO_TRACE();

  if (evolve_log) {
    // Evolving logN, but most calculations use N
    N = exp(logN);
  }

  mesh->communicate(N);

  if (neumann_boundary_average_z) {
    // Take Z (usually toroidal) average and apply as X (radial) boundary condition
    if (mesh->firstX()) {
      for (int j = mesh->ystart; j <= mesh->yend; j++) {
        BoutReal Navg = 0.0; // Average N in Z
        for (int k = 0; k < mesh->LocalNz; k++) {
          Navg += N(mesh->xstart, j, k);
        }
        Navg /= mesh->LocalNz;

        // Apply boundary condition
        for (int k = 0; k < mesh->LocalNz; k++) {
          N(mesh->xstart - 1, j, k) = 2. * Navg - N(mesh->xstart, j, k);
          N(mesh->xstart - 2, j, k) = N(mesh->xstart - 1, j, k);
        }
      }
    }

    if (mesh->lastX()) {
      for (int j = mesh->ystart; j <= mesh->yend; j++) {
        BoutReal Navg = 0.0; // Average N in Z
        for (int k = 0; k < mesh->LocalNz; k++) {
          Navg += N(mesh->xend, j, k);
        }
        Navg /= mesh->LocalNz;

        for (int k = 0; k < mesh->LocalNz; k++) {
          N(mesh->xend + 1, j, k) = 2. * Navg - N(mesh->xend, j, k);
          N(mesh->xend + 2, j, k) = N(mesh->xend + 1, j, k);
        }
      }
    }
  }

  auto& species = state["species"][name];
  set(species["density"], floor(N, 0.0)); // Density in state always >= 0
  set(species["AA"], AA);                 // Atomic mass
  if (charge != 0.0) {                    // Don't set charge for neutral species
    set(species["charge"], charge);
  }

  if (low_n_diffuse) {
    // Calculate a diffusion coefficient which can be used in N, P and NV equations

    auto* coord = mesh->getCoordinates();

    Field3D low_n_coeff =
        SQ(coord->dy) * coord->g_22
        * log(density_floor / clamp(N, 1e-3 * density_floor, density_floor));
    low_n_coeff.applyBoundary("neumann");
    set(species["low_n_coeff"], low_n_coeff);
  }

  // The particle source needs to be known in other components
  // (e.g when electromagnetic terms are enabled)
  // So evaluate them here rather than in finally()
  if (source_time_dependent) {
    // Evaluate the source_prefactor function at the current time in seconds and scale source with it
    BoutReal time = get<BoutReal>(state["time"]);
    BoutReal source_prefactor = source_prefactor_function ->generate(bout::generator::Context().set("x",0,"y",0,"z",0,"t",time*time_normalisation));
    final_source = source * source_prefactor;
  } else {
    final_source = source;
  }
  final_source.allocate(); // Ensure unique memory storage.
  add(species["density_source"], final_source);
}

void EvolveDensity::finally(const Options& state) {
  AUTO_TRACE();

  auto& species = state["species"][name];

  // Get density boundary conditions
  // but retain densities which fall below zero
  N.setBoundaryTo(get<Field3D>(species["density"]));

  if ((fabs(charge) > 1e-5) and state.isSection("fields") and state["fields"].isSet("phi")) {
    // Electrostatic potential set and species is charged -> include ExB flow

    Field3D phi = get<Field3D>(state["fields"]["phi"]);

    ddt(N) = -Div_n_bxGrad_f_B_XPPM(N, phi, bndry_flux, poloidal_flows,
                                    true); // ExB drift
  } else {
    ddt(N) = 0.0;
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

    ddt(N) -= FV::Div_par_mod<hermes::Limiter>(N, V, fastest_wave, flow_ylow);

    if (state.isSection("fields") and state["fields"].isSet("Apar_flutter")) {
      // Magnetic flutter term
      const Field3D Apar_flutter = get<Field3D>(state["fields"]["Apar_flutter"]);
      // Note: Using -Apar_flutter rather than reversing sign in front,
      //       so that upwinding is handled correctly
      ddt(N) -= Div_n_g_bxGrad_f_B_XZ(N, V, -Apar_flutter);
    }
  }

  if (low_n_diffuse) {
    // Diffusion which kicks in at very low density, in order to
    // help prevent negative density regions

    Field3D low_n_coeff = get<Field3D>(species["low_n_coeff"]);
    ddt(N) += FV::Div_par_K_Grad_par(low_n_coeff, N);
  }

  if (low_n_diffuse_perp) {
    ddt(N) += Div_Perp_Lap_FV_Index(density_floor / floor(N, 1e-3 * density_floor), N,
                                    bndry_flux);
  }

  if (low_p_diffuse_perp) {
    Field3D Plim = floor(get<Field3D>(species["pressure"]), 1e-3 * pressure_floor);
    ddt(N) += Div_Perp_Lap_FV_Index(pressure_floor / Plim, N, true);
  }

  if (hyper_z > 0.) {
    auto* coord = N.getCoordinates();
    ddt(N) -= hyper_z * SQ(SQ(coord->dz)) * D4DZ4(N);
  }

  // Collect the external source from above with all the sources from
  // elsewhere (collisions, reactions, etc) for diagnostics
  Sn = get<Field3D>(species["density_source"]);
  ddt(N) += Sn;

  // Scale time derivatives
  if (state.isSet("scale_timederivs")) {
    ddt(N) *= get<Field3D>(state["scale_timederivs"]);
  }

  if (evolve_log) {
    ddt(logN) = ddt(N) / N;
  }

#if CHECKLEVEL >= 1
  for (auto& i : N.getRegion("RGN_NOBNDRY")) {
    if (!std::isfinite(ddt(N)[i])) {
      throw BoutException("ddt(N{}) non-finite at {}. Sn={}\n", name, i, Sn[i]);
    }
  }
#endif

  if (diagnose) {
    // Save flows if they are set

    if (species.isSet("particle_flow_xlow")) {
      flow_xlow = get<Field3D>(species["particle_flow_xlow"]);
    }
    if (species.isSet("particle_flow_ylow")) {
      flow_ylow += get<Field3D>(species["particle_flow_ylow"]);
    }
  }
}

void EvolveDensity::outputVars(Options& state) {
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);

  if (evolve_log) {
    // Save density to output files
    state[std::string("N") + name].force(N);
  }
  state[std::string("N") + name].setAttributes({{"time_dimension", "t"},
                                                {"units", "m^-3"},
                                                {"conversion", Nnorm},
                                                {"standard_name", "density"},
                                                {"long_name", name + " number density"},
                                                {"species", name},
                                                {"source", "evolve_density"}});

  if (diagnose) {
    set_with_attrs(
        state[std::string("ddt(N") + name + std::string(")")], ddt(N),
        {{"time_dimension", "t"},
         {"units", "m^-3 s^-1"},
         {"conversion", Nnorm * Omega_ci},
         {"long_name", std::string("Rate of change of ") + name + " number density"},
         {"species", name},
         {"source", "evolve_density"}});

    set_with_attrs(state[std::string("SN") + name], Sn,
                   {{"time_dimension", "t"},
                    {"units", "m^-3 s^-1"},
                    {"conversion", Nnorm * Omega_ci},
                    {"standard_name", "total density source"},
                    {"long_name", name + " total number density source"},
                    {"species", name},
                    {"source", "evolve_density"}});

    set_with_attrs(state[std::string("S") + name + std::string("_src")], final_source,
                   {{"time_dimension", "t"},
                    {"units", "m^-3 s^-1"},
                    {"conversion", Nnorm * Omega_ci},
                    {"standard_name", "external density source"},
                    {"long_name", name + " external number density source"},
                    {"species", name},
                    {"source", "evolve_density"}});

    // If fluxes have been set then add them to the output
    auto rho_s0 = get<BoutReal>(state["rho_s0"]);

    if (flow_xlow.isAllocated()) {
      set_with_attrs(state[fmt::format("pf{}_tot_xlow", name)], flow_xlow,
                   {{"time_dimension", "t"},
                    {"units", "s^-1"},
                    {"conversion", rho_s0 * SQ(rho_s0) * Nnorm * Omega_ci},
                    {"standard_name", "particle flow"},
                    {"long_name", name + " particle flow in X. Note: May be incomplete."},
                    {"species", name},
                    {"source", "evolve_density"}});
    }
    if (flow_ylow.isAllocated()) {
      set_with_attrs(state[fmt::format("pf{}_tot_ylow", name)], flow_ylow,
                   {{"time_dimension", "t"},
                    {"units", "s^-1"},
                    {"conversion", rho_s0 * SQ(rho_s0) * Nnorm * Omega_ci},
                    {"standard_name", "particle flow"},
                    {"long_name", name + " particle flow in Y. Note: May be incomplete."},
                    {"species", name},
                    {"source", "evolve_density"}});
    }
  }
}

