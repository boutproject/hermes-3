
#include <bout/constants.hxx>
#include <bout/fv_ops.hxx>
#include <bout/invert_pardiv.hxx>
#include <bout/output_bout_types.hxx>
#include <derivs.hxx>
#include <difops.hxx>
#include <initialprofiles.hxx>

#include "../include/div_ops.hxx"
#include "../include/evolve_energy.hxx"
#include "../include/hermes_utils.hxx"
#include "../include/hermes_build_config.hxx"

using bout::globals::mesh;

EvolveEnergy::EvolveEnergy(std::string name, Options& alloptions, Solver* solver)
    : name(name) {
  AUTO_TRACE();

  auto& options = alloptions[name];

  adiabatic_index =
      options["adiabatic_index"]
          .doc("Ratio of specific heats γ = Cp/Cv [5/3 for monatomic ideal gas]")
          .withDefault(5. / 3);
  Cv = 1. / (adiabatic_index - 1.);

  evolve_log = options["evolve_log"]
                   .doc("Evolve the logarithm of energy?")
                   .withDefault<bool>(false);

  density_floor = options["density_floor"].doc("Minimum density floor").withDefault(1e-5);
  if (evolve_log) {
    // Evolve logarithm of energy
    solver->add(logE, std::string("logE") + name);
    // Save the pressure to the restart file
    // so the simulation can be restarted evolving energy
    // get_restart_datafile()->addOnce(E, std::string("E") + name);

    if (!alloptions["hermes"]["restarting"]) {
      // Set logE from E input options
      initial_profile(std::string("E") + name, E);
      logE = log(E);
    } else {
      // Ignore these settings
      Options::root()[std::string("E") + name].setConditionallyUsed();
    }
  } else {
    // Evolve the energy in time
    solver->add(E, std::string("E") + name);
  }

  bndry_flux = options["bndry_flux"]
                   .doc("Allow flows through radial boundaries")
                   .withDefault<bool>(true);

  poloidal_flows =
      options["poloidal_flows"].doc("Include poloidal ExB flow").withDefault<bool>(true);

  thermal_conduction = options["thermal_conduction"]
                           .doc("Include parallel heat conduction?")
                           .withDefault<bool>(true);

  kappa_coefficient = options["kappa_coefficient"]
                          .doc("Numerical coefficient in parallel heat conduction. "
                               "Default is 3.16 for electrons, 3.9 otherwise")
                          .withDefault((name == "e") ? 3.16 : 3.9);

  kappa_limit_alpha = options["kappa_limit_alpha"]
                          .doc("Flux limiter factor. < 0 means no limit. Typical is 0.2 "
                               "for electrons, 1 for ions.")
                          .withDefault(-1.0);

  hyper_z = options["hyper_z"].doc("Hyper-diffusion in Z").withDefault(-1.0);

  diagnose = options["diagnose"]
                 .doc("Save additional output diagnostics")
                 .withDefault<bool>(false);

  enable_precon = options["precondition"]
                      .doc("Enable preconditioner? (Note: solver may not use it)")
                      .withDefault<bool>(true);

  const Options& units = alloptions["units"];
  const BoutReal Nnorm = units["inv_meters_cubed"];
  const BoutReal Tnorm = units["eV"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();

  // Try to read the pressure source from the mesh
  source = 0.0;
  mesh->get(source, std::string("P") + name + "_src"); // Units of Pascals per second
  source *= Cv;                                        // Convert to W/m^3

  // Allow the user to override the source
  source = alloptions[std::string("E") + name]["source"]
               .doc(std::string("Source term in ddt(E") + name
                    + std::string("). Units [W/m^3]"))
               .withDefault(source)
           / (SI::qe * Nnorm * Tnorm * Omega_ci);

  if (alloptions[std::string("E") + name]["source_only_in_core"]
          .doc("Zero the source outside the closed field-line region?")
          .withDefault<bool>(false)) {
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

  neumann_boundary_average_z =
      alloptions[std::string("E") + name]["neumann_boundary_average_z"]
          .doc("Apply neumann boundary with Z average?")
          .withDefault<bool>(false);
}

void EvolveEnergy::transform(Options& state) {
  AUTO_TRACE();

  if (evolve_log) {
    // Evolving logE, but most calculations use E
    E = exp(logE);
  }

  mesh->communicate(E);

  auto& species = state["species"][name];
  N = getNoBoundary<Field3D>(species["density"]);
  const Field3D V = getNoBoundary<Field3D>(species["velocity"]);
  const BoutReal AA = get<BoutReal>(species["AA"]);

  // Calculate pressure
  // E = Cv * P + (1/2) m n v^2
  P.allocate();
  BOUT_FOR(i, P.getRegion("RGN_NOBNDRY")) {
    P[i] = (E[i] - 0.5 * AA * N[i] * SQ(V[i])) / Cv;
    if (P[i] < 0.0) {
      P[i] = 0.0;
    }
  }

  if (neumann_boundary_average_z) {
    // Take Z (usually toroidal) average and apply as X (radial) boundary condition
    if (mesh->firstX()) {
      for (int j = mesh->ystart; j <= mesh->yend; j++) {
        BoutReal Pavg = 0.0; // Average P in Z
        for (int k = 0; k < mesh->LocalNz; k++) {
          Pavg += P(mesh->xstart, j, k);
        }
        Pavg /= mesh->LocalNz;

        // Apply boundary condition
        for (int k = 0; k < mesh->LocalNz; k++) {
          P(mesh->xstart - 1, j, k) = 2. * Pavg - P(mesh->xstart, j, k);
          P(mesh->xstart - 2, j, k) = P(mesh->xstart - 1, j, k);
        }
      }
    }

    if (mesh->lastX()) {
      for (int j = mesh->ystart; j <= mesh->yend; j++) {
        BoutReal Pavg = 0.0; // Average P in Z
        for (int k = 0; k < mesh->LocalNz; k++) {
          Pavg += P(mesh->xend, j, k);
        }
        Pavg /= mesh->LocalNz;

        for (int k = 0; k < mesh->LocalNz; k++) {
          P(mesh->xend + 1, j, k) = 2. * Pavg - P(mesh->xend, j, k);
          P(mesh->xend + 2, j, k) = P(mesh->xend + 1, j, k);
        }
      }
    }
  }

  // Calculate temperature
  T = P / floor(N, density_floor);
  P = N * T; // Ensure consistency

  set(species["pressure"], P);
  set(species["temperature"], T);
}

void EvolveEnergy::finally(const Options& state) {
  AUTO_TRACE();

  /// Get the section containing this species
  const auto& species = state["species"][name];

  // Get updated pressure and temperature with boundary conditions
  P = get<Field3D>(species["pressure"]);
  T = get<Field3D>(species["temperature"]);
  N = get<Field3D>(species["density"]);
  const Field3D V = get<Field3D>(species["velocity"]);

  Field3D Pfloor = P;

  if (state.isSection("fields") and state["fields"].isSet("phi")) {
    // Electrostatic potential set -> include ExB flow

    Field3D phi = get<Field3D>(state["fields"]["phi"]);

    ddt(E) = -Div_n_bxGrad_f_B_XPPM(E, phi, bndry_flux, poloidal_flows, true);
  } else {
    ddt(E) = 0.0;
  }

  if (species.isSet("velocity")) {
    Field3D V = get<Field3D>(species["velocity"]);

    // Typical wave speed used for numerical diffusion
    Field3D fastest_wave;
    if (state.isSet("fastest_wave")) {
      fastest_wave = get<Field3D>(state["fastest_wave"]);
    } else {
      const BoutReal AA = get<BoutReal>(species["AA"]);
      fastest_wave = sqrt(T / AA);
    }

    ddt(E) -= FV::Div_par<hermes::Limiter>(E + P, V, fastest_wave);
  }

  if (species.isSet("low_n_coeff")) {
    // Low density parallel diffusion
    Field3D low_n_coeff = get<Field3D>(species["low_n_coeff"]);
    ddt(E) += FV::Div_par_K_Grad_par(low_n_coeff * T, N)
              + FV::Div_par_K_Grad_par(low_n_coeff, P);
  }

  // Parallel heat conduction
  if (thermal_conduction) {

    // Calculate ion collision times
    const Field3D tau = 1. / floor(get<Field3D>(species["collision_frequency"]), 1e-10);
    const BoutReal AA = get<BoutReal>(species["AA"]); // Atomic mass

    // Parallel heat conduction
    // Braginskii expression for parallel conduction
    // kappa ~ n * v_th^2 * tau
    //
    // Note: Coefficient is slightly different for electrons (3.16) and ions (3.9)
    kappa_par = kappa_coefficient * Pfloor * tau / AA;

    if (kappa_limit_alpha > 0.0) {
      /*
       * Flux limiter, as used in SOLPS.
       *
       * Calculate the heat flux from Spitzer-Harm and flux limit
       *
       * Typical value of alpha ~ 0.2 for electrons
       *
       * R.Schneider et al. Contrib. Plasma Phys. 46, No. 1-2, 3 – 191 (2006)
       * DOI 10.1002/ctpp.200610001
       */

      // Spitzer-Harm heat flux
      Field3D q_SH = kappa_par * Grad_par(T);
      // Free-streaming flux
      Field3D q_fl = kappa_limit_alpha * N * T * sqrt(T / AA);

      // This results in a harmonic average of the heat fluxes
      kappa_par = kappa_par / (1. + abs(q_SH / floor(q_fl, 1e-10)));

      // Values of kappa on cell boundaries are needed for fluxes
      mesh->communicate(kappa_par);
    }

    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        auto i = indexAt(kappa_par, r.ind, mesh->ystart, jz);
        auto im = i.ym();
        kappa_par[im] = kappa_par[i];
      }
    }
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        auto i = indexAt(kappa_par, r.ind, mesh->yend, jz);
        auto ip = i.yp();
        kappa_par[ip] = kappa_par[i];
      }
    }

    // Note: Flux through boundary turned off, because sheath heat flux
    // is calculated and removed separately
    ddt(E) += FV::Div_par_K_Grad_par(kappa_par, T, false);
  }

  if (hyper_z > 0.) {
    auto* coord = N.getCoordinates();
    ddt(E) -= hyper_z * SQ(SQ(coord->dz)) * D4DZ4(E);
  }

  //////////////////////
  // Other sources

  Se = source;
  if (species.isSet("energy_source")) {
    Se += get<Field3D>(species["energy_source"]); // For diagnostic output
  }
  if (species.isSet("energy_source")) {
    Se += V * get<Field3D>(species["momentum_source"]);
  }
  ddt(E) += Se;

  // Scale time derivatives
  if (state.isSet("scale_timederivs")) {
    ddt(E) *= get<Field3D>(state["scale_timederivs"]);
  }

  if (evolve_log) {
    ddt(logE) = ddt(E) / E;
  }

#if CHECKLEVEL >= 1
  for (auto& i : E.getRegion("RGN_NOBNDRY")) {
    if (!std::isfinite(ddt(E)[i])) {
      throw BoutException("ddt(E{}) non-finite at {}. Se={}\n", name, i, Se[i]);
    }
  }
#endif
}

void EvolveEnergy::outputVars(Options& state) {
  AUTO_TRACE();
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto rho_s0 = get<BoutReal>(state["rho_s0"]);

  BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

  if (evolve_log) {
    state[std::string("E") + name].force(E);
  }

  state[std::string("E") + name].setAttributes({{"time_dimension", "t"},
                                                {"units", "J/m^3"},
                                                {"conversion", Pnorm},
                                                {"standard_name", "energy density"},
                                                {"long_name", name + " energy density"},
                                                {"species", name},
                                                {"source", "evolve_energy"}});

  set_with_attrs(state[std::string("P") + name], P,
                 {{"time_dimension", "t"},
                  {"units", "Pa"},
                  {"conversion", Pnorm},
                  {"standard_name", "pressure"},
                  {"long_name", name + " pressure"},
                  {"species", name},
                  {"source", "evolve_energy"}});

  if (diagnose) {
    if (thermal_conduction) {
      set_with_attrs(state[std::string("kappa_par_") + name], kappa_par,
                     {{"time_dimension", "t"},
                      {"units", "W / m / eV"},
                      {"conversion", Pnorm * Omega_ci * SQ(rho_s0)},
                      {"long_name", name + " heat conduction coefficient"},
                      {"species", name},
                      {"source", "evolve_energy"}});
    }
    set_with_attrs(state[std::string("T") + name], T,
                   {{"time_dimension", "t"},
                    {"units", "eV"},
                    {"conversion", Tnorm},
                    {"standard_name", "temperature"},
                    {"long_name", name + " temperature"},
                    {"species", name},
                    {"source", "evolve_energy"}});

    set_with_attrs(
        state[std::string("ddt(E") + name + std::string(")")], ddt(E),
        {{"time_dimension", "t"},
         {"units", "J/m^3/s"},
         {"conversion", Pnorm * Omega_ci},
         {"long_name", std::string("Rate of change of ") + name + " energy density"},
         {"species", name},
         {"source", "evolve_energy"}});

    set_with_attrs(state[std::string("SE") + name], Se,
                   {{"time_dimension", "t"},
                    {"units", "W/m^3"},
                    {"conversion", Pnorm * Omega_ci},
                    {"standard_name", "energy source"},
                    {"long_name", name + " energy source"},
                    {"species", name},
                    {"source", "evolve_energy"}});

    set_with_attrs(state[std::string("E") + name + std::string("_src")], source,
                   {{"time_dimension", "t"},
                    {"units", "W/m^3"},
                    {"conversion", Pnorm * Omega_ci},
                    {"standard_name", "energy source"},
                    {"long_name", name + " energy source"},
                    {"species", name},
                    {"source", "evolve_energy"}});
  }
}

void EvolveEnergy::precon(const Options& state, BoutReal gamma) {
  if (!(enable_precon and thermal_conduction)) {
    return; // Disabled
  }

  static std::unique_ptr<InvertParDiv> inv;
  if (!inv) {
    // Initialise parallel inversion class
    inv = InvertParDiv::create();
    inv->setCoefA(1.0);
  }
  const auto& species = state["species"][name];
  const Field3D N = get<Field3D>(species["density"]);

  // Set the coefficient in Div_par( B * Grad_par )
  Field3D coef = -gamma * kappa_par / floor(N, density_floor);

  if (state.isSet("scale_timederivs")) {
    coef *= get<Field3D>(state["scale_timederivs"]);
  }

  inv->setCoefB(coef);
  Field3D dT = ddt(P);
  dT.applyBoundary("neumann");
  ddt(P) = inv->solve(dT);
}
