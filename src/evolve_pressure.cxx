
#include <bout/constants.hxx>
#include <bout/fv_ops.hxx>
#include <bout/field_factory.hxx>
#include <bout/derivs.hxx>
#include <bout/difops.hxx>
#include <bout/output_bout_types.hxx>
#include <bout/initialprofiles.hxx>
#include <bout/invert_pardiv.hxx>

#include "../include/div_ops.hxx"
#include "../include/evolve_pressure.hxx"
#include "../include/hermes_utils.hxx"
#include "../include/hermes_build_config.hxx"

using bout::globals::mesh;

EvolvePressure::EvolvePressure(std::string name, Options& alloptions, Solver* solver)
    : name(name) {
  AUTO_TRACE();

  auto& options = alloptions[name];

  evolve_log = options["evolve_log"].doc("Evolve the logarithm of pressure?").withDefault<bool>(false);

  density_floor = options["density_floor"].doc("Minimum density floor").withDefault(1e-7);

  low_n_diffuse_perp = options["low_n_diffuse_perp"]
                           .doc("Perpendicular diffusion at low density")
                           .withDefault<bool>(false);

  temperature_floor = options["temperature_floor"].doc("Low temperature scale for low_T_diffuse_perp")
    .withDefault<BoutReal>(0.1) / get<BoutReal>(alloptions["units"]["eV"]);

  low_T_diffuse_perp = options["low_T_diffuse_perp"].doc("Add cross-field diffusion at low temperature?")
    .withDefault<bool>(false);

  pressure_floor = density_floor * (1./get<BoutReal>(alloptions["units"]["eV"]));

  low_p_diffuse_perp = options["low_p_diffuse_perp"]
                           .doc("Perpendicular diffusion at low pressure")
                           .withDefault<bool>(false);

  if (evolve_log) {
    // Evolve logarithm of pressure
    solver->add(logP, std::string("logP") + name);
    // Save the pressure to the restart file
    // so the simulation can be restarted evolving pressure
    //get_restart_datafile()->addOnce(P, std::string("P") + name);

    if (!alloptions["hermes"]["restarting"]) {
      // Set logN from N input options
      initial_profile(std::string("P") + name, P);
      logP = log(P);
    } else {
      // Ignore these settings
      Options::root()[std::string("P") + name].setConditionallyUsed();
    }
  } else {
    // Evolve the pressure in time
    solver->add(P, std::string("P") + name);
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
    .doc("Numerical coefficient in parallel heat conduction. Default is 3.16 for electrons, 3.9 otherwise")
    .withDefault((name == "e") ? 3.16 : 3.9);

  kappa_limit_alpha = options["kappa_limit_alpha"]
    .doc("Flux limiter factor. < 0 means no limit. Typical is 0.2 for electrons, 1 for ions.")
    .withDefault(-1.0);

  p_div_v = options["p_div_v"]
                .doc("Use p*Div(v) form? Default, false => v * Grad(p) form")
                .withDefault<bool>(false);

  hyper_z = options["hyper_z"].doc("Hyper-diffusion in Z").withDefault(-1.0);

  hyper_z_T = options["hyper_z_T"]
    .doc("4th-order dissipation of temperature")
    .withDefault<BoutReal>(-1.0);

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

  auto& p_options = alloptions[std::string("P") + name];
  source_normalisation = SI::qe * Nnorm * Tnorm * Omega_ci;   // [Pa/s] or [W/m^3] if converted to energy
  time_normalisation = 1./Omega_ci;   // [s]
  
  // Try to read the pressure source from the mesh
  // Units of Pascals per second
  source = 0.0;
  mesh->get(source, std::string("P") + name + "_src");
  // Allow the user to override the source
  source = p_options["source"]
               .doc(std::string("Source term in ddt(P") + name
                    + std::string("). Units [Pa/s], note P = 2/3 E"))
               .withDefault(source)
           / (source_normalisation);

  source_time_dependent = p_options["source_time_dependent"]
    .doc("Use a time-dependent source?")
    .withDefault<bool>(false);

  // If time dependent, parse the function with respect to time from the input file
  if (source_time_dependent) {
    auto str = p_options["source_prefactor"]
      .doc("Time-dependent function of multiplier on ddt(P" + name + std::string(") source."))
      .as<std::string>();
      source_prefactor_function = FieldFactory::get()->parse(str, &p_options);
  }

  if (p_options["source_only_in_core"]
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

  neumann_boundary_average_z = p_options["neumann_boundary_average_z"]
    .doc("Apply neumann boundary with Z average?")
    .withDefault<bool>(false);

  numerical_viscous_heating = options["numerical_viscous_heating"]
    .doc("Include heating due to numerical viscosity?")
    .withDefault<bool>(false);

  if (numerical_viscous_heating) {
    fix_momentum_boundary_flux = options["fix_momentum_boundary_flux"]
      .doc("Fix Y boundary momentum flux to boundary midpoint value?")
      .withDefault<bool>(false);
  }
}

void EvolvePressure::transform(Options& state) {
  AUTO_TRACE();

  if (evolve_log) {
    // Evolving logP, but most calculations use P
    P = exp(logP);
  }

  mesh->communicate(P);

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

  auto& species = state["species"][name];

  // Calculate temperature
  // Not using density boundary condition
  N = getNoBoundary<Field3D>(species["density"]);

  Field3D Pfloor = floor(P, 0.0);
  T = Pfloor / floor(N, density_floor);
  Pfloor = N * T; // Ensure consistency

  set(species["pressure"], Pfloor);
  set(species["temperature"], T);
}

void EvolvePressure::finally(const Options& state) {
  AUTO_TRACE();

  /// Get the section containing this species
  const auto& species = state["species"][name];

  // Get updated pressure and temperature with boundary conditions
  // Note: Retain pressures which fall below zero
  P.clearParallelSlices();
  P.setBoundaryTo(get<Field3D>(species["pressure"]));
  Field3D Pfloor = floor(P, 0.0); // Restricted to never go below zero

  T = get<Field3D>(species["temperature"]);
  N = get<Field3D>(species["density"]);

  if (species.isSet("charge") and (fabs(get<BoutReal>(species["charge"])) > 1e-5) and
      state.isSection("fields") and state["fields"].isSet("phi")) {
    // Electrostatic potential set and species is charged -> include ExB flow

    Field3D phi = get<Field3D>(state["fields"]["phi"]);

    ddt(P) = -Div_n_bxGrad_f_B_XPPM(P, phi, bndry_flux, poloidal_flows, true);
  } else {
    ddt(P) = 0.0;
  }

  if (species.isSet("velocity")) {
    Field3D V = get<Field3D>(species["velocity"]);

    // Typical wave speed used for numerical diffusion
    Field3D fastest_wave;
    if (state.isSet("fastest_wave")) {
      fastest_wave = get<Field3D>(state["fastest_wave"]);
    } else {
      BoutReal AA = get<BoutReal>(species["AA"]);
      fastest_wave = sqrt(T / AA);
    }

    if (p_div_v) {
      // Use the P * Div(V) form
      ddt(P) -= FV::Div_par_mod<hermes::Limiter>(P, V, fastest_wave, flow_ylow);

      // Work done. This balances energetically a term in the momentum equation
      ddt(P) -= (2. / 3) * Pfloor * Div_par(V);

    } else {
      // Use V * Grad(P) form
      // Note: A mixed form has been tried (on 1D neon example)
      //       -(4/3)*FV::Div_par(P,V) + (1/3)*(V * Grad_par(P) - P * Div_par(V))
      //       Caused heating of charged species near sheath like p_div_v
      ddt(P) -= (5. / 3) * FV::Div_par_mod<hermes::Limiter>(P, V, fastest_wave, flow_ylow);

      ddt(P) += (2. / 3) * V * Grad_par(P);
    }
    flow_ylow *= 5. / 2; // Energy flow

    if (state.isSection("fields") and state["fields"].isSet("Apar_flutter")) {
      // Magnetic flutter term
      const Field3D Apar_flutter = get<Field3D>(state["fields"]["Apar_flutter"]);
      ddt(P) -= (5. / 3) * Div_n_g_bxGrad_f_B_XZ(P, V, -Apar_flutter);
      ddt(P) += (2. / 3) * V * bracket(P, Apar_flutter, BRACKET_ARAKAWA);
    }

    if (numerical_viscous_heating || diagnose) {
      // Viscous heating coming from numerical viscosity
      Field3D Nlim = floor(N, density_floor);
      const BoutReal AA = get<BoutReal>(species["AA"]); // Atomic mass
      Sp_nvh = (2. / 3) * AA * FV::Div_par_fvv_heating(Nlim, V, fastest_wave, flow_ylow_kinetic, fix_momentum_boundary_flux);
      flow_ylow_kinetic *= AA;
      flow_ylow += flow_ylow_kinetic;
      if (numerical_viscous_heating) {
        ddt(P) += Sp_nvh;
      }
    }
  }

  if (species.isSet("low_n_coeff")) {
    // Low density parallel diffusion
    Field3D low_n_coeff = get<Field3D>(species["low_n_coeff"]);
    ddt(P) += FV::Div_par_K_Grad_par(low_n_coeff * T, N) + FV::Div_par_K_Grad_par(low_n_coeff, P);
  }

  if (low_n_diffuse_perp) {
    ddt(P) += Div_Perp_Lap_FV_Index(density_floor / floor(N, 1e-3 * density_floor), P, true);
  }

  if (low_T_diffuse_perp) {
    ddt(P) += 1e-4 * Div_Perp_Lap_FV_Index(floor(temperature_floor / floor(T, 1e-3 * temperature_floor) - 1.0, 0.0),
                                           T, false);
  }

  if (low_p_diffuse_perp) {
    Field3D Plim = floor(P, 1e-3 * pressure_floor);
    ddt(P) += Div_Perp_Lap_FV_Index(pressure_floor / Plim, P, true);
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
    flow_ylow_conduction;
    ddt(P) += (2. / 3) * Div_par_K_Grad_par_mod(kappa_par, T, flow_ylow_conduction, false);
    flow_ylow += flow_ylow_conduction;

    if (state.isSection("fields") and state["fields"].isSet("Apar_flutter")) {
      // Magnetic flutter term. The operator splits into 4 pieces:
      // Div(k b b.Grad(T)) = Div(k b0 b0.Grad(T)) + Div(k d0 db.Grad(T))
      //                    + Div(k db b0.Grad(T)) + Div(k db db.Grad(T))
      // The first term is already calculated above.
      // Here we add the terms containing db
      const Field3D Apar_flutter = get<Field3D>(state["fields"]["Apar_flutter"]);
      Field3D db_dot_T = bracket(T, Apar_flutter, BRACKET_ARAKAWA);
      Field3D b0_dot_T = Grad_par(T);
      mesh->communicate(db_dot_T, b0_dot_T);
      db_dot_T.applyBoundary("neumann");
      b0_dot_T.applyBoundary("neumann");
      ddt(P) += (2. / 3) * (Div_par(kappa_par * db_dot_T) -
                            Div_n_g_bxGrad_f_B_XZ(kappa_par, db_dot_T + b0_dot_T, Apar_flutter));
    }
  }

  if (hyper_z > 0.) {
    ddt(P) -= hyper_z * D4DZ4_Index(P);
  }

  if (hyper_z_T > 0.) {
    ddt(P) -= hyper_z_T * D4DZ4_Index(T);
  }

  //////////////////////
  // Other sources

  if (source_time_dependent) {
    // Evaluate the source_prefactor function at the current time in seconds and scale source with it
    BoutReal time = get<BoutReal>(state["time"]);
    BoutReal source_prefactor = source_prefactor_function ->generate(bout::generator::Context().set("x",0,"y",0,"z",0,"t",time*time_normalisation));
    final_source = source * source_prefactor;
  } else {
    final_source = source;
  }

  Sp = final_source;
  if (species.isSet("energy_source")) {
    Sp += (2. / 3) * get<Field3D>(species["energy_source"]); // For diagnostic output
  }
#if CHECKLEVEL >= 1
  if (species.isSet("pressure_source")) {
    throw BoutException("Components must evolve `energy_source` rather then `pressure_source`");
  }
#endif
  ddt(P) += Sp;

  // Term to force evolved P towards N * T
  // This is active when P < 0 or when N < density_floor
  ddt(P) += N * T - P;

  // Scale time derivatives
  if (state.isSet("scale_timederivs")) {
    ddt(P) *= get<Field3D>(state["scale_timederivs"]);
  }

  if (evolve_log) {
    ddt(logP) = ddt(P) / P;
  }

#if CHECKLEVEL >= 1
  for (auto& i : P.getRegion("RGN_NOBNDRY")) {
    if (!std::isfinite(ddt(P)[i])) {
      throw BoutException("ddt(P{}) non-finite at {}. Sp={}\n", name, i, Sp[i]);
    }
  }
#endif

  if (diagnose) {
    // Save flows of energy if they are set

    if (species.isSet("energy_flow_xlow")) {
      flow_xlow = get<Field3D>(species["energy_flow_xlow"]);
    }
    if (species.isSet("energy_flow_ylow")) {
      flow_ylow += get<Field3D>(species["energy_flow_ylow"]);
    }
  }
}

void EvolvePressure::outputVars(Options& state) {
  AUTO_TRACE();
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto rho_s0 = get<BoutReal>(state["rho_s0"]);

  BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

  if (evolve_log) {
    state[std::string("P") + name].force(P);
  }

  state[std::string("P") + name].setAttributes({{"time_dimension", "t"},
                                                {"units", "Pa"},
                                                {"conversion", Pnorm},
                                                {"standard_name", "pressure"},
                                                {"long_name", name + " pressure"},
                                                {"species", name},
                                                {"source", "evolve_pressure"}});

  if (diagnose) {
    if (thermal_conduction) {
      set_with_attrs(state[std::string("kappa_par_") + name], kappa_par,
                     {{"time_dimension", "t"},
                      {"units", "W / m / eV"},
                      {"conversion", (Pnorm * Omega_ci * SQ(rho_s0) )/ Tnorm},
                      {"long_name", name + " heat conduction coefficient"},
                      {"species", name},
                      {"source", "evolve_pressure"}});
    }
    set_with_attrs(state[std::string("T") + name], T,
                   {{"time_dimension", "t"},
                    {"units", "eV"},
                    {"conversion", Tnorm},
                    {"standard_name", "temperature"},
                    {"long_name", name + " temperature"},
                    {"species", name},
                    {"source", "evolve_pressure"}});

    set_with_attrs(state[std::string("ddt(P") + name + std::string(")")], ddt(P),
                   {{"time_dimension", "t"},
                    {"units", "Pa s^-1"},
                    {"conversion", Pnorm * Omega_ci},
                    {"long_name", std::string("Rate of change of ") + name + " pressure"},
                    {"species", name},
                    {"source", "evolve_pressure"}});

    set_with_attrs(state[std::string("SP") + name], Sp,
                   {{"time_dimension", "t"},
                    {"units", "Pa s^-1"},
                    {"conversion", Pnorm * Omega_ci},
                    {"standard_name", "pressure source"},
                    {"long_name", name + " pressure source"},
                    {"species", name},
                    {"source", "evolve_pressure"}});

    set_with_attrs(state[std::string("P") + name + std::string("_src")], final_source,
                   {{"time_dimension", "t"},
                    {"units", "Pa s^-1"},
                    {"conversion", Pnorm * Omega_ci},
                    {"standard_name", "pressure source"},
                    {"long_name", name + " pressure source"},
                    {"species", name},
                    {"source", "evolve_pressure"}});

    if (flow_xlow.isAllocated()) {
      set_with_attrs(state[fmt::format("ef{}_tot_xlow", name)], flow_xlow,
                   {{"time_dimension", "t"},
                    {"units", "W"},
                    {"conversion", rho_s0 * SQ(rho_s0) * Pnorm * Omega_ci},
                    {"standard_name", "power"},
                    {"long_name", name + " power through X cell face. Note: May be incomplete."},
                    {"species", name},
                    {"source", "evolve_pressure"}});
    }
    if (flow_ylow.isAllocated()) {
      set_with_attrs(state[fmt::format("ef{}_tot_ylow", name)], flow_ylow,
                   {{"time_dimension", "t"},
                    {"units", "W"},
                    {"conversion", rho_s0 * SQ(rho_s0) * Pnorm * Omega_ci},
                    {"standard_name", "power"},
                    {"long_name", name + " power through Y cell face. Note: May be incomplete."},
                    {"species", name},
                    {"source", "evolve_pressure"}});
                    
      set_with_attrs(state[fmt::format("ef{}_cond_ylow", name)], flow_ylow_conduction,
                   {{"time_dimension", "t"},
                    {"units", "W"},
                    {"conversion", rho_s0 * SQ(rho_s0) * Pnorm * Omega_ci},
                    {"standard_name", "power"},
                    {"long_name", name + " conduction through Y cell face. Note: May be incomplete."},
                    {"species", name},
                    {"source", "evolve_pressure"}});

      set_with_attrs(state[fmt::format("ef{}_kin_ylow", name)], flow_ylow_kinetic,
                   {{"time_dimension", "t"},
                    {"units", "W"},
                    {"conversion", rho_s0 * SQ(rho_s0) * Pnorm * Omega_ci},
                    {"standard_name", "power"},
                    {"long_name", name + " kinetic energy flow through Y cell face. Note: May be incomplete."},
                    {"species", name},
                    {"source", "evolve_pressure"}});
    }

    if (numerical_viscous_heating) {
      set_with_attrs(state[std::string("E") + name + std::string("_nvh")], Sp_nvh * 3/.2,
                   {{"time_dimension", "t"},
                    {"units", "W"},
                    {"conversion", Pnorm * Omega_ci},
                    {"standard_name", "energy source"},
                    {"long_name", name + " energy source from numerical viscous heating"},
                    {"species", name},
                    {"source", "evolve_pressure"}});
    }
  }
}

void EvolvePressure::precon(const Options &state, BoutReal gamma) {
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
  Field3D coef = -(2. / 3) * gamma * kappa_par / floor(N, density_floor);

  if (state.isSet("scale_timederivs")) {
    coef *= get<Field3D>(state["scale_timederivs"]);
  }

  inv->setCoefB(coef);
  Field3D dT = ddt(P);
  dT.applyBoundary("neumann");
  ddt(P) = inv->solve(dT);
}
