
#include <bout/constants.hxx>
#include <bout/fv_ops.hxx>
#include <derivs.hxx>
#include <difops.hxx>
#include <initialprofiles.hxx>

#include "../include/div_ops.hxx"
#include "../include/evolve_pressure.hxx"
#include "../include/hermes_utils.hxx"

using bout::globals::mesh;

EvolvePressure::EvolvePressure(std::string name, Options& alloptions, Solver* solver)
    : name(name) {
  AUTO_TRACE();

  auto& options = alloptions[name];

  evolve_log = options["evolve_log"].doc("Evolve the logarithmof pressure?").withDefault<bool>(false);

  density_floor = options["density_floor"].doc("Minimum density floor").withDefault(1e-5);
  if (evolve_log) {
    // Evolve logarithm of density
    solver->add(logP, std::string("logP") + name);
    // Save the density to the restart file
    // so the simulation can be restarted evolving density
    get_restart_datafile()->addOnce(P, std::string("P") + name);
    // Save density to output files
    bout::globals::dump.addRepeat(P, std::string("P") + name);

    if (!alloptions["hermes"]["restarting"]) {
      // Set logN from N input options
      initial_profile(std::string("P") + name, P);
      logP = log(P);
    } else {
      // Ignore these settings
      Options::root()[std::string("P") + name].setConditionallyUsed();
    }
  } else {
    // Evolve the density in time
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

  if (options["diagnose"]
          .doc("Save additional output diagnostics")
          .withDefault<bool>(false)) {
    if (thermal_conduction) {
      bout::globals::dump.addRepeat(kappa_par, std::string("kappa_par_") + name);
    }
    bout::globals::dump.addRepeat(T, std::string("T") + name);

    bout::globals::dump.addRepeat(ddt(P), std::string("ddt(P") + name + std::string(")"));
    bout::globals::dump.addRepeat(Sp, std::string("SP") + name);
    Sp = 0.0;
  }

  const Options& units = alloptions["units"];
  const BoutReal Nnorm = units["inv_meters_cubed"];
  const BoutReal Tnorm = units["eV"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();

  source = alloptions[std::string("P") + name]["source"]
               .doc(std::string("Source term in ddt(P") + name
                    + std::string("). Units [N/m^2/s]"))
               .withDefault(Field3D(0.0))
           / (SI::qe * Nnorm * Tnorm * Omega_ci);
}

void EvolvePressure::transform(Options& state) {
  AUTO_TRACE();

  if (evolve_log) {
    // Evolving logP, but most calculations use P
    P = exp(logP);
  } else {
    // Check for low pressures, ensure that Pi >= 0
    P = floor(P, 0.0);
  }

  mesh->communicate(P);

  auto& species = state["species"][name];

  set(species["pressure"], P);

  // Calculate temperature
  // Not using density boundary condition
  N = getNoBoundary<Field3D>(species["density"]);
  T = P / floor(N, density_floor);

  set(species["temperature"], T);
}

void EvolvePressure::finally(const Options& state) {
  AUTO_TRACE();

  /// Get the section containing this species
  const auto& species = state["species"][name];

  // Get updated pressure and temperature with boundary conditions
  P = get<Field3D>(species["pressure"]);
  T = get<Field3D>(species["temperature"]);

  if (state.isSection("fields") and state["fields"].isSet("phi")) {
    // Electrostatic potential set -> include ExB flow

    Field3D phi = get<Field3D>(state["fields"]["phi"]);

    ddt(P) = -Div_n_bxGrad_f_B_XPPM(P, phi, bndry_flux, poloidal_flows, true);
  } else {
    ddt(P) = 0.0;
  }

  if (species.isSet("velocity")) {
    Field3D V = get<Field3D>(species["velocity"]);

    // Typical wave speed used for numerical diffusion
    Field3D sound_speed;
    if (state.isSet("sound_speed")) {
      sound_speed = get<Field3D>(state["sound_speed"]);
    } else {
      Field3D T = get<Field3D>(species["temperature"]);
      BoutReal AA = get<BoutReal>(species["AA"]);
      sound_speed = sqrt(T / AA);
    }

    if (p_div_v) {
      // Use the P * Div(V) form
      ddt(P) -= FV::Div_par(P, V, sound_speed);

      // Work done. This balances energetically a term in the momentum equation
      ddt(P) -= (2. / 3) * P * Div_par(V);

    } else {
      // Use V * Grad(P) form
      // Note: A mixed form has been tried (on 1D neon example)
      //       -(4/3)*FV::Div_par(P,V) + (1/3)*(V * Grad_par(P) - P * Div_par(V))
      //       Caused heating of charged species near sheath like p_div_v
      ddt(P) -= (5. / 3) * FV::Div_par(P, V, sound_speed);

      ddt(P) += (2. / 3) * V * Grad_par(P);
    }
  }

  if (species.isSet("low_n_coeff")) {
    // Low density parallel diffusion
    Field3D low_n_coeff = get<Field3D>(species["low_n_coeff"]);
    ddt(P) += FV::Div_par_K_Grad_par(low_n_coeff * T, N) + FV::Div_par_K_Grad_par(low_n_coeff, P);
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
    kappa_par = kappa_coefficient * P * tau / AA;

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

      Field3D N = get<Field3D>(species["density"]);

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
    ddt(P) += (2. / 3) * FV::Div_par_K_Grad_par(kappa_par, T, false);
  }

  if (hyper_z > 0.) {
    auto* coord = N.getCoordinates();
    ddt(P) -= hyper_z * SQ(SQ(coord->dz)) * D4DZ4(P);
  }

  //////////////////////
  // Other sources

  Sp = source;
  if (species.isSet("energy_source")) {
    Sp += (2. / 3) * get<Field3D>(species["energy_source"]); // For diagnostic output
  }
  ddt(P) += Sp;

#if CHECK > 1
  bout::checkFinite(ddt(P), std::string("ddt P") + name, "RGN_NOBNDRY");
#endif

  if (evolve_log) {
    ddt(logP) = ddt(P) / P;
  }
}
