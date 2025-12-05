#include <cmath>
#include <string>
#include <vector>

#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/constants.hxx>
#include <bout/difops.hxx>
#include <bout/field2d.hxx>
#include <bout/field3d.hxx>
#include <bout/field_factory.hxx>
#include <bout/fv_ops.hxx>
#include <bout/globals.hxx>
#include <bout/options.hxx>
#include <bout/output.hxx>
#include <bout/solver.hxx>
#include <bout/sys/range.hxx>
#include <bout/utils.hxx>
#include <fmt/format.h>

#include "../include/braginskii_conduction.hxx"
#include "../include/component.hxx"
#include "../include/div_ops.hxx"
#include "../include/hermes_utils.hxx"

using bout::globals::mesh;

BraginskiiConduction::BraginskiiConduction(const std::string& name, Options& alloptions,
                                           Solver* solver) {
  AUTO_TRACE();

  // Get settings for each species
  for (auto& kv : alloptions.getChildren()) {
    auto& options = alloptions[kv.first];
    // Check if the component is a species which undergoes energy/pressure evolution
    if (options.isValue() || !options["type"].isValue()
        || (options["type"].as<std::string>().find("evolve_pressure") == std::string::npos
            && options["type"].as<std::string>().find("evolve_energy")
                   == std::string::npos)) {
      continue;
    }
    if (!options["thermal_conduction"]
             .doc("Include parallel heat conduction?")
             .withDefault<bool>(true)) {
      continue;
    }
    const std::string name = kv.first;

    all_diagnose[name] = options["diagnose"]
                             .doc("Save additional output diagnostics")
                             .withDefault<bool>(false);

    BoutReal default_kappa; // default conductivity, changes depending on species
    switch (identifySpeciesType(name)) {
    case SpeciesType::ion:
      default_kappa = 3.9;
      break;
    case SpeciesType::electron:
      // Hermes-3 electron collision time is in Fitzpatrick form (3.187 in
      // https://farside.ph.utexas.edu/teaching/plasma/Plasma/node41.html) This means that
      // the Braginskii prefactor of 3.16 needs to be divided by sqrt(2) to be consistent.
      default_kappa = 3.16 / sqrt(2);
      break;
    case SpeciesType::neutral:
      default_kappa = 2.5;
      break;
    default:
      throw BoutException("Unhandled species type in default_kappa switch");
    }

    all_kappa_coefficient[name] =
        options["kappa_coefficient"]
            .doc("Numerical coefficient in parallel heat conduction. Default is "
                 "3.16/sqrt(2) for electrons, 2.5 for neutrals and 3.9 otherwise")
            .withDefault(default_kappa);

    all_kappa_limit_alpha[name] =
        options["kappa_limit_alpha"]
            .doc("Flux limiter factor. < 0 means no limit. Typical is 0.2 "
                 "for electrons, 1 for ions.")
            .withDefault(-1.0);
    all_conduction_collisions_mode[name] =
        options["conduction_collisions_mode"]
            .doc("Can be multispecies: all collisions, or "
                 "braginskii: self collisions and ie")
            .withDefault<std::string>("multispecies");
  }
}

void BraginskiiConduction::transform(Options& state) {
  AUTO_TRACE();

  for (auto& kv : state["species"].getChildren()) {
    const std::string& name = kv.first;
    // Skip any species for which energy is not being evolved
    if (all_conduction_collisions_mode.count(name) == 0) {
      continue;
    }
    /// Get the section containing this species
    auto& species = state["species"][name];
    std::string const conduction_collisions_mode = all_conduction_collisions_mode[name];

    // Braginskii mode: plasma - self collisions and ei, neutrals - CX, IZ
    if (all_collision_names.count(name) == 0) { /// Calculate only once - at the beginning
      std::vector<std::string> collision_names;
      const auto species_type = identifySpeciesType(species.name());

      if (conduction_collisions_mode == "braginskii") {
        for (const auto& collision : species["collision_frequencies"].getChildren()) {

          std::string collision_name = collision.second.name();

          if (species_type == SpeciesType::neutral) {
            throw BoutException("\tBraginskii conduction collisions mode not available "
                                "for neutrals, choose multispecies or afn");
          } else if (species_type == SpeciesType::electron) {
            if ( /// Electron-electron collisions
                (collisionSpeciesMatch(collision_name, species.name(), "e", "coll",
                                       "exact"))) {
              collision_names.push_back(collision_name);
            }

          } else if (species_type == SpeciesType::ion) {
            if ( /// Self-collisions
                (collisionSpeciesMatch(collision_name, species.name(), species.name(),
                                       "coll", "exact"))) {
              collision_names.push_back(collision_name);
            }
          }
        }

        // Multispecies mode: all collisions and CX are included
      } else if (conduction_collisions_mode == "multispecies") {
        for (const auto& collision : species["collision_frequencies"].getChildren()) {

          std::string const collision_name = collision.second.name();

          if ( /// Charge exchange
              (collisionSpeciesMatch(collision_name, species.name(), "", "cx", "partial"))
              or
              /// Any collision (en, in, ee, ii, nn)
              (collisionSpeciesMatch(collision_name, species.name(), "", "coll",
                                     "partial"))) {
            collision_names.push_back(collision_name);
          }
        }

      } else if (conduction_collisions_mode == "afn") {
        for (const auto& collision : species["collision_frequencies"].getChildren()) {

          std::string const collision_name = collision.second.name();

          if (species_type != SpeciesType::neutral) {
            throw BoutException("\tAFN conduction collisions mode not available for ions "
                                "or electrons, choose braginskii or multispecies");
          }
          if ( /// Charge exchange
              (collisionSpeciesMatch(collision_name, species.name(), "+", "cx",
                                     "partial"))
              or
              /// Ionisation
              (collisionSpeciesMatch(collision_name, species.name(), "+", "iz",
                                     "partial"))) {
            collision_names.push_back(collision_name);
          }
        }

      } else {
        throw BoutException("\tConduction_collisions_mode incorrect", species.name());
      }

      if (collision_names.empty()) {
        throw BoutException("\tNo collisions found for {:s} in evolve_pressure for "
                            "selected collisions mode",
                            species.name());
      }

      /// Write chosen collisions to log file
      output_info.write("\t{:s} conduction collisionality mode: '{:s}' using ",
                        species.name(), conduction_collisions_mode);
      for (const auto& collision : collision_names) {
        output_info.write("{:s} ", collision);
      }

      output_info.write("\n");
      all_collision_names[name] = collision_names;
    }

    const auto& collision_names = all_collision_names[name];
    Field3D& nu = all_nu[name];
    Field3D& kappa_par = all_kappa_par[name];
    Field3D& flow_ylow_conduction = all_flow_ylow_conduction[name];
    const BoutReal kappa_coefficient = all_kappa_coefficient[name];
    const BoutReal kappa_limit_alpha = all_kappa_limit_alpha[name];

    /// Collect the collisionalities based on list of names
    nu = 0;
    for (const auto& collision_name : collision_names) {
      nu += GET_VALUE(Field3D, species["collision_frequencies"][collision_name]);
    }

    // Get updated pressure and temperature with boundary conditions
    // Note: Retain pressures which fall below zero
    // FIXME: Can I be sure this is the same as the value used in EvolvePressure::finally?
    // FIXME: We end up applying these operations twice: here and in
    // EvolvePressure::finally

    Field3D P = species["pressure"];
    P.clearParallelSlices();
    P.setBoundaryTo(get<Field3D>(species["pressure"]));
    Field3D const Pfloor = floor(P, 0.0); // Restricted to never go below zero
    Field3D const T = get<Field3D>(species["temperature"]);
    Field3D const N = get<Field3D>(species["density"]);

    // Calculate ion collision times
    const Field3D tau = 1. / softFloor(nu, 1e-10);
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
      Field3D const q_SH = kappa_par * Grad_par(T);
      // Free-streaming flux
      Field3D const q_fl = kappa_limit_alpha * N * T * sqrt(T / AA);

      // This results in a harmonic average of the heat fluxes
      kappa_par /= (1. + abs(q_SH / softFloor(q_fl, 1e-10)));

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
    set(species["kappa_par"], kappa_par);
    add(species["energy_source"],
        Div_par_K_Grad_par_mod(kappa_par, T, flow_ylow_conduction, false));
    add(species["energy_flow_ylow"], flow_ylow_conduction);

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
      add(species["energy_source"],
          Div_par(kappa_par * db_dot_T)
              - Div_n_g_bxGrad_f_B_XZ(kappa_par, db_dot_T + b0_dot_T, Apar_flutter));
    }
  }
}

void BraginskiiConduction::outputVars(Options& state) {
  AUTO_TRACE();
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto rho_s0 = get<BoutReal>(state["rho_s0"]);

  BoutReal const Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

  for (const auto& [name, diagnose] : all_diagnose) {
    if (!diagnose) {
      continue;
    }
    const Field3D& kappa_par = all_kappa_par[name];
    const Field3D& nu = all_nu[name];
    const Field3D& flow_ylow_conduction = all_flow_ylow_conduction[name];
    set_with_attrs(state[std::string("kappa_par_") + name], kappa_par,
                   {{"time_dimension", "t"},
                    {"units", "W / m / eV"},
                    {"conversion", (Pnorm * Omega_ci * SQ(rho_s0)) / Tnorm},
                    {"long_name", name + " heat conduction coefficient"},
                    {"species", name},
                    {"source", "evolve_pressure"}});

    set_with_attrs(state[std::string("K") + name + std::string("_cond")], nu,
                   {{"time_dimension", "t"},
                    {"units", "s^-1"},
                    {"conversion", Omega_ci},
                    {"long_name", "collision frequency for conduction"},
                    {"species", name},
                    {"source", "evolve_pressure"}});

    set_with_attrs(state[fmt::format("ef{}_cond_ylow", name)], flow_ylow_conduction,
                   {{"time_dimension", "t"},
                    {"units", "W"},
                    {"conversion", rho_s0 * SQ(rho_s0) * Pnorm * Omega_ci},
                    {"standard_name", "power"},
                    {"long_name", name + " conduction through Y cell face."},
                    {"species", name},
                    {"source", "evolve_pressure"}});
  }
}
