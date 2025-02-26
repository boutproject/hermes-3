
#include "../include/sheath_boundary_penalty.hxx"

#include <bout/constants.hxx>
#include <bout/globals.hxx>
#include <bout/mesh.hxx>
using bout::globals::mesh;

SheathBoundaryPenalty::SheathBoundaryPenalty(std::string name, Options& alloptions,
                                             Solver*) {
  AUTO_TRACE();

  Options& options = alloptions[name];

  diagnose = options["diagnose"]
                .doc("Save penalty term diagnostics?")
                .withDefault<bool>(false);

  gamma_e = options["gamma_e"]
                .doc("Electron sheath heat transmission coefficient")
                .withDefault(3.5);

  gamma_i =
      options["gamma_i"].doc("Ion sheath heat transmission coefficient").withDefault(3.5);

  penalty_timescale = options["penalty_timescale"]
                          .doc("Timescale of penalisation [seconds]")
                          .withDefault(1e-6)
                      / alloptions["units"]["seconds"].as<BoutReal>();

  std::string mask_name = options["mask_name"]
                              .doc("Name of the mesh variable containing penalty mask")
                              .withDefault<std::string>("penalty_mask");

  if (mesh->get(penalty_mask, mask_name) != 0) {
    throw BoutException("Could not read penalty mask variable '{}'", mask_name);
  }
  penalty_mask.applyBoundary("neumann");

  // Find every cell that has penalty_mask > 0
  // so we can efficiently iterate over them later
  Region<Ind3D>::RegionIndices indices;
  BOUT_FOR_SERIAL(i, penalty_mask.getRegion("RGN_NOBNDRY")) {
    if (penalty_mask[i] > 1e-5) {
      // Add this cell to the iteration
      indices.push_back(i);
    }
  }
  penalty_region = Region<Ind3D>(indices);
}

void SheathBoundaryPenalty::transform(Options& state) {
  AUTO_TRACE();

  Options& allspecies = state["species"];

  // Electrons
  auto& electrons = allspecies["e"];
  auto Ne = getNoBoundary<Field3D>(electrons["density"]);
  auto Te = getNoBoundary<Field3D>(electrons["temperature"]);
  const BoutReal Me = get<BoutReal>(electrons["AA"]);

  Field3D phi;
  bool has_phi = IS_SET_NOBOUNDARY(state["fields"]["phi"]);
  if (has_phi) {
    phi = getNoBoundary<Field3D>(state["fields"]["phi"]);
  }

  {
    Field3D Pe = electrons.isSet("pressure")
                     ? getNoBoundary<Field3D>(electrons["pressure"])
                     : Ne * Te;

    Field3D Ve = electrons.isSet("velocity")
                     ? getNoBoundary<Field3D>(electrons["velocity"])
                     : zeroFrom(Ne);

    Field3D NVe = electrons.isSet("momentum")
                      ? getNoBoundary<Field3D>(electrons["momentum"])
                      : Me * Ne * Ve;

    Field3D density_source = electrons.isSet("density_source")
                                 ? getNonFinal<Field3D>(electrons["density_source"])
                                 : zeroFrom(Ne);

    Field3D momentum_source = electrons.isSet("momentum_source")
                                  ? getNonFinal<Field3D>(electrons["momentum_source"])
                                  : zeroFrom(Ne);

    Field3D energy_source = electrons.isSet("energy_source")
                                ? getNonFinal<Field3D>(electrons["energy_source"])
                                : zeroFrom(Ne);

    BOUT_FOR(i, penalty_region) {
      BoutReal mask = penalty_mask[i]; // 1 in boundary

      BoutReal nfloor = BOUTMAX(Ne[i], 1e-5);
      density_source[i] = (1 - mask) * density_source[i]
                          - mask * BOUTMAX(Ne[i] - 1e-5, 0.0) / penalty_timescale;
      momentum_source[i] = (1 - mask) * momentum_source[i]
                           - mask * Me * nfloor * Ve[i] / penalty_timescale;
      energy_source[i] = (1 - mask) * energy_source[i]
                         - mask * gamma_e * nfloor * Te[i] / penalty_timescale;

      if (has_phi) {
        // Surface penalty terms, to impose sheath current
        // The gradient of the mask gives the direction into the wall
        BoutReal dmask_yup = penalty_mask[i.yp()] - mask;
        if (std::fabs(dmask_yup) > 1e-5) {
          const auto iyp = i.yp();
          const BoutReal tesheath = 0.5 * (Te[i] + Te[iyp]);
          const BoutReal vesheath = 0.5 * (Ve[i] + Ve[iyp]);
          const BoutReal phisheath = 0.5 * (phi[i] + phi[iyp]);

          const BoutReal Cse =
              sqrt(tesheath / (TWOPI * Me)) * exp(-phisheath / BOUTMAX(tesheath, 1e-5));

          momentum_source[i] += mask * std::fabs(dmask_yup) * Me * nfloor
                                * (SIGN(dmask_yup) * Cse - vesheath) / penalty_timescale;
        }

        BoutReal dmask_ydown = mask - penalty_mask[i.ym()];
        if (std::fabs(dmask_ydown) > 1e-5) {
          const auto iym = i.ym();
          const BoutReal tesheath = 0.5 * (Te[i] + Te[iym]);
          const BoutReal vesheath = 0.5 * (Ve[i] + Ve[iym]);
          const BoutReal phisheath = 0.5 * (phi[i] + phi[iym]);

          const BoutReal Cse =
              sqrt(tesheath / (TWOPI * Me)) * exp(-phisheath / BOUTMAX(tesheath, 1e-5));

          momentum_source[i] += mask * std::fabs(dmask_ydown) * Me * nfloor
                                * (SIGN(dmask_ydown) * Cse - vesheath)
                                / penalty_timescale;
        }
      }
    }

    set(electrons["density_source"], density_source);
    set(electrons["momentum_source"], momentum_source);
    set(electrons["energy_source"], energy_source);
  }

  for (auto& kv : allspecies.getChildren()) {
    if (kv.first == "e") {
      continue; // Skip electrons
    }
    Options& species = allspecies[kv.first]; // Note: Need non-const

    // Characteristics of this species
    const BoutReal Mi = GET_VALUE(BoutReal, species["AA"]);
    const BoutReal Zi = GET_VALUE(BoutReal, species["charge"]);

    if (fabs(Zi) < 1e-5) {
      // Skip neutrals
      continue;
    }

    // Boundary conditions should already have been applied
    Field3D Ni = GET_VALUE(Field3D, species["density"]);
    Field3D Ti = GET_VALUE(Field3D, species["temperature"]);

    Field3D Pi =
        species.isSet("pressure") ? GET_VALUE(Field3D, species["pressure"]) : Ni * Ti;

    Field3D Vi = species.isSet("velocity") ? GET_VALUE(Field3D, species["velocity"])
                                           : zeroFrom(Ni);

    Field3D NVi = species.isSet("momentum") ? GET_VALUE(Field3D, species["momentum"])
                                            : Mi * Ni * Vi;

    // Get the particle source to modify
    Field3D density_source = species.isSet("density_source")
                                 ? getNonFinal<Field3D>(species["density_source"])
                                 : zeroFrom(Ni);

    Field3D momentum_source = species.isSet("momentum_source")
                                  ? getNonFinal<Field3D>(species["momentum_source"])
                                  : zeroFrom(Ni);

    Field3D energy_source = species.isSet("energy_source")
                                ? getNonFinal<Field3D>(species["energy_source"])
                                : zeroFrom(Ni);

    // Save the sources as diagnostics and for recycling
    Field3D density_penalty {zeroFrom(Ni)};
    Field3D momentum_penalty {zeroFrom(Ni)};
    Field3D energy_penalty {zeroFrom(Ni)};

    BOUT_FOR(i, penalty_region) {
      BoutReal mask = penalty_mask[i]; // 1 in boundary, 0 in the plasma domain

      // Volumetric penalty terms
      BoutReal nfloor = BOUTMAX(Ni[i], 1e-5);
      density_penalty[i] = - mask * density_source[i]
                           - mask * BOUTMAX(Ni[i] - 1e-5, 0.0) / penalty_timescale;

      momentum_penalty[i] = - mask * momentum_source[i]
                            - mask * Mi * nfloor * Vi[i] / penalty_timescale;

      energy_penalty[i] = - mask * energy_source[i]
                          - mask * gamma_i * nfloor * Ti[i] / penalty_timescale;

      // Surface penalty terms.
      // The gradient of the mask gives the direction
      BoutReal dmask_yup = penalty_mask[i.yp()] - mask;
      if (std::fabs(dmask_yup) > 1e-5) {
        const auto iyp = i.yp();
        const BoutReal tisheath = 0.5 * (Ti[i] + Ti[iyp]);
        const BoutReal tesheath = 0.5 * (Te[i] + Te[iyp]);
        const BoutReal visheath = 0.5 * (Vi[i] + Vi[iyp]);

        const BoutReal Cs = sqrt((tesheath + tisheath) / Mi);

        momentum_penalty[i] += mask * std::fabs(dmask_yup) * Mi * nfloor
                              * (SIGN(dmask_yup) * Cs - visheath) / penalty_timescale;
      }

      BoutReal dmask_ydown = mask - penalty_mask[i.ym()];
      if (std::fabs(dmask_ydown) > 1e-5) {
        const auto iym = i.ym();
        const BoutReal tisheath = 0.5 * (Ti[i] + Ti[iym]);
        const BoutReal tesheath = 0.5 * (Te[i] + Te[iym]);
        const BoutReal visheath = 0.5 * (Vi[i] + Vi[iym]);

        const BoutReal Cs = sqrt((tesheath + tisheath) / Mi);

        momentum_penalty[i] += mask * std::fabs(dmask_ydown) * Mi * nfloor
                              * (SIGN(dmask_ydown) * Cs - visheath) / penalty_timescale;
      }
    }

    add(species["density_source"], density_penalty);
    add(species["momentum_source"], momentum_penalty);
    add(species["energy_source"], energy_penalty);

    set(species["density_penalty"], density_penalty);
    set(species["momentum_penalty"], momentum_penalty);
    set(species["energy_penalty"], energy_penalty);

    if (diagnose) {
      // Store penalty term diagnostics to used in outputVars
      auto& diagnostic_species = diagnostics[kv.first];
      set(diagnostic_species["density_penalty"], density_penalty);
      set(diagnostic_species["momentum_penalty"], momentum_penalty);
      set(diagnostic_species["energy_penalty"], energy_penalty);
    }
  }
}

void SheathBoundaryPenalty::outputVars(Options& state) {
  AUTO_TRACE();

  set_with_attrs(state["penalty_mask"], penalty_mask,
                 {{"units", ""},
                  {"long_name", "Penalty mask"},
                  {"standard_name", "Penalty mask"},
                  {"source", "sheath_boundary_penalty"}});

  if (diagnose) {
    // Normalisations
    auto Nnorm = get<BoutReal>(state["Nnorm"]);
    auto Tnorm = get<BoutReal>(state["Tnorm"]);
    BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation
    auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
    auto Cs0 = get<BoutReal>(state["Cs0"]);

    for (const auto& kv : diagnostics.getChildren()) {
      // Save diagnostics for this species
      set_with_attrs(state[std::string("S") + kv.first + std::string("_penalty")],
                     getNonFinal<Field3D>(kv.second["density_penalty"]),
                     {{"time_dimension", "t"},
                      {"units", "m^-3 s^-1"},
                      {"conversion", Nnorm * Omega_ci},
                      {"standard_name", "particle source"},
                      {"long_name", kv.first + " particle source penalty term"},
                      {"species", kv.first},
                      {"source", "sheath_boundary_penalty"}});

      set_with_attrs(state[std::string("F") + kv.first + std::string("_penalty")],
                     getNonFinal<Field3D>(kv.second["momentum_penalty"]),
                     {{"time_dimension", "t"},
                      {"units", "kg m^-2 s^-2"},
                      {"conversion", SI::Mp * Nnorm * Cs0 * Omega_ci},
                      {"standard_name", "momentum source"},
                      {"long_name", kv.first + " momentum source penalty term"},
                      {"species", kv.first},
                      {"source", "sheath_boundary_penalty"}});

      set_with_attrs(state[std::string("R") + kv.first + std::string("_penalty")],
                     getNonFinal<Field3D>(kv.second["energy_penalty"]),
                     {{"time_dimension", "t"},
                      {"units", "W m^-3"},
                      {"conversion", Pnorm * Omega_ci},
                      {"standard_name", "energy source"},
                      {"long_name", kv.first + " energy source penalty term"},
                      {"species", kv.first},
                      {"source", "sheath_boundary_penalty"}});
    }
  }
}
