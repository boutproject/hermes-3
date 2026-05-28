
#include "../include/sheath_boundary_penalty.hxx"

#include <bout/constants.hxx>
#include <bout/globals.hxx>
#include <bout/mesh.hxx>
using bout::globals::mesh;

SheathBoundaryPenalty::SheathBoundaryPenalty(std::string name, Options& alloptions,
                                             Solver*)
  : Component({
      readIfSet("fields:phi"),
      readIfSet("species:{all_species}:charge"),
      readIfSet("species:{all_species}:AA"),
      // Electrons
      readIfSet("species:e:density", Regions::Interior),
      readIfSet("species:e:temperature", Regions::Interior),
      readIfSet("species:e:pressure", Regions::Interior),
      readIfSet("species:e:velocity", Regions::Interior),
      readIfSet("species:e:momentum", Regions::Interior),
      readWrite("species:e:density_source"),
      readWrite("species:e:momentum_source"),
      readWrite("species:e:energy_source"),
      // Ions
      readIfSet("species:{ions}:adiabatic"),
      readIfSet("species:{ions}:density", Regions::Interior),
      readIfSet("species:{ions}:temperature", Regions::Interior),
      readIfSet("species:{ions}:pressure", Regions::Interior),
      readIfSet("species:{ions}:velocity", Regions::Interior),
      readIfSet("species:{ions}:momentum", Regions::Interior),
      readWrite("species:{ions}:density_source"),
      readWrite("species:{ions}:momentum_source"),
      readWrite("species:{ions}:energy_source"),
      writeFinal("species:{ions}:density_penalty"),
      writeFinal("species:{ions}:momentum_penalty"),
      writeFinal("species:{ions}:energy_penalty"),
    }) {
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

  surface_terms = options["surface_terms"]
    .doc("Include surface terms?")
    .withDefault<bool>(false);

  std::string mask_name = options["mask_name"]
                              .doc("Name of the mesh variable containing penalty mask")
                              .withDefault<std::string>("penalty_mask");

  if (mesh->get(penalty_mask, mask_name) != 0) {
    throw BoutException("Could not read penalty mask variable '{}'", mask_name);
  }
  penalty_mask.applyBoundary("neumann");

  // Find every cell that has penalty_mask > 0
  // so we can efficiently iterate over them later
  {
    Region<Ind3D>::RegionIndices indices;
    BOUT_FOR_SERIAL(i, penalty_mask.getRegion("RGN_NOBNDRY")) {
      if (penalty_mask[i] > 1e-5) {
	// Add this cell to the iteration
	indices.push_back(i);
      }
    }
    penalty_region = Region<Ind3D>(indices);
  }

  if (surface_terms) {
    // Calculate surface terms using field-aligned coordinates
    penalty_mask_fa = toFieldAligned(penalty_mask);
    Region<Ind3D>::RegionIndices indices;
    BOUT_FOR_SERIAL(i, penalty_mask.getRegion("RGN_NOBNDRY")) {
      if (penalty_mask[i] > 1e-5) {
	// Add this cell to the iteration
	indices.push_back(i);
      }
    }
    penalty_region_fa = Region<Ind3D>(indices);

    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
	auto i = indexAt(penalty_mask_fa, r.ind, mesh->ystart, jz);
	penalty_mask_fa[i.ym()] = penalty_mask_fa[i];
      }
    }
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
	auto i = indexAt(penalty_mask_fa, r.ind, mesh->yend, jz);
	penalty_mask_fa[i.yp()] = penalty_mask_fa[i];
      }
    }
  }
}

void SheathBoundaryPenalty::transform_impl(GuardedOptions& state) {

  GuardedOptions allspecies = state["species"];

  // Electrons
  auto electrons = allspecies["e"];
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
                     ? GET_NOBOUNDARY(Field3D, electrons["pressure"])
                     : Ne * Te;

    Field3D Ve = electrons.isSet("velocity")
                     ? GET_NOBOUNDARY(Field3D, electrons["velocity"])
                     : zeroFrom(Ne);

    Field3D NVe = electrons.isSet("momentum")
                      ? GET_NOBOUNDARY(Field3D, electrons["momentum"])
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
    }

    if (surface_terms and has_phi) {
      // Surface penalty terms, to impose sheath current
      // The gradient of the mask gives the direction into the wall

      auto Ne_fa = toFieldAligned(Ne);
      auto Te_fa = toFieldAligned(Te);
      auto Ve_fa = toFieldAligned(Ve);
      auto phi_fa = toFieldAligned(phi);

      auto momentum_source_fa = toFieldAligned(momentum_source);
      BOUT_FOR(i, penalty_region_fa) {
	const auto iyp = i.yp();
	const BoutReal mask = penalty_mask_fa[i];
	const BoutReal dmask_yup = penalty_mask_fa[iyp] - mask;
        if (std::fabs(dmask_yup) > 1e-5) {
	  const BoutReal nfloor = BOUTMAX(Ne_fa[i], 1e-5);
          const BoutReal tesheath = 0.5 * (Te_fa[i] + Te_fa[iyp]);
          const BoutReal vesheath = 0.5 * (Ve_fa[i] + Ve_fa[iyp]);
          const BoutReal phisheath = 0.5 * (phi_fa[i] + phi_fa[iyp]);

          const BoutReal Cse =
              sqrt(tesheath / (TWOPI * Me)) * exp(-phisheath / BOUTMAX(tesheath, 1e-5));

          momentum_source_fa[i] += mask * std::fabs(dmask_yup) * Me * nfloor
	    * (SIGN(dmask_yup) * Cse - vesheath) / penalty_timescale;
        }

	const auto iym = i.ym();
        const BoutReal dmask_ydown = mask - penalty_mask_fa[iym];
        if (std::fabs(dmask_ydown) > 1e-5) {
	  const BoutReal nfloor = BOUTMAX(Ne_fa[i], 1e-5);
          const BoutReal tesheath = 0.5 * (Te_fa[i] + Te_fa[iym]);
          const BoutReal vesheath = 0.5 * (Ve_fa[i] + Ve_fa[iym]);
          const BoutReal phisheath = 0.5 * (phi_fa[i] + phi_fa[iym]);

          const BoutReal Cse =
              sqrt(tesheath / (TWOPI * Me)) * exp(-phisheath / BOUTMAX(tesheath, 1e-5));

          momentum_source_fa[i] += mask * std::fabs(dmask_ydown) * Me * nfloor
	    * (SIGN(dmask_ydown) * Cse - vesheath)
	    / penalty_timescale;
        }
      }
      momentum_source = fromFieldAligned(momentum_source_fa);
    }

    set(electrons["density_source"], density_source);
    set(electrons["momentum_source"], momentum_source);
    set(electrons["energy_source"], energy_source);
  }

  for (auto& kv : allspecies.getChildren()) {
    if (kv.first == "e") {
      continue; // Skip electrons
    }
    GuardedOptions species = allspecies[kv.first]; // Note: Need non-const

    // Characteristics of this species
    if (!IS_SET(species["charge"])) {
      // Skip neutrals
      continue;
    }
    const BoutReal Zi = GET_VALUE(BoutReal, species["charge"]);
    if (fabs(Zi) < 1e-5) {
      // Skip neutrals
      continue;
    }

    const BoutReal Mi = GET_VALUE(BoutReal, species["AA"]);

    Field3D Ni = GET_NOBOUNDARY(Field3D, species["density"]);
    Field3D Ti = GET_NOBOUNDARY(Field3D, species["temperature"]);

    Field3D Pi =
        species.isSet("pressure") ? GET_NOBOUNDARY(Field3D, species["pressure"]) : Ni * Ti;

    Field3D Vi = species.isSet("velocity") ? GET_NOBOUNDARY(Field3D, species["velocity"])
                                           : zeroFrom(Ni);

    Field3D NVi = species.isSet("momentum") ? GET_NOBOUNDARY(Field3D, species["momentum"])
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
    }

    if (surface_terms) {
      // Surface penalty terms.

      auto Ni_fa = toFieldAligned(Ni);
      auto Ti_fa = toFieldAligned(Ti);
      auto Te_fa = toFieldAligned(Te);
      auto Vi_fa = toFieldAligned(Vi);

      auto momentum_penalty_fa = toFieldAligned(momentum_penalty);
      BOUT_FOR(i, penalty_region_fa) {
        const auto iyp = i.yp();
	const BoutReal mask = penalty_mask_fa[i];
	// The gradient of the mask gives the direction
	const BoutReal dmask_yup = penalty_mask_fa[iyp] - mask;
	if (std::fabs(dmask_yup) > 1e-5) {
	  const BoutReal nfloor = BOUTMAX(Ni_fa[i], 1e-5);
	  const BoutReal tisheath = 0.5 * (Ti_fa[i] + Ti_fa[iyp]);
	  const BoutReal tesheath = 0.5 * (Te_fa[i] + Te_fa[iyp]);
	  const BoutReal visheath = 0.5 * (Vi_fa[i] + Vi_fa[iyp]);

	  const BoutReal Cs = sqrt((tesheath + tisheath) / Mi);

	  momentum_penalty_fa[i] += mask * std::fabs(dmask_yup) * Mi * nfloor
	    * (SIGN(dmask_yup) * Cs - visheath) / penalty_timescale;
	}

	const auto iym = i.ym();
	const BoutReal dmask_ydown = mask - penalty_mask_fa[iym];
	if (std::fabs(dmask_ydown) > 1e-5) {
	  const BoutReal nfloor = BOUTMAX(Ni_fa[i], 1e-5);
	  const BoutReal tisheath = 0.5 * (Ti_fa[i] + Ti_fa[iym]);
	  const BoutReal tesheath = 0.5 * (Te_fa[i] + Te_fa[iym]);
	  const BoutReal visheath = 0.5 * (Vi_fa[i] + Vi_fa[iym]);

	  const BoutReal Cs = sqrt((tesheath + tisheath) / Mi);

	  momentum_penalty_fa[i] += mask * std::fabs(dmask_ydown) * Mi * nfloor
	    * (SIGN(dmask_ydown) * Cs - visheath) / penalty_timescale;
	}
      }
      momentum_penalty = fromFieldAligned(momentum_penalty_fa);
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
