
#include "../include/sheath_boundary_penalty.hxx"
#include "../include/component.hxx"
#include "../include/guarded_options.hxx"
#include "../include/permissions.hxx"

#include <algorithm>
#include <cmath>
#include <string>

#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/constants.hxx>
#include <bout/field.hxx>
#include <bout/field3d.hxx>
#include <bout/globals.hxx>
#include <bout/mesh.hxx>
#include <bout/region.hxx>
#include <bout/solver.hxx>
#include <bout/sys/range.hxx>
#include <bout/utils.hxx>

using bout::globals::mesh;

Region<Ind3D> SheathBoundaryPenalty::buildPenaltyRegion(const Field3D& mask,
                                                        BoutReal threshold) {
  Region<Ind3D>::RegionIndices indices;
  BOUT_FOR_SERIAL(i, mask.getRegion("RGN_NOBNDRY")) {
    if (mask[i] > threshold) {
      indices.push_back(i);
    }
  }
  return Region<Ind3D>(indices);
}

SheathBoundaryPenalty::PenaltyMaskData
SheathBoundaryPenalty::preparePenaltyMask(Field3D mask, BoutReal threshold) {
  mask.applyBoundary("neumann");
  return {mask, buildPenaltyRegion(mask, threshold)};
}

SheathBoundaryPenalty::PenaltyMaskData
SheathBoundaryPenalty::prepareFieldAlignedPenaltyMask(Field3D mask_fa, Mesh& localmesh,
                                                      BoutReal threshold) {
  auto region = buildPenaltyRegion(mask_fa, threshold);

  for (RangeIterator r = localmesh.iterateBndryLowerY(); !r.isDone(); r++) {
    for (int jz = 0; jz < localmesh.LocalNz; jz++) {
      auto i = indexAt(mask_fa, r.ind, localmesh.ystart, jz);
      mask_fa[i.ym()] = mask_fa[i];
    }
  }
  for (RangeIterator r = localmesh.iterateBndryUpperY(); !r.isDone(); r++) {
    for (int jz = 0; jz < localmesh.LocalNz; jz++) {
      auto i = indexAt(mask_fa, r.ind, localmesh.yend, jz);
      mask_fa[i.yp()] = mask_fa[i];
    }
  }

  return {mask_fa, region};
}

SheathBoundaryPenalty::PenaltySourceData
SheathBoundaryPenalty::calculateVolumetricPenalty(
    const PenaltyMaskData& penalty_data, const Field3D& Ni, const Field3D& Ti,
    const Field3D& Vi, BoutReal Mi, BoutReal gamma_i, BoutReal penalty_timescale,
    const Field3D& density_source, const Field3D& momentum_source,
    const Field3D& energy_source, BoutReal density_floor) {
  Field3D density_penalty{zeroFrom(Ni)};
  Field3D momentum_penalty{zeroFrom(Ni)};
  Field3D energy_penalty{zeroFrom(Ni)};

  BOUT_FOR(i, penalty_data.region) {
    const BoutReal mask = penalty_data.mask[i];
    const BoutReal nfloor = std::max(Ni[i], density_floor);

    density_penalty[i] =
        -mask * density_source[i]
        - mask * std::max(Ni[i] - density_floor, 0.0) / penalty_timescale;
    momentum_penalty[i] =
        -mask * momentum_source[i] - mask * Mi * nfloor * Vi[i] / penalty_timescale;
    energy_penalty[i] =
        -mask * energy_source[i] - mask * gamma_i * nfloor * Ti[i] / penalty_timescale;
  }

  return {density_penalty, momentum_penalty, energy_penalty};
}

Field3D SheathBoundaryPenalty::calculateElectronSurfaceMomentumPenalty(
    const PenaltyMaskData& penalty_fa_data, const Field3D& Ne_fa, const Field3D& Te_fa,
    const Field3D& Ve_fa, const Field3D& phi_fa, BoutReal Me, BoutReal penalty_timescale,
    BoutReal density_floor, BoutReal mask_threshold) {
  Field3D momentum_penalty_fa{zeroFrom(Ne_fa)};
  BOUT_FOR(i, penalty_fa_data.region) {
    const auto iyp = i.yp();
    const BoutReal mask = penalty_fa_data.mask[i];
    const BoutReal dmask_yup = penalty_fa_data.mask[iyp] - mask;
    if (std::fabs(dmask_yup) > mask_threshold) {
      const BoutReal nfloor = std::max(Ne_fa[i], density_floor);
      const BoutReal tesheath = 0.5 * (Te_fa[i] + Te_fa[iyp]);
      const BoutReal vesheath = 0.5 * (Ve_fa[i] + Ve_fa[iyp]);
      const BoutReal phisheath = 0.5 * (phi_fa[i] + phi_fa[iyp]);
      const BoutReal Cse = sqrt(tesheath / (TWOPI * Me))
                           * exp(-phisheath / BOUTMAX(tesheath, density_floor));

      momentum_penalty_fa[i] += mask * std::fabs(dmask_yup) * Me * nfloor
                                * (SIGN(dmask_yup) * Cse - vesheath) / penalty_timescale;
    }

    const auto iym = i.ym();
    const BoutReal dmask_ydown = mask - penalty_fa_data.mask[iym];
    if (std::fabs(dmask_ydown) > mask_threshold) {
      const BoutReal nfloor = BOUTMAX(Ne_fa[i], density_floor);
      const BoutReal tesheath = 0.5 * (Te_fa[i] + Te_fa[iym]);
      const BoutReal vesheath = 0.5 * (Ve_fa[i] + Ve_fa[iym]);
      const BoutReal phisheath = 0.5 * (phi_fa[i] + phi_fa[iym]);
      const BoutReal Cse = sqrt(tesheath / (TWOPI * Me))
                           * exp(-phisheath / BOUTMAX(tesheath, density_floor));

      momentum_penalty_fa[i] += mask * std::fabs(dmask_ydown) * Me * nfloor
                                * (SIGN(dmask_ydown) * Cse - vesheath)
                                / penalty_timescale;
    }
  }

  return momentum_penalty_fa;
}

Field3D SheathBoundaryPenalty::calculateIonSurfaceMomentumPenalty(
    const PenaltyMaskData& penalty_fa_data, const Field3D& Ni_fa, const Field3D& Ti_fa,
    const Field3D& Te_fa, const Field3D& Vi_fa, BoutReal Mi, BoutReal penalty_timescale,
    BoutReal density_floor, BoutReal mask_threshold) {
  Field3D momentum_penalty_fa{zeroFrom(Ni_fa)};
  BOUT_FOR(i, penalty_fa_data.region) {
    const auto iyp = i.yp();
    const BoutReal mask = penalty_fa_data.mask[i];
    const BoutReal dmask_yup = penalty_fa_data.mask[iyp] - mask;
    if (std::fabs(dmask_yup) > mask_threshold) {
      const BoutReal nfloor = BOUTMAX(Ni_fa[i], density_floor);
      const BoutReal tisheath = 0.5 * (Ti_fa[i] + Ti_fa[iyp]);
      const BoutReal tesheath = 0.5 * (Te_fa[i] + Te_fa[iyp]);
      const BoutReal visheath = 0.5 * (Vi_fa[i] + Vi_fa[iyp]);
      const BoutReal Cs = sqrt((tesheath + tisheath) / Mi);

      momentum_penalty_fa[i] += mask * std::fabs(dmask_yup) * Mi * nfloor
                                * (SIGN(dmask_yup) * Cs - visheath) / penalty_timescale;
    }

    const auto iym = i.ym();
    const BoutReal dmask_ydown = mask - penalty_fa_data.mask[iym];
    if (std::fabs(dmask_ydown) > mask_threshold) {
      const BoutReal nfloor = BOUTMAX(Ni_fa[i], density_floor);
      const BoutReal tisheath = 0.5 * (Ti_fa[i] + Ti_fa[iym]);
      const BoutReal tesheath = 0.5 * (Te_fa[i] + Te_fa[iym]);
      const BoutReal visheath = 0.5 * (Vi_fa[i] + Vi_fa[iym]);
      const BoutReal Cs = sqrt((tesheath + tisheath) / Mi);

      momentum_penalty_fa[i] += mask * std::fabs(dmask_ydown) * Mi * nfloor
                                * (SIGN(dmask_ydown) * Cs - visheath) / penalty_timescale;
    }
  }

  return momentum_penalty_fa;
}

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

  diagnose =
      options["diagnose"].doc("Save penalty term diagnostics?").withDefault<bool>(false);

  gamma_e = options["gamma_e"]
                .doc("Electron sheath heat transmission coefficient")
                .withDefault(3.5);

  gamma_i =
      options["gamma_i"].doc("Ion sheath heat transmission coefficient").withDefault(3.5);

  penalty_timescale = options["penalty_timescale"]
                          .doc("Timescale of penalisation [seconds]")
                          .withDefault(1e-6)
                      / alloptions["units"]["seconds"].as<BoutReal>();

  surface_terms =
      options["surface_terms"].doc("Include surface terms?").withDefault<bool>(false);

  std::string mask_name = options["mask_name"]
                              .doc("Name of the mesh variable containing penalty mask")
                              .withDefault<std::string>("penalty_mask");

  Field3D penalty_mask;
  if (mesh->get(penalty_mask, mask_name) != 0) {
    throw BoutException("Could not read penalty mask variable '{}'", mask_name);
  }
  penalty_data = preparePenaltyMask(penalty_mask);

  if (surface_terms) {
    // Calculate surface terms using field-aligned coordinates
    penalty_fa_data =
        prepareFieldAlignedPenaltyMask(toFieldAligned(penalty_data.mask), *mesh);
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
  const bool has_phi = IS_SET_NOBOUNDARY(state["fields"]["phi"]);
  if (has_phi) {
    phi = getNoBoundary<Field3D>(state["fields"]["phi"]);
  }

  Field3D Te_fa;
  if (surface_terms) {
    // This is re-used in each species below
    Te_fa = toFieldAligned(Te);
  }

  {
    const Field3D Pe = electrons.isSet("pressure")
                           ? GET_NOBOUNDARY(Field3D, electrons["pressure"])
                           : Ne * Te;

    const Field3D Ve = electrons.isSet("velocity")
                           ? GET_NOBOUNDARY(Field3D, electrons["velocity"])
                           : zeroFrom(Ne);

    const Field3D NVe = electrons.isSet("momentum")
                            ? GET_NOBOUNDARY(Field3D, electrons["momentum"])
                            : Me * Ne * Ve;

    const Field3D density_source = electrons.isSet("density_source")
                                       ? getNonFinal<Field3D>(electrons["density_source"])
                                       : zeroFrom(Ne);

    const Field3D momentum_source =
        electrons.isSet("momentum_source")
            ? getNonFinal<Field3D>(electrons["momentum_source"])
            : zeroFrom(Ne);

    const Field3D energy_source = electrons.isSet("energy_source")
                                      ? getNonFinal<Field3D>(electrons["energy_source"])
                                      : zeroFrom(Ne);

    auto electron_penalty = calculateVolumetricPenalty(
        penalty_data, Ne, Te, Ve, Me, gamma_e, penalty_timescale, density_source,
        momentum_source, energy_source);

    if (surface_terms and has_phi) {
      // Surface penalty terms, to impose sheath current
      // The gradient of the mask gives the direction into the wall

      auto Ne_fa = toFieldAligned(Ne);
      auto Ve_fa = toFieldAligned(Ve);
      auto phi_fa = toFieldAligned(phi);

      auto momentum_penalty_fa = calculateElectronSurfaceMomentumPenalty(
          penalty_fa_data, Ne_fa, Te_fa, Ve_fa, phi_fa, Me, penalty_timescale);
      electron_penalty.momentum += fromFieldAligned(momentum_penalty_fa);
    }

    add(electrons["density_source"], electron_penalty.density);
    add(electrons["momentum_source"], electron_penalty.momentum);
    add(electrons["energy_source"], electron_penalty.energy);

    if (diagnose) {
      // Store penalty term diagnostics to used in outputVars
      auto& diagnostic_species = diagnostics["e"];
      set(diagnostic_species["density_penalty"], electron_penalty.density);
      set(diagnostic_species["momentum_penalty"], electron_penalty.momentum);
      set(diagnostic_species["energy_penalty"], electron_penalty.energy);
    }
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

    const Field3D Ni = GET_NOBOUNDARY(Field3D, species["density"]);
    const Field3D Ti = GET_NOBOUNDARY(Field3D, species["temperature"]);

    const Field3D Pi = species.isSet("pressure")
                           ? GET_NOBOUNDARY(Field3D, species["pressure"])
                           : Ni * Ti;

    const Field3D Vi = species.isSet("velocity")
                           ? GET_NOBOUNDARY(Field3D, species["velocity"])
                           : zeroFrom(Ni);

    const Field3D NVi = species.isSet("momentum")
                            ? GET_NOBOUNDARY(Field3D, species["momentum"])
                            : Mi * Ni * Vi;

    // Get the particle source to modify
    const Field3D density_source = species.isSet("density_source")
                                       ? getNonFinal<Field3D>(species["density_source"])
                                       : zeroFrom(Ni);

    const Field3D momentum_source = species.isSet("momentum_source")
                                        ? getNonFinal<Field3D>(species["momentum_source"])
                                        : zeroFrom(Ni);

    const Field3D energy_source = species.isSet("energy_source")
                                      ? getNonFinal<Field3D>(species["energy_source"])
                                      : zeroFrom(Ni);

    auto ion_penalty = calculateVolumetricPenalty(penalty_data, Ni, Ti, Vi, Mi, gamma_i,
                                                  penalty_timescale, density_source,
                                                  momentum_source, energy_source);

    if (surface_terms) {
      // Surface penalty terms.

      auto Ni_fa = toFieldAligned(Ni);
      auto Ti_fa = toFieldAligned(Ti);
      auto Vi_fa = toFieldAligned(Vi);

      auto momentum_penalty_fa = calculateIonSurfaceMomentumPenalty(
          penalty_fa_data, Ni_fa, Ti_fa, Te_fa, Vi_fa, Mi, penalty_timescale);
      ion_penalty.momentum += fromFieldAligned(momentum_penalty_fa);
    }

    add(species["density_source"], ion_penalty.density);
    add(species["momentum_source"], ion_penalty.momentum);
    add(species["energy_source"], ion_penalty.energy);

    set(species["density_penalty"], ion_penalty.density);
    set(species["momentum_penalty"], ion_penalty.momentum);
    set(species["energy_penalty"], ion_penalty.energy);

    if (diagnose) {
      // Store penalty term diagnostics to used in outputVars
      auto& diagnostic_species = diagnostics[kv.first];
      set(diagnostic_species["density_penalty"], ion_penalty.density);
      set(diagnostic_species["momentum_penalty"], ion_penalty.momentum);
      set(diagnostic_species["energy_penalty"], ion_penalty.energy);
    }
  }
}

void SheathBoundaryPenalty::outputVars(Options& state) {
  set_with_attrs(state["penalty_mask"], penalty_data.mask,
                 {{"units", ""},
                  {"long_name", "Penalty mask"},
                  {"standard_name", "Penalty mask"},
                  {"source", "sheath_boundary_penalty"}});

  if (diagnose) {
    // Normalisations
    auto Nnorm = get<BoutReal>(state["Nnorm"]);
    auto Tnorm = get<BoutReal>(state["Tnorm"]);
    const BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation
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
