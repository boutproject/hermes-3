#include <iterator>

#include <bout/constants.hxx>
#include <bout/field.hxx>
#include <bout/output_bout_types.hxx>

#include "../include/braginskii_heat_exchange.hxx"
#include "../include/hermes_utils.hxx"

BraginskiiHeatExchange::BraginskiiHeatExchange(std::string name, Options& alloptions,
                                               Solver*) {
  AUTO_TRACE();
  diagnose = alloptions[name]["diagnose"]
                 .doc("Output additional diagnostics?")
                 .withDefault<bool>(false);
}

void BraginskiiHeatExchange::transform(Options& state) {
  AUTO_TRACE();

  Options& allspecies = state["species"];

  // Iterate through all species
  // To avoid double counting, this needs to iterate over pairs
  // i.e. the diagonal and above
  //
  // Iterators kv1 and kv2 over the species map
  //
  //               kv2 ->
  //             species1  species2  species3
  // kv1   species1     X         X         X
  //  ||   species2               X         X
  //  \/   species3                         X
  //
  const std::map<std::string, Options>& children = allspecies.getChildren();
  for (auto kv1 = std::begin(children); kv1 != std::end(children); ++kv1) {
    Options& species1 = allspecies[kv1->first];
    // If collisions were not calculated for this species, skip it.
    if (not species1.isSection("collision_frequencies"))
      continue;

    const Field3D density1 = GET_NOBOUNDARY(Field3D, species1["density"]),
                  temperature1 = species1.isSet("temperature")
                                     ? GET_NOBOUNDARY(Field3D, species1["temperature"])
                                     : 0.0;

    const BoutReal A1 = GET_VALUE(BoutReal, species1["AA"]),
                   Z1 = species1.isSet("charge") ? GET_VALUE(BoutReal, species1["charge"])
                                                 : 0;

    // Copy the iterator, so we don't iterate over the
    // lower half of the matrix, but start at the diagonal
    for (std::map<std::string, Options>::const_iterator kv2 = kv1;
         kv2 != std::end(children); ++kv2) {
      // Can't have heat exchange with oneself
      if (kv1->first == kv2->first)
        continue;

      Options& species2 = allspecies[kv2->first];

      // At least one of the species must have a velocity for there to be friction.
      if (!(species1.isSet("temperature") or species2.isSet("temperature")))
        continue;

      const std::string coll_name =
          kv1->first + std::string("_") + kv2->first + std::string("_coll");
      // If collisions were not calculated between these two species, skip
      if (not species1["collision_frequencies"].isSet(coll_name))
        continue;

      const BoutReal A2 = GET_VALUE(BoutReal, species2["AA"]);

      const Field3D nu = GET_VALUE(Field3D, species1["collision_frequencies"][coll_name]),
                    temperature2 = species2.isSet("temperature")
                                       ? GET_NOBOUNDARY(Field3D, species2["temperature"])
                                       : 0.0;

      const Field3D Q12 =
          3 * (A1 / (A1 + A2)) * nu * density1 * (temperature2 - temperature1);

      add(species1["energy_source"], Q12);
      subtract(species2["energy_source"], Q12);

      // Diagnostics
      set(energy_channels[species1.name()][species2.name()], Q12);
      set(energy_channels[species2.name()][species1.name()], -Q12);
    }
  }
}

void BraginskiiHeatExchange::outputVars(Options& state) {
  AUTO_TRACE();

  if (!diagnose) {
    return; // Don't save diagnostics
  }

  // Normalisations
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation
  auto Cs0 = get<BoutReal>(state["Cs0"]);

  for (const auto& [A, section] : energy_channels) {
    for (const auto& [B, child] : section) {
      std::string AB = A + B;

      // Collisional energy transfer channels (i.e. thermal equilibration)
      if ((energy_channels.isSection(A)) and (energy_channels[A].isSet(B))) {

        set_with_attrs(state[std::string("E") + AB + std::string("_coll")],
                       getNonFinal<Field3D>(energy_channels[A][B]),
                       {{"time_dimension", "t"},
                        {"units", "W / m^3"},
                        {"conversion", Pnorm * Omega_ci},
                        {"standard_name", AB + "collisional energy transfer source"},
                        {"long_name", AB + "collisional energy transfer source"},
                        {"species", A},
                        {"source", "collisions"}});
      }
    }
  }
}
