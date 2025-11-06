#include <iterator>
#include <map>
#include <string>

#include <bout/bout_types.hxx>
#include <bout/constants.hxx>
#include <bout/field.hxx>
#include <bout/field3d.hxx>
#include <bout/options.hxx>
#include <fmt/format.h>

#include "../include/braginskii_heat_exchange.hxx"
#include "../include/component.hxx"

BraginskiiHeatExchange::BraginskiiHeatExchange(const std::string& name,
                                               Options& alloptions, Solver*)
    // FIXME: Not all species are actually read or written; only those with collision
    // rates and temperatures
    : Component({readOnly("species:{all_species}:{input_vars}"),
                 readIfSet("species:{all_species}:{optional_vars}"),
                 readWrite("species:{all_species}:{output_vars}")}) {
  AUTO_TRACE();
  diagnose = alloptions[name]["diagnose"]
                 .doc("Output additional diagnostics?")
                 .withDefault<bool>(false);
  state_variable_access.substitute("input_vars", {"AA", "density"});
  // FIXME: We don't access the self-collision rate
  state_variable_access.substitute(
      "optional_vars",
      {"charge", "collision_frequencies:{all_species}_{all_species2}_coll",
       "temperature"});
  state_variable_access.substitute("output_vars", {"momentum_source", "energy_source"});
}

void BraginskiiHeatExchange::transform_impl(GuardedOptions& state) {
  AUTO_TRACE();

  GuardedOptions allspecies = state["species"];

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
  const std::map<std::string, GuardedOptions> children = allspecies.getChildren();
  for (auto kv1 = std::begin(children); kv1 != std::end(children); ++kv1) {
    GuardedOptions species1 = allspecies[kv1->first];
    // If collisions were not calculated for this species, skip it.
    if (not species1.isSection("collision_frequencies")) {
      continue;
    }

    const Field3D density1 = GET_NOBOUNDARY(Field3D, species1["density"]);
    const Field3D temperature1 = species1.isSet("temperature")
                                     ? GET_NOBOUNDARY(Field3D, species1["temperature"])
                                     : 0.0;

    const BoutReal A1 = GET_VALUE(BoutReal, species1["AA"]);

    // Copy the iterator, so we don't iterate over the
    // lower half of the matrix, but start at the diagonal
    for (std::map<std::string, GuardedOptions>::const_iterator kv2 = kv1;
         kv2 != std::end(children); ++kv2) {
      // Can't have heat exchange with oneself
      if (kv1->first == kv2->first) {
        continue;
      }

      GuardedOptions species2 = allspecies[kv2->first];

      // At least one of the species must have a temperature for there to be heat
      // exchange.
      if (!(species1.isSet("temperature") or species2.isSet("temperature"))) {
        continue;
      }

      const std::string coll_name = fmt::format("{}_{}_coll", kv1->first, kv2->first);
      // If collisions were not calculated between these two species, skip
      if (not species1["collision_frequencies"].isSet(coll_name)) {
        continue;
      }

      const BoutReal A2 = GET_VALUE(BoutReal, species2["AA"]);

      const Field3D nu = GET_VALUE(Field3D, species1["collision_frequencies"][coll_name]);
      const Field3D temperature2 = species2.isSet("temperature")
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
  BoutReal const Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

  for (const auto& [A, section] : energy_channels) {
    for (const auto& [B, child] : section) {
      const std::string AB = A + B;

      // Collisional energy transfer channels (i.e. thermal equilibration)
      if ((energy_channels.isSection(A)) and (energy_channels[A].isSet(B))) {

        set_with_attrs(state[fmt::format("E{}_coll", AB)],
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
