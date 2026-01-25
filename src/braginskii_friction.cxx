#include <iterator>
#include <map>
#include <string>

#include <bout/bout_types.hxx>
#include <bout/constants.hxx>
#include <bout/field.hxx>
#include <bout/field3d.hxx>
#include <bout/options.hxx>
#include <fmt/format.h>

#include "../include/braginskii_friction.hxx"
#include "../include/component.hxx"

BraginskiiFriction::BraginskiiFriction(const std::string& name, Options& alloptions,
                                       Solver*)
    // FIXME: Not all species actually have collisions calculated
    : Component({readOnly("species:{all_species}:density"),
                 readIfSet("species:{all_species}:velocity", Regions::Interior),
                 readOnly("species:{all_species}:AA"),
                 readIfSet("species:{all_species}:charge"),
                 readIfSet("species:{all_species}:collision_frequencies:{all_species}_{"
                           "all_species2}_coll"),
                 readWrite("species:{all_species}:momentum_source")}) {
  AUTO_TRACE();
  Options& options = alloptions[name];
  frictional_heating = options["frictional_heating"]
                           .doc("Include R dot v heating term as energy source?")
                           .withDefault<bool>(true);
  diagnose =
      options["diagnose"].doc("Output additional diagnostics?").withDefault<bool>(false);

  if (frictional_heating) {
    setPermissions(readWrite("species:{all_species}:energy_source"));
  }
}

BoutReal momentumCoefficient(BoutReal Zi) {
  return Zi == 1 ? 0.51 : Zi == 2 ? 0.44 : Zi == 3 ? 0.40 : 0.38;
}

BoutReal momentumCoefficient(const std::string& name1, BoutReal Z1,
                             const std::string& name2, BoutReal Z2) {
  if (name1 == "e") {
    return momentumCoefficient(Z2);
  }
  if (name2 == "e") {
    return momentumCoefficient(Z1);
  }
  return 1.;
}

void BraginskiiFriction::transform_impl(GuardedOptions& state) {
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
    const Field3D velocity1 =
        species1.isSet("velocity") ? GET_NOBOUNDARY(Field3D, species1["velocity"]) : 0.0;
    ;

    const BoutReal A1 = GET_VALUE(BoutReal, species1["AA"]);
    const BoutReal Z1 =
        species1.isSet("charge") ? GET_VALUE(BoutReal, species1["charge"]) : 0;

    // Copy the iterator, so we don't iterate over the
    // lower half of the matrix, but start at the diagonal
    for (auto kv2 = kv1; kv2 != std::end(children); ++kv2) {
      // Can't have friction with oneself
      if (kv1->first == kv2->first) {
        continue;
      }

      GuardedOptions species2 = allspecies[kv2->first];

      // At least one of the species must have a velocity for there to be friction.
      if (!(isSetFinalNoBoundary(species1["velocity"])
            or isSetFinalNoBoundary(species2["velocity"]))) {
        continue;
      }

      const std::string coll_name = fmt::format("{}_{}_coll", kv1->first, kv2->first);
      // If collisions were not calculated between these two species, skip
      if (not species1["collision_frequencies"].isSet(coll_name)) {
        continue;
      }

      const Field3D nu = GET_VALUE(Field3D, species1["collision_frequencies"][coll_name]);
      const Field3D density2 = GET_VALUE(Field3D, species2["density"]);
      const Field3D velocity2 = species2.isSet("velocity")
                                    ? GET_NOBOUNDARY(Field3D, species2["velocity"])
                                    : 0.0;
      const BoutReal A2 = GET_VALUE(BoutReal, species2["AA"]);
      const BoutReal Z2 =
          species2.isSet("charge") ? GET_VALUE(BoutReal, species2["charge"]) : 0.;
      const BoutReal momentum_coefficient =
          momentumCoefficient(kv1->first, Z1, kv2->first, Z2);
      const Field3D F12 =
          momentum_coefficient * A1 * nu * density1 * (velocity2 - velocity1);

      add(species1["momentum_source"], F12);
      subtract(species2["momentum_source"], F12);

      // Diagnostics
      set(momentum_channels[species1.name()][species2.name()], F12);
      set(momentum_channels[species2.name()][species1.name()], -F12);

      if (frictional_heating) {
        // Heating due to friction and energy transfer
        //
        // In the pressure (thermal energy) equation we have a term
        // that transfers translational kinetic energy to thermal
        // energy, and an energy transfer between species:
        //
        // d/dt(3/2p_1) = ...  - F_12 v_1 + W_12
        //
        // The energy transfer term W_12 is chosen to make the
        // pressure change frame invariant:
        //
        // W_12 = (m_1 v_1 + m_2 v_2) / (m_1 + m_2) * F_12
        //
        // The sum of these two terms is:
        //
        // - F_12 v_1 + W_12 = m_2 (v_2  - v_1) / (m_1 + m_2) * F_12
        //
        // Note:
        //  1) This term is always positive: Collisions don't lead to cooling
        //  2) In the limit that m_2 << m_1 (e.g. electron-ion collisions),
        //     the lighter species is heated more than the heavy species.
        Field3D const species1_source = (A2 / (A1 + A2)) * (velocity2 - velocity1) * F12;
        Field3D const species2_source = (A1 / (A1 + A2)) * (velocity2 - velocity1) * F12;

        add(species1["energy_source"], species1_source);
        add(species2["energy_source"], species2_source);

        // Diagnostics
        set(friction_energy_channels[species1.name()][species2.name()], species1_source);
        set(friction_energy_channels[species2.name()][species1.name()], species2_source);
      }
    }
  }
}

void BraginskiiFriction::outputVars(Options& state) {
  AUTO_TRACE();

  if (!diagnose) {
    return; // Don't save diagnostics
  }

  // Normalisations
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  BoutReal const Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation
  auto Cs0 = get<BoutReal>(state["Cs0"]);

  /// Iterate through the first species in each collision pair
  for (const auto& [A, section] : momentum_channels.getChildren()) {
    for (const auto& [B, child] : section.getChildren()) {
      const std::string AB = A + B;

      // Frictional energy sources (both species heat through friction)
      if ((friction_energy_channels.isSection(A))
          and (friction_energy_channels[A].isSet(B))) {

        set_with_attrs(state[std::string("E") + AB + std::string("_coll_friction")],
                       getNonFinal<Field3D>(friction_energy_channels[A][B]),
                       {{"time_dimension", "t"},
                        {"units", "W / m^3"},
                        {"conversion", Pnorm * Omega_ci},
                        {"standard_name", AB + "frictional energy source"},
                        {"long_name", AB + "frictional energy source"},
                        {"species", A},
                        {"source", "collisions"}});
      }

      // Momentum exchange channels
      if ((momentum_channels.isSection(A)) and (momentum_channels[A].isSet(B))) {

        set_with_attrs(state[std::string("F") + AB + std::string("_coll")],
                       getNonFinal<Field3D>(friction_energy_channels[A][B]),
                       {{"time_dimension", "t"},
                        {"units", "kg m^-2 s^-2"},
                        {"conversion", SI::Mp * Nnorm * Cs0 * Omega_ci},
                        {"standard_name", AB + "collisional momentum transfer"},
                        {"long_name", AB + "collisional momentum transfer"},
                        {"species", A},
                        {"source", "collisions"}});
      }
    }
  }
}
