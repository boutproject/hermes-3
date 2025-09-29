#include <iterator>

#include <bout/constants.hxx>
#include <bout/field.hxx>
#include <bout/output_bout_types.hxx>

#include "../include/hermes_utils.hxx"
#include "../include/braginskii_friction.hxx"

BraginskiiFriction::BraginskiiFriction(std::string name, Options& alloptions, Solver*) {
  AUTO_TRACE();
  Options& options = alloptions[name];
  frictional_heating = options["frictional_heating"]
    .doc("Include R dot v heating term as energy source?")
    .withDefault<bool>(true);
  diagnose =
      options["diagnose"].doc("Output additional diagnostics?").withDefault<bool>(false);
}

BoutReal momentumCoefficient(BoutReal Zi) {
  return Zi == 1 ? 0.51 : Zi == 2 ? 0.44 : Zi == 3 ? 0.40 : 0.38;
}

BoutReal momentumCoefficient(std::string name1, BoutReal Z1, std::string name2,
                             BoutReal Z2) {
  if (name1 == "e") {
    return momentumCoefficient(Z2);
  } else if (name2 == "e") {
    return momentumCoefficient(Z1);
  } else {
    return 1.;
  }
}

void BraginskiiFriction::transform(Options& state) {
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
    if (not species1.isSection("collision_frequencies")) continue;

    const Field3D density1 = GET_NOBOUNDARY(Field3D, species1["density"]),
      velocity1 = species1.isSet("velocity") ? GET_NOBOUNDARY(Field3D, species1["velocity"]) : 0.0;;

    const BoutReal A1 = GET_VALUE(BoutReal, species1["AA"]), Z1 = species1.isSet("charge") ? GET_VALUE(BoutReal, species1["charge"]) : 0;

    // Copy the iterator, so we don't iterate over the
    // lower half of the matrix, but start at the diagonal
    for (std::map<std::string, Options>::const_iterator kv2 = kv1;
         kv2 != std::end(children); ++kv2) {
      // Can't have friction with oneself
      if (kv1->first == kv2->first) continue;

      Options& species2 = allspecies[kv2->first];

      // At least one of the species must have a velocity for there to be friction.
      if (!(isSetFinalNoBoundary(species1["velocity"]) or
            isSetFinalNoBoundary(species2["velocity"]))) continue;
      
      const std::string coll_name = kv1->first + std::string("_") + kv2->first + std::string("_coll");
      // If collisions were not calculated between these two species, skip
      if (not species1["collision_frequencies"].isSet(coll_name)) continue;

      const Field3D nu = GET_VALUE(Field3D, species1["collision_frequencies"][coll_name]),
        density2 = GET_VALUE(Field3D, species2["density"]),
        velocity2 = species2.isSet("velocity") ? GET_NOBOUNDARY(Field3D, species2["velocity"]) : 0.0;
      const BoutReal A2 = GET_VALUE(BoutReal, species2["AA"]), Z2 = species2.isSet("charge") ? GET_VALUE(BoutReal, species2["charge"]) : 0., momentum_coefficient = momentumCoefficient(kv1->first, Z1, kv2->first, Z2);
      const Field3D F12 = momentum_coefficient * A1 * nu * density1 * (velocity2 - velocity1);

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
        Field3D species1_source = (A2 / (A1 + A2)) * (velocity2 - velocity1) * F12;
        Field3D species2_source = (A1 / (A1 + A2)) * (velocity2 - velocity1) * F12;

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
  BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation
  auto Cs0 = get<BoutReal>(state["Cs0"]);

  /// Iterate through the first species in each collision pair
  const std::map<std::string, Options>& level1 = momentum_channels.getChildren();
  for (auto s1 = std::begin(level1); s1 != std::end(level1); ++s1) {
    const Options& section = momentum_channels[s1->first];

    /// Iterate through the second species in each collision pair
    const std::map<std::string, Options>& level2 = section.getChildren();
    for (auto s2 = std::begin(level2); s2 != std::end(level2); ++s2) {

      std::string A = s1->first;
      std::string B = s2->first;
      std::string AB = A + B;

      // Frictional energy sources (both species heat through friction)
      if ((friction_energy_channels.isSection(A)) and (friction_energy_channels[A].isSet(B))) {

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
