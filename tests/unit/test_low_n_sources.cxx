#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "low_n_sources.hxx"
#include "species_parser.hxx"
#include "test_extras.hxx"
#include <bout/field_factory.hxx>
#include <bout/output.hxx>
#include <algorithm>

using ReactionDefs = std::vector<hermes::LowNSources::ReactionDef>;

namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

using namespace bout::globals;

// Fixed species names for all tests
const std::string atom_species = "d";
const std::string ion_species = "d+";
const std::string mol_species = "d2";

// Helper functions
namespace {

using LowNSourcesTest = FakeMeshFixture;

/**
 * @brief Generate an Options object with particular species density floors and other required params.
 *
 * @param atom_density_floor
 * @param ion_density_floor
 * @param molecule_density_floor
 * @return Options
 */
Options make_test_options(BoutReal atom_density_floor, BoutReal ion_density_floor,
                          BoutReal molecule_density_floor) {
  Options options;
  options["units"]["eV"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["inv_meters_cubed"] = 1.0;

  options["hermes"]["components"] = "e, low_n_sources, d, d+, d2";

  options[ion_species]["AA"] = 2.0;
  options[atom_species]["AA"] = 2.0;
  options[mol_species]["AA"] = 4.0;
  options["e"]["AA"] = 1.0 / 1836.0;

  options[ion_species]["velocity"] = 1.0;
  options[atom_species]["velocity"] = 2.0;
  options[mol_species]["velocity"] = 3.0;

  options[ion_species]["temperature"] = 3.0;
  options[atom_species]["temperature"] = 2.0;
  options[mol_species]["temperature"] = 1.0;

  options[atom_species]["density_floor"] = atom_density_floor;
  options[ion_species]["density_floor"] = ion_density_floor;
  options[mol_species]["density_floor"] = molecule_density_floor;

  options["low_n_sources"]["low_n_floor_fac"] = 100.0;
  options["low_n_sources"]["low_n_timescale"] = 1.0;

  // Tests all assume 'standard' rate strategt for now.
  options["low_n_sources"]["strategy"] = "standard";

  return options;
}

/**
 * @brief Set up a state with fixed, uniform atom, ion, molecule and electron densities
 *
 * @param state
 * @param n_atom
 * @param n_ion
 * @param n_molecule
 * @param n_electron
 */
void set_state(Options& state, BoutReal n_atom, BoutReal n_ion, BoutReal n_molecule,
               BoutReal n_electron = 0.0) {
  state["units"]["eV"] = 1.0;
  state["units"]["seconds"] = 1.0;
  state["units"]["inv_meters_cubed"] = 1.0;

  Options& species = state["species"];
  species[atom_species]["density"] =
      FieldFactory::get()->create3D(std::to_string(n_atom), &state, mesh);
  species[atom_species]["temperature"] =
      FieldFactory::get()->create3D("1.0", &state, mesh);
  species[atom_species]["velocity"] = FieldFactory::get()->create3D("2.0", &state, mesh);
  species[atom_species]["AA"] = 2.0;

  species[ion_species]["density"] =
      FieldFactory::get()->create3D(std::to_string(n_ion), &state, mesh);
  species[ion_species]["temperature"] =
      FieldFactory::get()->create3D("1.0", &state, mesh);
  species[ion_species]["velocity"] = FieldFactory::get()->create3D("1.0", &state, mesh);
  species[ion_species]["AA"] = 2.0;

  species[mol_species]["density"] =
      FieldFactory::get()->create3D(std::to_string(n_molecule), &state, mesh);
  species[mol_species]["temperature"] =
      FieldFactory::get()->create3D("1.0", &state, mesh);
  species[mol_species]["velocity"] = FieldFactory::get()->create3D("3.0", &state, mesh);
  species[mol_species]["AA"] = 4.0;

  species["e"]["density"] =
      FieldFactory::get()->create3D(std::to_string(n_electron), &state, mesh);
  species["e"]["temperature"] = FieldFactory::get()->create3D("1.0", &state, mesh);
  species["e"]["AA"] = 5.446e-4;
}

// Get the first non boundary value from a Field3D
BoutReal first_interior_value(const Field3D& field) {
  auto region = field.getRegion("RGN_NOBNDRY");
  return field[*region.begin()];
}

} // namespace

/**
 * @brief Test that LowNSources generates the expected reaction definitions
 * @todo Move to anonymous namespace?
 *
 * @param species_masses Map of species_name -> species_mass (as strings)
 * @param expected Expected vector of reaction definitions
 */
void test_reaction_def_creation(const std::map<std::string, std::string>& species_masses,
                                const ReactionDefs& expected) {

  // Set up species components
  // "AA" is currently the only way to identify a component as a species; must be set
  Options options;
  // PseudoReaction inherits Reaction, which reads unit normalisations at construction.
  options["units"]["eV"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["inv_meters_cubed"] = 1.0;
  std::string component_names_str = "e, low_n_sources";
  for (const auto& [species, mass] : species_masses) {
    component_names_str += ", " + species;
    options[species]["AA"] = mass;
  }
  options["hermes"]["components"] = component_names_str;

  // Construct the component (disable output to prevent values being echoed to stdout)
  output_info.disable();
  hermes::LowNSources component("low_n_sources", options);
  output_info.enable();

  // Retrieve generated pairs
  ReactionDefs result = component.get_reaction_defs();

  // Sort result, expected pairs before comparing
  std::sort(result.begin(), result.end());
  ReactionDefs expected_sorted = expected;
  std::sort(expected_sorted.begin(), expected_sorted.end());

  EXPECT_EQ(result, expected_sorted);
}

/**
 * @brief Test that LowNSources generates the expected reaction definitions for a variety of different species configurations
 *
 */
TEST(LowNSourcesConstructionTest, ConstructorBuildsOrderedPairsByElement) {

  // Ions, atom molecules for Deuterium, Tritium, Helium
  // => recombination, ionization, association, external_src for each isotope/element
  std::map<std::string, std::string> species_masses1 = {
      {"d+", "2.0"}, {"d", "2.0"},   {"d2", "4.0"}, {"t+", "3.0"}, {"t", "3.0"},
      {"t2", "6.0"}, {"he+", "2.0"}, {"he", "2.0"}, {"he2", "4.0"}};
  ReactionDefs expected_pairs1({{"d+", "d", PseudoReactionType::recombination},
                                {"d", "d+", PseudoReactionType::ionization},
                                {"d", "d2", PseudoReactionType::association},
                                {"d+", "d", PseudoReactionType::external_src},
                                {"t+", "t", PseudoReactionType::recombination},
                                {"t", "t+", PseudoReactionType::ionization},
                                {"t", "t2", PseudoReactionType::association},
                                {"t+", "t", PseudoReactionType::external_src},
                                {"he+", "he", PseudoReactionType::recombination},
                                {"he", "he+", PseudoReactionType::ionization},
                                {"he", "he2", PseudoReactionType::association},
                                {"he+", "he", PseudoReactionType::external_src}});
  test_reaction_def_creation(species_masses1, expected_pairs1);

  // Ions only; Deuterium, Tritium, Helium
  // => external sources for each species
  std::map<std::string, std::string> species_masses2 = {
      {"d+", "2.0"}, {"t+", "3.0"}, {"he+", "2.0"}};
  ReactionDefs expected_pairs2({{"d+", "d+", PseudoReactionType::external_src},
                                {"t+", "t+", PseudoReactionType::external_src},
                                {"he+", "he+", PseudoReactionType::external_src}});
  test_reaction_def_creation(species_masses2, expected_pairs2);

  // Ions, atoms and molecules for Deuterium, ions and atoms only for Tritium
  // => no association reaction for Tritium
  std::map<std::string, std::string> species_masses3 = {
      {"h+", "1.0"}, {"h", "1.0"}, {"h2", "2.0"}, {"t+", "3.0"}, {"t", "3.0"}};
  ReactionDefs expected_pairs3({{"h+", "h", PseudoReactionType::recombination},
                                {"h", "h+", PseudoReactionType::ionization},
                                {"h", "h2", PseudoReactionType::association},
                                {"h+", "h", PseudoReactionType::external_src},
                                {"t+", "t", PseudoReactionType::recombination},
                                {"t", "t+", PseudoReactionType::ionization},
                                {"t+", "t", PseudoReactionType::external_src}});
  test_reaction_def_creation(species_masses3, expected_pairs3);

  // Possibly pathological cases
  // 1. no atomic species, only ions and molecules
  // => external source for each ion species
  std::map<std::string, std::string> species_masses4 = {
      {"h+", "1.0"}, {"h2", "2.0"}, {"d+", "2.0"}, {"d2", "4.0"}};
  ReactionDefs expected_pairs4({{"h+", "h+", PseudoReactionType::external_src},
                                {"d+", "d+", PseudoReactionType::external_src}});
  test_reaction_def_creation(species_masses4, expected_pairs4);

  // 2. Helium atoms and ions, Deuterium ions only
  // => external ion source only for Deuterium
  // => external ion source + recombination + ionization for Helium
  std::map<std::string, std::string> species_masses5 = {
      {"he+", "2.0"}, {"he", "2.0"}, {"d+", "2.0"}};
  ReactionDefs expected_pairs5({{"he+", "he", PseudoReactionType::recombination},
                                {"he", "he+", PseudoReactionType::ionization},
                                {"he+", "he", PseudoReactionType::external_src},
                                {"d+", "d+", PseudoReactionType::external_src}});
  test_reaction_def_creation(species_masses5, expected_pairs5);
}

/**
 * @brief Test that sources are zero when densities are all above thresholds
 */
TEST_F(LowNSourcesTest, NoSrcsAboveThreshold) {
  // Set the same floor for all species
  const BoutReal nfloor = 1e-7;
  Options options = make_test_options(nfloor, nfloor, nfloor);
  hermes::LowNSources component("low_n_sources", options);

  // Set all densities well above the threshold and call transform
  const BoutReal nth =
      nfloor * options["low_n_sources"]["low_n_floor_fac"].as<BoutReal>();
  Options state;
  set_state(state, 10 * nth, 10 * nth, 10 * nth);
  component.transform(state);

  // Check that all density source terms are zero
  for (const auto& species : {atom_species, ion_species, mol_species}) {
    EXPECT_TRUE(IsFieldEqual(get<Field3D>(state["species"][species]["density_source"]),
                             0.0, "RGN_NOBNDRY"));
  }
}

/**
 * @brief Check that the (ion) density source increases in magnitude the further you get below the density threshold.
 * @details Sets up two states with varying degrees of ion deficits and one with no deficit, then compare the source sizes.
 * Also implicitly checks that densities of exactly zero trigger the sources.
 *
 */
TEST_F(LowNSourcesTest, LargerDensityDeficitGivesLargerSource) {
  // Configuration with single density floor for all species
  const BoutReal nfloor = 1e-7;
  Options options = make_test_options(nfloor, nfloor, nfloor);
  hermes::LowNSources component("low_n_sources", options);

  // Set (uniform) atomic, molecular densities well above their respective thresholds
  const BoutReal nth =
      nfloor * options["low_n_sources"]["low_n_floor_fac"].as<BoutReal>();
  const BoutReal natom = 10 * nth;
  const BoutReal nmol = 10 * nth;

  // Low ion density case
  Options state_low_nion;
  const BoutReal low_nion = 0.5 * nth;
  set_state(state_low_nion, natom, low_nion, nmol);
  component.transform(state_low_nion);
  const Field3D src_low_nion =
      get<Field3D>(state_low_nion["species"][ion_species]["density_source"]);

  // Very low ion density case
  Options state_vlow_nion;
  const BoutReal vlow_nion = 0.1 * nth;
  set_state(state_vlow_nion, natom, vlow_nion, nmol);
  component.transform(state_vlow_nion);
  const Field3D src_vlow_nion =
      get<Field3D>(state_vlow_nion["species"][ion_species]["density_source"]);

  // The source term should be larger for the state with the very low ion density
  EXPECT_GT(first_interior_value(src_vlow_nion), first_interior_value(src_low_nion));

  // Zero ion density case
  Options state_zero_nion;
  const BoutReal zero_nion = 0.0;
  set_state(state_zero_nion, natom, zero_nion, nmol);
  component.transform(state_zero_nion);
  const Field3D src_zero_nion =
      get<Field3D>(state_zero_nion["species"][ion_species]["density_source"]);

  // The logic should still apply at ZERO ion density (although the boost factor is limited to avoid div zero)
  EXPECT_GT(first_interior_value(src_zero_nion), first_interior_value(src_vlow_nion));
}

/**
 * @brief Check that source sizes are as expected when the atom and molecule densities are low
 *  (i.e. when the external source is off, and particle number is conserved.)
 */
TEST_F(LowNSourcesTest, ParticleConservingSourceMagnitudes) {
  // Set up a component; same floor for all species
  const BoutReal nfloor = 1e-7;
  Options options = make_test_options(nfloor, nfloor, nfloor);
  hermes::LowNSources component("low_n_sources", options);

  // Same density threshold for all species
  const BoutReal nth =
      nfloor * options["low_n_sources"]["low_n_floor_fac"].as<BoutReal>();

  // Set up state with ions above threshold, atoms and molecules below
  const BoutReal natom = 0.5 * nth;
  const BoutReal nmol = 0.1 * nth;
  const BoutReal nion = 2 * nth;
  Options state;
  set_state(state, natom, nion, nmol);

  // Expected source sizes.
  // Note that target species population changes are +1 for all pseudo-reactions, so density_src == rate
  const BoutReal low_n_timescale =
      options["low_n_sources"]["low_n_timescale"].as<BoutReal>();
  // Recombination restores atoms, boosted the further molecules fall below their own
  // threshold (boost factor clamped to avoid divide-by-zero, as in the production code)
  const BoutReal mol_thresh_ratio = nmol / nth;
  const BoutReal recomb_boost_fac = 1.0 / std::clamp(mol_thresh_ratio, 0.01, 1.0);
  const BoutReal recomb_rate = ((nth - natom) / low_n_timescale) * recomb_boost_fac;
  // Ions are depleted by recombination
  const BoutReal expected_ion_source = -recomb_rate;
  // Molecules are boosted by association; source only depends on molecule deficit
  const BoutReal expected_mol_source = (nth - nmol) / low_n_timescale;
  // Atoms are boosted by recombination, but depleted by association - need to subtract the molecule source
  //  (association consumes two atoms to make one molecule, so the atom loss is 2x the molecule source)
  const BoutReal atom_loss_to_association = 2 * expected_mol_source;
  BoutReal expected_atom_source = recomb_rate - atom_loss_to_association;

  // Transform and extract sources
  component.transform(state);
  const Field3D ion_source =
      get<Field3D>(state["species"][ion_species]["density_source"]);
  const Field3D atom_source =
      get<Field3D>(state["species"][atom_species]["density_source"]);
  const Field3D mol_source =
      get<Field3D>(state["species"][mol_species]["density_source"]);

  // Check source sizes at first interior point (sources are spatially uniform by construction)
  EXPECT_DOUBLE_EQ(first_interior_value(ion_source), expected_ion_source);
  EXPECT_DOUBLE_EQ(first_interior_value(mol_source), expected_mol_source);
  EXPECT_DOUBLE_EQ(first_interior_value(atom_source), expected_atom_source);
}

/**
 * @brief Check the sizes of momentum and energy sources generated by LowNSources.
 * Only trigger association, for now, to keep things simple.
 *
 */
TEST_F(LowNSourcesTest, MomentumEnergySources) {
  // Set up a component; same floor for all species
  const BoutReal nfloor = 1e-7;
  Options options = make_test_options(nfloor, nfloor, nfloor);
  hermes::LowNSources component("low_n_sources", options);

  // Same density threshold for all species
  const BoutReal nth =
      nfloor * options["low_n_sources"]["low_n_floor_fac"].as<BoutReal>();

  // Set up state with only molecules below threshold; only association will be triggered
  const BoutReal natom = 2 * nth;
  const BoutReal nmol = 0.5 * nth;
  const BoutReal nion = 2 * nth;
  Options state;
  set_state(state, natom, nion, nmol);

  // Calc expected source sizes
  BoutReal expected_rate =
      (nth - nmol) / options["low_n_sources"]["low_n_timescale"].as<BoutReal>();
  const int atom_pop_change = -2;
  const int mol_pop_change = 1;
  const BoutReal atom_mass = get<BoutReal>(state["species"][atom_species]["AA"]);
  const Field3D atom_vel = get<Field3D>(state["species"][atom_species]["velocity"]);
  const Field3D atom_temp = get<Field3D>(state["species"][atom_species]["temperature"]);
  const Field3D mol_vel = get<Field3D>(state["species"][mol_species]["velocity"]);
  const Field3D mol_temp = get<Field3D>(state["species"][mol_species]["temperature"]);

  BoutReal expected_atom_density_src = atom_pop_change * expected_rate;
  BoutReal expected_atom_momentum_src =
      expected_atom_density_src * atom_mass * first_interior_value(atom_vel);
  BoutReal expected_atom_energy_src =
      expected_atom_density_src * 1.5 * first_interior_value(atom_temp);
  BoutReal expected_mol_density_src = mol_pop_change * expected_rate;
  BoutReal expected_mol_momentum_src = -expected_atom_momentum_src;
  BoutReal expected_mol_energy_src = -expected_atom_energy_src;

  // Transform and extract atomic and molecular density, momentum, energy sources
  component.transform(state);
  const Field3D atom_density_src =
      get<Field3D>(state["species"][atom_species]["density_source"]);
  const Field3D atom_momentum_src =
      get<Field3D>(state["species"][atom_species]["momentum_source"]);
  const Field3D atom_energy_src =
      get<Field3D>(state["species"][atom_species]["energy_source"]);
  const Field3D mol_density_src =
      get<Field3D>(state["species"][mol_species]["density_source"]);
  const Field3D mol_momentum_src =
      get<Field3D>(state["species"][mol_species]["momentum_source"]);
  const Field3D mol_energy_src =
      get<Field3D>(state["species"][mol_species]["energy_source"]);

  // None of the sources should be zero
  EXPECT_FALSE(IsFieldEqual(atom_density_src, 0.0, "RGN_NOBNDRY"));
  EXPECT_FALSE(IsFieldEqual(atom_momentum_src, 0.0, "RGN_NOBNDRY"));
  EXPECT_FALSE(IsFieldEqual(atom_energy_src, 0.0, "RGN_NOBNDRY"));
  EXPECT_FALSE(IsFieldEqual(mol_density_src, 0.0, "RGN_NOBNDRY"));
  EXPECT_FALSE(IsFieldEqual(mol_momentum_src, 0.0, "RGN_NOBNDRY"));
  EXPECT_FALSE(IsFieldEqual(mol_energy_src, 0.0, "RGN_NOBNDRY"));

  // Check source sizes at first interior point (sources are spatially uniform by construction)
  EXPECT_DOUBLE_EQ(first_interior_value(atom_density_src), expected_atom_density_src);
  EXPECT_DOUBLE_EQ(first_interior_value(atom_momentum_src), expected_atom_momentum_src);
  EXPECT_DOUBLE_EQ(first_interior_value(atom_energy_src), expected_atom_energy_src);
  EXPECT_DOUBLE_EQ(first_interior_value(mol_density_src), expected_mol_density_src);
  EXPECT_DOUBLE_EQ(first_interior_value(mol_momentum_src), expected_mol_momentum_src);
  EXPECT_DOUBLE_EQ(first_interior_value(mol_energy_src), expected_mol_energy_src);
}

/**
 * @brief Check that source sizes are as expected when the ion density is low, particularly the effect of the external ion source.
 * @details Three scenarios:
 * 1. nion is low, natom is above its threshold (no external source, ionisation boosts ions)
 * 2. nion is low, natom is somewhat below its threshold (external source of ions enabled, ionisation off, recombination boosts atoms)
 * 3. nion is low, natom is further below its threshold (stronger external source of ions, ionisation off, recombination boosts atoms)
 */
TEST_F(LowNSourcesTest, EffectOfExternalIonSource) {
  // Set up a component; same floor for all species
  const BoutReal nfloor = 1e-7;
  Options options = make_test_options(nfloor, nfloor, nfloor);
  hermes::LowNSources component("low_n_sources", options);

  // Same density threshold for all species
  const BoutReal nth =
      nfloor * options["low_n_sources"]["low_n_floor_fac"].as<BoutReal>();

  // Ion density is always below threshold, molecules always above.
  const BoutReal nion = 1e-2 * nth;
  const BoutReal nmol = 10 * nth;

  // ==========================================================================
  // Transform a state where the atomic density is above its threshold and extract sources
  const BoutReal natom = 10 * nth;
  Options no_natom_deficit_result;
  set_state(no_natom_deficit_result, natom, nion, nmol);
  component.transform(no_natom_deficit_result);
  const Field3D no_natom_deficit_ion_src =
      get<Field3D>(no_natom_deficit_result["species"][ion_species]["density_source"]);
  const Field3D no_natom_deficit_atom_src =
      get<Field3D>(no_natom_deficit_result["species"][atom_species]["density_source"]);

  // Only (pseudo-)ionisation should be active
  // => Expect -ve atom source
  EXPECT_LT(first_interior_value(no_natom_deficit_atom_src), 0.0);
  // => Expect +ve ion source with the same magnitude
  EXPECT_DOUBLE_EQ(first_interior_value(no_natom_deficit_ion_src),
                   -1 * first_interior_value(no_natom_deficit_atom_src));

  // ==========================================================================
  // Transform a state where the atomic density is 10x lower than the threshold and extract sources
  const BoutReal small_deficit_natom = 1e-1 * nth;
  Options small_natom_deficit_result;
  set_state(small_natom_deficit_result, small_deficit_natom, nion, nmol);
  component.transform(small_natom_deficit_result);
  const Field3D small_natom_deficit_ion_src =
      get<Field3D>(small_natom_deficit_result["species"][ion_species]["density_source"]);
  const Field3D small_natom_deficit_atom_src =
      get<Field3D>(small_natom_deficit_result["species"][atom_species]["density_source"]);

  // Both (pseudo-)recombination and the external ion source should be active
  // => Expect +ve atom source
  EXPECT_GT(first_interior_value(small_natom_deficit_atom_src), 0.0);

  // Expect +ve ion source too, because external source has been boosted to facilitate recombination
  EXPECT_GT(first_interior_value(small_natom_deficit_ion_src), 0.0);

  // ==========================================================================
  // Transform a state where the atomic density is 100x lower than the threshold and extract sources
  const BoutReal large_deficit_natom = 1e-2 * nth;
  Options large_natom_deficit_result;
  set_state(large_natom_deficit_result, large_deficit_natom, nion, nmol);
  component.transform(large_natom_deficit_result);
  const Field3D large_natom_deficit_ion_src =
      get<Field3D>(large_natom_deficit_result["species"][ion_species]["density_source"]);
  const Field3D large_natom_deficit_atom_src =
      get<Field3D>(large_natom_deficit_result["species"][atom_species]["density_source"]);

  // Both (pseudo-)recombination and the external ion source should be active.
  // The external source should be boosted relative to the previous case.
  // => Expect +ve atom source
  EXPECT_GT(first_interior_value(large_natom_deficit_atom_src), 0.0);

  // => Expect ion source to be larger than in small atomic deficit case
  EXPECT_GT(first_interior_value(large_natom_deficit_ion_src),
            first_interior_value(small_natom_deficit_ion_src));
}
