#include "test_reaction_parser.hxx"

/// Test the various species filters for the reaction "d2 + d+ -> d2+ + d"
TEST_F(ReactionParserMolCXTest, SpeciesFilters) {
  check_filter_gives(species_filter::consumed, {"d2", "d+"});
  check_filter_gives(species_filter::electron, {});
  check_filter_gives(species_filter::heavy, {"d", "d+", "d2", "d2+"});
  check_filter_gives(species_filter::ion, {"d+", "d2+"});
  check_filter_gives(species_filter::neutral, {"d", "d2"});
  check_filter_gives(species_filter::reactants, {"d2", "d+"});
  check_filter_gives(species_filter::produced, {"d", "d2+"});
  check_filter_gives(species_filter::products, {"d", "d2+"});
}

/// Test the various species filters for the reaction "t + e -> t+ + 2e"
TEST_F(ReactionParserTIznTest, SpeciesFilters) {
  check_filter_gives(species_filter::consumed, {"t"});
  check_filter_gives(species_filter::electron, {"e"});
  check_filter_gives(species_filter::heavy, {"t", "t+"});
  check_filter_gives(species_filter::ion, {"t+"});
  check_filter_gives(species_filter::neutral, {"t"});
  check_filter_gives(species_filter::reactants, {"t", "e"});
  check_filter_gives(species_filter::produced, {"e", "t+"});
  check_filter_gives(species_filter::products, {"e", "t+"});
}

/// Test population changes are as expected for the reaction "t + e -> t+ + 2e"
TEST_F(ReactionParserTIznTest, PopChange) {
  check_pop_change("t", -1);
  check_pop_change("t+", 1);
  check_pop_change("e", 1);
}

/// Reaction strings that do contain a separator, but aren't long enough to be valid reaction should throw
TEST(ReactionParserInvalidTest, TooShort) {
  EXPECT_THROW(ReactionParser{"->e"}, BoutException);
}

/// Reaction strings with an invalid separator should throw
TEST(ReactionParserInvalidTest, WrongSeparator) {
  // Invalid separator
  EXPECT_THROW(ReactionParser{"t + e >- t+ + 2e"}, BoutException);
}

/// Reaction strings with no reactants should throw
TEST(ReactionParserInvalidTest, NoReactants) {
  EXPECT_THROW(ReactionParser{"-> t+ + 2e"}, BoutException);
}

/// Reaction strings with no products should throw
TEST(ReactionParserInvalidTest, NoProducts) {
  EXPECT_THROW(ReactionParser{"t + e ->"}, BoutException);
}
