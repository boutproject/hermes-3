#include "test_reaction_parser.hxx"

TEST(ReactionParserTest, MolecularIonFilter) {
  const ReactionParser parser{"d2 + d+ -> d2+ + d"};
  EXPECT_EQ(parser.get_species(species_filter::ion),
            (std::vector<std::string>{"d+", "d2+"}));
  EXPECT_EQ(parser.get_species(species_filter::neutral),
            (std::vector<std::string>{"d", "d2"}));
}

TEST_F(ReactionParserTIznTest, SpeciesFilters) {
  check_reactants_are({"t", "e"});
  check_products_are({"e", "t+"});
  check_heavy_species_are({"t", "t+"});
  check_heavy_products_are({"t+"});
}

TEST_F(ReactionParserTIznTest, PopChange) {
  check_pop_change("t", -1);
  check_pop_change("t+", 1);
  check_pop_change("e", 1);
}

TEST(ReactionParserInvalidTest, TooShort) {
  // Valid separator, but not long enough to be a valid reaction string
  EXPECT_THROW(ReactionParser{"->e"}, BoutException);
}

TEST(ReactionParserInvalidTest, WrongSeparator) {
  // Invalid separator
  EXPECT_THROW(ReactionParser{"t + e >- t+ + 2e"}, BoutException);
}

TEST(ReactionParserInvalidTest, NoReactants) {
  // Valid separator, but no reactants
  EXPECT_THROW(ReactionParser{"-> t+ + 2e"}, BoutException);
}

TEST(ReactionParserInvalidTest, NoProducts) {
  // Valid separator, but no products
  EXPECT_THROW(ReactionParser{"t + e ->"}, BoutException);
}
