#include "test_reaction_parser.hxx"

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