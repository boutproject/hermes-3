#ifndef TEST_REACTION_PARSER_H
#define TEST_REACTION_PARSER_H

#include <algorithm>

#include <gtest/gtest.h>

#include "reaction_parser.hxx"

class ReactionParserTest : public ::testing::Test {
public:
  ReactionParserTest(std::string reaction_str) : parser(reaction_str) {};

protected:
  ReactionParser parser;

  /// Check that a species filter returns the expected names, independently of order.
  void check_filter_gives(species_filter type,
                          std::vector<std::string> expected_species) {
    std::sort(expected_species.begin(), expected_species.end());
    auto filtered_species = parser.get_species(type);
    std::sort(filtered_species.begin(), filtered_species.end());
    ASSERT_EQ(filtered_species, expected_species)
        << "Species filter: " << static_cast<int>(type) << '\n'
        << "Reaction string: " << parser.get_reaction_str();
  }

  /// Helper function to check that the population change determined for a particular species matches an expected value
  void check_pop_change(std::string sp, int pop_change) {
    int result = parser.pop_change(sp);
    ASSERT_EQ(result, pop_change)
        << std::endl
        << "Expected population change of species [" << sp << "] to be [" << pop_change
        << "], but ReactionParser computed [" << result << "]" << std::endl
        << "Reaction string was (" << parser.get_reaction_str() << ") ";
  }
};

class ReactionParserTIznTest : public ReactionParserTest {
public:
  ReactionParserTIznTest() : ReactionParserTest("t + e -> t+ + 2e") {}
};

class ReactionParserMolCXTest : public ReactionParserTest {
public:
  ReactionParserMolCXTest() : ReactionParserTest("d2 + d+ -> d2+ + d") {}
};

#endif
