#ifndef TEST_REACTION_PARSER_H
#define TEST_REACTION_PARSER_H

#include <algorithm>

#include <gtest/gtest.h>

#include "reaction_parser.hxx"

std::string str_vec_to_str(std::vector<std::string> vec) {
  std::stringstream result;
  std::copy(vec.begin(), vec.end(), std::ostream_iterator<std::string>(result, ","));
  return result.str();
}

class ReactionParserTest : public ::testing::Test {
public:
  ReactionParserTest(std::string reaction_str) : parser(reaction_str){};

protected:
  ReactionParser parser;

  void check_heavy_products_are(std::vector<std::string> expected_species) {
    std::sort(expected_species.begin(), expected_species.end());
    std::vector<std::string> heavy_products =
        parser.get_species(species_filter::heavy, species_filter::products);
    std::sort(heavy_products.begin(), heavy_products.end());
    ASSERT_TRUE(std::equal(heavy_products.begin(), heavy_products.end(),
                           expected_species.begin(), expected_species.end()))
        << filter_check_err_msg("heavy products", expected_species, heavy_products);
  }
  void check_heavy_species_are(std::vector<std::string> expected_species) {
    std::sort(expected_species.begin(), expected_species.end());
    std::vector<std::string> heavy_species = parser.get_species(species_filter::heavy);
    std::sort(heavy_species.begin(), heavy_species.end());
    ASSERT_TRUE(std::equal(heavy_species.begin(), heavy_species.end(),
                           expected_species.begin(), expected_species.end()))
        << filter_check_err_msg("heavy species", expected_species, heavy_species);
  }

  void check_pop_change(std::string sp, int pop_change) {
    int result = parser.get_stoich().at(sp);
    ASSERT_EQ(result, pop_change)
        << std::endl
        << "Expected population change of species [" << sp << "] to be [" << pop_change
        << "], but ReactionParser computed [" << result << "]" << std::endl
        << "Reaction string was (" << parser.get_reaction_str() << ") ";
  }

  void check_products_are(std::vector<std::string> expected_species) {
    std::sort(expected_species.begin(), expected_species.end());
    std::vector<std::string> products = parser.get_species(species_filter::products);
    std::sort(products.begin(), products.end());
    ASSERT_TRUE(std::equal(products.begin(), products.end(), expected_species.begin(),
                           expected_species.end()))
        << filter_check_err_msg("products", expected_species, products);
  }

  void check_reactants_are(std::vector<std::string> expected_species) {
    std::sort(expected_species.begin(), expected_species.end());
    std::vector<std::string> reactants = parser.get_species(species_filter::reactants);
    std::sort(reactants.begin(), reactants.end());

    ASSERT_TRUE(std::equal(reactants.begin(), reactants.end(), expected_species.begin(),
                           expected_species.end()))
        << filter_check_err_msg("reactants", expected_species, reactants);
  }

private:
  std::string filter_check_err_msg(std::string lbl,
                                   std::vector<std::string> expected_species,
                                   std::vector<std::string> parsed_species) {
    std::stringstream ss;
    ss << std::endl
       << "Expected " << lbl << " to be" << std::endl
       << "  [" << str_vec_to_str(expected_species) << "]," << std::endl
       << "but ReactionParser found " << std::endl
       << "  [" << str_vec_to_str(parsed_species) << "]" << std::endl
       << "Reaction string was (" << parser.get_reaction_str() << ") ";
    return ss.str();
  }
};

class ReactionParserTIznTest : public ReactionParserTest {
public:
  ReactionParserTIznTest() : ReactionParserTest("t + e -> t+ + 2e") {}
};

#endif