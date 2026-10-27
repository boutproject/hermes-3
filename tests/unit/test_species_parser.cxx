#include <bout/utils.hxx>
#include <gtest/gtest.h>

#include "species_parser.hxx"

/// Test fixture
class SpeciesParserTest : public ::testing::Test {};

/// @brief Check extraction of the underlying element symbol.
TEST_F(SpeciesParserTest, GetElement) {
  EXPECT_EQ(SpeciesParser("e").get_element(), "e");
  EXPECT_EQ(SpeciesParser("He+2").get_element(), "he");
  EXPECT_EQ(SpeciesParser("d2+").get_element(), "d");
  EXPECT_EQ(SpeciesParser("t2").get_element(), "t");
}

/// @brief Check species types identified directly from species strings.
TEST_F(SpeciesParserTest, GetTypeFromString) {
  EXPECT_EQ(SpeciesParser::get_type("e"), SpeciesType::electron);
  EXPECT_EQ(SpeciesParser::get_type("d+"), SpeciesType::ion);
  EXPECT_EQ(SpeciesParser::get_type("d2+"), SpeciesType::molecular_ion);
  EXPECT_EQ(SpeciesParser::get_type("d"), SpeciesType::neutral);
  EXPECT_EQ(SpeciesParser::get_type("d2"), SpeciesType::molecule);
  // Recognise legacy generic ion name
  EXPECT_EQ(SpeciesParser::get_type("i"), SpeciesType::ion);
}

/// @brief Check species type predicates work as expected.
TEST_F(SpeciesParserTest, TypePredicates) {
  // Test cases
  const struct {
    const char* name;
    bool is_ion;
    bool is_neutral;
  } test_cases[] = {{"e", false, false}, {"e-", false, false}, {"i", true, false},
                    {"d+", true, false}, {"d2+", true, false}, {"d2-", true, false},
                    {"d", false, true},  {"d2", false, true}};
  // Check that predicates return the expected result for each test case.
  for (const auto& test_case : test_cases) {
    const SpeciesParser species(test_case.name);
    const SpeciesType type = species.get_type();
    EXPECT_EQ(type, SpeciesParser::get_type(test_case.name)) << test_case.name;
    EXPECT_EQ(isIonSpeciesType(type), test_case.is_ion);
    EXPECT_EQ(isNeutralSpeciesType(type), test_case.is_neutral);
  }
}

/// @brief Check that constructor throws for a variety of invalid strings
TEST_F(SpeciesParserTest, InvalidStrings) {
  // Empty string
  ASSERT_THROW(SpeciesParser(""), BoutException);
  // Only numbers
  ASSERT_THROW(SpeciesParser("123"), BoutException);
  // Only a charge sign without a base species
  ASSERT_THROW(SpeciesParser("+"), BoutException);
  // Only a charge sign without a base species
  ASSERT_THROW(SpeciesParser("+2"), BoutException);
  // Charge sign before element
  ASSERT_THROW(SpeciesParser("-h"), BoutException);
  // Invalid characters
  ASSERT_THROW(SpeciesParser("%&?!"), BoutException);
  // Invalid characters
  ASSERT_THROW(SpeciesParser("he|2"), BoutException);
  // Space before charge specifier
  ASSERT_THROW(SpeciesParser("he +"), BoutException);
  // Space before charge specifier
  ASSERT_THROW(SpeciesParser("H -"), BoutException);
  // Double charge sign
  ASSERT_THROW(SpeciesParser("d2++"), BoutException);
  // Invalid molecule count of zero
  ASSERT_THROW(SpeciesParser("d0"), BoutException);
  // Invalid molecule count (leading zero)
  ASSERT_THROW(SpeciesParser("d01"), BoutException);
  // Base species can't be more than 2 characters
  ASSERT_THROW(SpeciesParser("abc"), BoutException);
  // Invalid charge specifier
  ASSERT_THROW(SpeciesParser("he+0"), BoutException);
}

/// @brief Check ionisation of a highly-ionised species.
TEST_F(SpeciesParserTest, IoniseHighlyIonised) {
  SpeciesParser highly("ne+8");
  SpeciesParser next = highly.ionised();
  EXPECT_EQ(next.get_base_species(), "ne");
  EXPECT_EQ(next.get_charge(), 9);
  EXPECT_EQ(next.get_str(), "ne+9");
}

/// @brief Check ionisation of a neutral species.
TEST_F(SpeciesParserTest, IoniseNeutral) {
  SpeciesParser neutral("h");
  SpeciesParser ionised = neutral.ionised();
  EXPECT_EQ(ionised.get_base_species(), "h");
  EXPECT_EQ(ionised.get_charge(), 1);
  EXPECT_EQ(ionised.get_str(), "h+");
}

/// @brief Check ionisation of a singly-ionised species.
TEST_F(SpeciesParserTest, IoniseSinglyIonised) {
  SpeciesParser singly("li+");
  SpeciesParser doubly = singly.ionised();
  EXPECT_EQ(doubly.get_base_species(), "li");
  EXPECT_EQ(doubly.get_charge(), 2);
  EXPECT_EQ(doubly.get_str(), "li+2");
}

/// @brief Check that some species strings are correctly identified as molecules/non-molecules.
TEST_F(SpeciesParserTest, IsMolecule) {
  EXPECT_FALSE(SpeciesParser("d").is_molecule());
  EXPECT_FALSE(SpeciesParser("d+").is_molecule());
  EXPECT_FALSE(SpeciesParser("e").is_molecule());

  EXPECT_TRUE(SpeciesParser("d2").is_molecule());
  EXPECT_TRUE(SpeciesParser("d2+").is_molecule());
  EXPECT_TRUE(SpeciesParser("he3+").is_molecule());
}

/// @brief Check whether species are identified as neutral from their charge.
TEST_F(SpeciesParserTest, IsNeutral) {
  EXPECT_TRUE(SpeciesParser("d").is_neutral());
  EXPECT_TRUE(SpeciesParser("d2").is_neutral());
  EXPECT_FALSE(SpeciesParser("d+").is_neutral());
  EXPECT_FALSE(SpeciesParser("d2+").is_neutral());
  EXPECT_FALSE(SpeciesParser("e").is_neutral());
  EXPECT_FALSE(SpeciesParser("t-").is_neutral());
}

/// @brief Check multiple ionisations
TEST_F(SpeciesParserTest, MultipleIonisations) {
  // Start with neutral He, ionise twice
  SpeciesParser he("he");
  EXPECT_EQ(he.get_charge(), 0);

  SpeciesParser hep1 = he.ionised();
  EXPECT_EQ(hep1.get_charge(), 1);
  EXPECT_EQ(hep1.get_base_species(), "he");
  EXPECT_EQ(hep1.get_str(), "he+");

  SpeciesParser hep2 = hep1.ionised();
  EXPECT_EQ(hep2.get_charge(), 2);
  EXPECT_EQ(hep2.get_base_species(), "he");
  EXPECT_EQ(hep2.get_str(), "he+2");
}

/// @brief Check multiple recombinations
TEST_F(SpeciesParserTest, MultipleRecombinations) {
  SpeciesParser hep2("he+2");
  EXPECT_EQ(hep2.get_charge(), 2);

  SpeciesParser hep1 = hep2.recombined();
  EXPECT_EQ(hep1.get_charge(), 1);
  EXPECT_EQ(hep1.get_base_species(), "he");
  EXPECT_EQ(hep1.get_str(), "he+");

  SpeciesParser he0 = hep1.recombined();
  EXPECT_EQ(he0.get_charge(), 0);
  EXPECT_EQ(he0.get_base_species(), "he");
  EXPECT_EQ(he0.get_str(), "he");
}

/// @brief Don't allow ionisation/recombination of electrons
TEST_F(SpeciesParserTest, NoIonisationRecombinationOfElectrons) {
  SpeciesParser electron("e");
  EXPECT_THROW(electron.ionised(), BoutException);
  EXPECT_THROW(electron.recombined(), BoutException);
}

/// @brief Check that element long names are set correctly.
TEST_F(SpeciesParserTest, LongNames) {
  // Species with preset long names
  SpeciesParser parser1("e");
  ASSERT_EQ(parser1.long_name(), "electron");
  SpeciesParser parser2("H");
  ASSERT_EQ(parser2.long_name(), "hydrogen");
  SpeciesParser parser3("h");
  ASSERT_EQ(parser3.long_name(), "hydrogen");
  SpeciesParser parser4("D");
  ASSERT_EQ(parser4.long_name(), "deuterium");
  SpeciesParser parser5("t");
  ASSERT_EQ(parser5.long_name(), "tritium");

  // Species without preset long names should return the base species name
  SpeciesParser parser6("li");
  ASSERT_EQ(parser6.long_name(), "li");
}

/// Check special case of electrons
TEST_F(SpeciesParserTest, ParseElectron) {
  // "e" parsed as electron
  SpeciesParser parser1("e");
  EXPECT_EQ(parser1.get_charge(), -1);
  EXPECT_EQ(parser1.get_base_species(), "e");
  EXPECT_EQ(parser1.get_str(), "e");

  // "e-" parsed as electron
  SpeciesParser parser2("e-");
  EXPECT_EQ(parser2.get_charge(), -1);
  EXPECT_EQ(parser2.get_base_species(), "e");
  EXPECT_EQ(parser2.get_str(), "e");

  // Charges other than -1 throw for electrons
  ASSERT_THROW(SpeciesParser("e-2"), BoutException);
  ASSERT_THROW(SpeciesParser("e+"), BoutException);
}

/// @brief Check that some properties of molecular species are parsed correctly.
TEST_F(SpeciesParserTest, ParseMolecules) {

  // A neutral molecule
  SpeciesParser parser1("d2");
  EXPECT_EQ(parser1.get_charge(), 0);
  EXPECT_EQ(parser1.get_base_species(), "d2");
  EXPECT_EQ(parser1.get_str(), "d2");

  // A molecular ion
  SpeciesParser parser2("d2+");
  EXPECT_EQ(parser2.get_charge(), 1);
  EXPECT_EQ(parser2.get_base_species(), "d2");
  EXPECT_EQ(parser2.get_str(), "d2+");

  // Another molecular ion
  SpeciesParser parser3("he3+");
  EXPECT_EQ(parser3.get_charge(), 1);
  EXPECT_EQ(parser3.get_base_species(), "he3");
  EXPECT_EQ(parser3.get_str(), "he3+");
}

/// @brief Check parsing of various negatively charged species.
TEST_F(SpeciesParserTest, ParseNegativelyCharged) {

  SpeciesParser parser1("He-1");
  EXPECT_EQ(parser1.get_base_species(), "he");
  EXPECT_EQ(parser1.get_charge(), -1);
  EXPECT_EQ(parser1.get_str(), "he-");

  SpeciesParser parser2("li-2");
  EXPECT_EQ(parser2.get_base_species(), "li");
  EXPECT_EQ(parser2.get_charge(), -2);
  EXPECT_EQ(parser2.get_str(), "li-2");

  SpeciesParser parser3("ne-5");
  EXPECT_EQ(parser3.get_base_species(), "ne");
  EXPECT_EQ(parser3.get_charge(), -5);
  EXPECT_EQ(parser3.get_str(), "ne-5");
}

/// @brief Check parsing of a neutral species with a prefix > 1
TEST_F(SpeciesParserTest, ParseNeutral) {
  SpeciesParser parser("2he");
  EXPECT_EQ(parser.get_base_species(), "he");
  EXPECT_EQ(parser.get_charge(), 0);
  // Species string should not include the prefix number and should be lower case
  EXPECT_EQ(parser.get_str(), "he");
}

/// @brief Extract every component in a prefixed molecular ion.
TEST_F(SpeciesParserTest, ParsePrefixedMolecularIon) {
  SpeciesParser parser("12He10-2");
  EXPECT_EQ(parser.get_element(), "he");
  EXPECT_EQ(parser.get_base_species(), "he10");
  EXPECT_TRUE(parser.is_molecule());
  EXPECT_EQ(parser.get_charge(), -2);
  EXPECT_EQ(parser.get_str(), "he10-2");
}

/// @brief Check parsing of various positively charged species.
TEST_F(SpeciesParserTest, ParsePositivelyCharged) {
  SpeciesParser parser1("h+");
  EXPECT_EQ(parser1.get_base_species(), "h");
  EXPECT_EQ(parser1.get_charge(), 1);
  EXPECT_EQ(parser1.get_str(), "h+");

  SpeciesParser parser2("He+2");
  EXPECT_EQ(parser2.get_base_species(), "he");
  EXPECT_EQ(parser2.get_charge(), 2);
  EXPECT_EQ(parser2.get_str(), "he+2");

  SpeciesParser parser3("ne+8");
  EXPECT_EQ(parser3.get_base_species(), "ne");
  EXPECT_EQ(parser3.get_charge(), 8);
  EXPECT_EQ(parser3.get_str(), "ne+8");
}

/// @brief Check recombination of a highly-ionised species.
TEST_F(SpeciesParserTest, RecombineHighlyIonised) {
  SpeciesParser highly("ne+9");
  SpeciesParser next = highly.recombined();
  EXPECT_EQ(next.get_base_species(), "ne");
  EXPECT_EQ(next.get_charge(), 8);
  EXPECT_EQ(next.get_str(), "ne+8");
}

/// @brief Check recombination of a singly-ionised species.
TEST_F(SpeciesParserTest, RecombineSinglyIonised) {
  SpeciesParser singly("h+");
  SpeciesParser neutral = singly.recombined();
  EXPECT_EQ(neutral.get_base_species(), "h");
  EXPECT_EQ(neutral.get_charge(), 0);
  EXPECT_EQ(neutral.get_str(), "h");
}
