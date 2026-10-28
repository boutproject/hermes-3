#include <bout/utils.hxx>
#include <gtest/gtest.h>

#include "species_parser.hxx"

/// Test fixture
class SpeciesParserTest : public ::testing::Test {};

/// @brief Check parsing of a neutral species with a prefix > 1
TEST_F(SpeciesParserTest, ParseNeutral) {
  SpeciesParser parser("2he");
  EXPECT_EQ(parser.get_base_species(), "he");
  EXPECT_EQ(parser.get_charge(), 0);
  // Species string should not include the prefix number and should be lower case
  EXPECT_EQ(parser.get_str(), "he");
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

/// @brief Check that element long names are parsed correctly.
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

/// @brief Check ionisation of a highly-ionised species.
TEST_F(SpeciesParserTest, IoniseHighlyIonised) {
  SpeciesParser highly("ne+8");
  SpeciesParser next = highly.ionised();
  EXPECT_EQ(next.get_base_species(), "ne");
  EXPECT_EQ(next.get_charge(), 9);
  EXPECT_EQ(next.get_str(), "ne+9");
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

/// @brief Check recombination of a singly-ionised species.
TEST_F(SpeciesParserTest, RecombineSinglyIonised) {
  SpeciesParser singly("h+");
  SpeciesParser neutral = singly.recombined();
  EXPECT_EQ(neutral.get_base_species(), "h");
  EXPECT_EQ(neutral.get_charge(), 0);
  EXPECT_EQ(neutral.get_str(), "h");
}

/// @brief Check recombination of a highly-ionised species.
TEST_F(SpeciesParserTest, RecombineHighlyIonised) {
  SpeciesParser highly("ne+9");
  SpeciesParser next = highly.recombined();
  EXPECT_EQ(next.get_base_species(), "ne");
  EXPECT_EQ(next.get_charge(), 8);
  EXPECT_EQ(next.get_str(), "ne+8");
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

//=============================== Error handling ==============================

/// @brief Don't allow ionisation/recombination of electrons
TEST_F(SpeciesParserTest, NoIonisationRecombinationOfElectrons) {
  SpeciesParser electron("e");
  EXPECT_THROW(electron.ionised(), BoutException);
  EXPECT_THROW(electron.recombined(), BoutException);
}

/// @brief Check that constructor throws for a variety of invalid strings
TEST_F(SpeciesParserTest, InvalidStrings) {
  ASSERT_THROW(SpeciesParser(""), BoutException);
  ASSERT_THROW(SpeciesParser("123"), BoutException);
  ASSERT_THROW(SpeciesParser("+"), BoutException);
  ASSERT_THROW(SpeciesParser("+2"), BoutException);
  ASSERT_THROW(SpeciesParser("-h"), BoutException);
  ASSERT_THROW(SpeciesParser("%&?!"), BoutException);
  ASSERT_THROW(SpeciesParser("he +"), BoutException);
  ASSERT_THROW(SpeciesParser("H -"), BoutException);
}
