#include <bout/utils.hxx>
#include <gtest/gtest.h>

#include "species_parser.hxx"

/// Test fixture
class SpeciesParserTest : public ::testing::Test {};

/// @brief Check parsing of a neutral species with a prefix > 1
TEST_F(SpeciesParserTest, ParseNeutral) {
  SpeciesParser parser("2he");
  EXPECT_EQ(parser.get_element(), "he");
  EXPECT_EQ(parser.get_charge(), 0);
  // Species string should not include the prefix number and should be lower case
  EXPECT_EQ(parser.get_str(), "he");
}

/// @brief Check parsing of a singly-ionised species.
TEST_F(SpeciesParserTest, ParseSinglyIonised) {
  SpeciesParser parser("h+");
  EXPECT_EQ(parser.get_element(), "h");
  EXPECT_EQ(parser.get_charge(), 1);
  EXPECT_EQ(parser.get_str(), "h+");
}

/// @brief Check parsing of a doubly-ionised species.
TEST_F(SpeciesParserTest, ParseDoublyIonised) {
  SpeciesParser parser("He+2");
  EXPECT_EQ(parser.get_element(), "he");
  EXPECT_EQ(parser.get_charge(), 2);
  EXPECT_EQ(parser.get_str(), "he+2");
}

// Check special case of electrons
TEST_F(SpeciesParserTest, ParseElectron) {
  SpeciesParser parser1("e");
  EXPECT_EQ(parser1.get_charge(), -1);
  EXPECT_EQ(parser1.get_element(), "e");
  EXPECT_EQ(parser1.get_str(), "e");

  SpeciesParser parser2("e-");
  EXPECT_EQ(parser2.get_charge(), -1);
  EXPECT_EQ(parser2.get_element(), "e");
  EXPECT_EQ(parser2.get_str(), "e");
}

/// @brief Check ionisation of a neutral species.
TEST_F(SpeciesParserTest, IoniseNeutral) {
  SpeciesParser neutral("h");
  SpeciesParser ionised = neutral.ionised();
  EXPECT_EQ(ionised.get_element(), "h");
  EXPECT_EQ(ionised.get_charge(), 1);
  EXPECT_EQ(ionised.get_str(), "h+");
}

/// @brief Check ionisation of a singly-ionised species.
TEST_F(SpeciesParserTest, IoniseSinglyIonised) {
  SpeciesParser singly("li+");
  SpeciesParser doubly = singly.ionised();
  EXPECT_EQ(doubly.get_element(), "li");
  EXPECT_EQ(doubly.get_charge(), 2);
  EXPECT_EQ(doubly.get_str(), "li+2");
}

/// @brief Check ionisation of a highly-ionised species.
TEST_F(SpeciesParserTest, IoniseHighlyIonised) {
  SpeciesParser highly("ne+8");
  SpeciesParser next = highly.ionised();
  EXPECT_EQ(next.get_element(), "ne");
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
  EXPECT_EQ(hep1.get_element(), "he");
  EXPECT_EQ(hep1.get_str(), "he+");

  SpeciesParser hep2 = hep1.ionised();
  EXPECT_EQ(hep2.get_charge(), 2);
  EXPECT_EQ(hep2.get_element(), "he");
  EXPECT_EQ(hep2.get_str(), "he+2");
}

/// @brief Check recombination of a singly-ionised species.
TEST_F(SpeciesParserTest, RecombineSinglyIonised) {
  SpeciesParser singly("h+");
  SpeciesParser neutral = singly.recombined();
  EXPECT_EQ(neutral.get_element(), "h");
  EXPECT_EQ(neutral.get_charge(), 0);
  EXPECT_EQ(neutral.get_str(), "h");
}

/// @brief Check recombination of a highly-ionised species.
TEST_F(SpeciesParserTest, RecombineHighlyIonised) {
  SpeciesParser highly("ne+9");
  SpeciesParser next = highly.recombined();
  EXPECT_EQ(next.get_element(), "ne");
  EXPECT_EQ(next.get_charge(), 8);
  EXPECT_EQ(next.get_str(), "ne+8");
}

/// @brief Check multiple recombinations
TEST_F(SpeciesParserTest, MultipleRecombinations) {
  SpeciesParser hep2("he+2");
  EXPECT_EQ(hep2.get_charge(), 2);

  SpeciesParser hep1 = hep2.recombined();
  EXPECT_EQ(hep1.get_charge(), 1);
  EXPECT_EQ(hep1.get_element(), "he");
  EXPECT_EQ(hep1.get_str(), "he+");

  SpeciesParser he0 = hep1.recombined();
  EXPECT_EQ(he0.get_charge(), 0);
  EXPECT_EQ(he0.get_element(), "he");
  EXPECT_EQ(he0.get_str(), "he");
}

//=============================== Error handling ==============================

/// @brief Don't allow ionisation/recombination of electrons
TEST_F(SpeciesParserTest, NoIonisationRecombinationOfElectrons) {
  SpeciesParser electron("e");
  EXPECT_THROW(electron.ionised(), BoutException);
  EXPECT_THROW(electron.recombined(), BoutException);
}

/// @brief Check that constructor throws for invalid strings
TEST_F(SpeciesParserTest, InvalidStrings) {
  EXPECT_THROW(SpeciesParser(""), BoutException);
  EXPECT_THROW(SpeciesParser("123"), BoutException);
  EXPECT_THROW(SpeciesParser("+"), BoutException);
  EXPECT_THROW(SpeciesParser("+2"), BoutException);
}