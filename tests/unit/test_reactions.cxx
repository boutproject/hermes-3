#include "test_adas_reactions.hxx"
#include "test_cx_reactions.hxx"
#include "test_izn_rec_reactions.hxx"

//======================== General reaction class tests =======================

/// @brief CXReaction constructor should throw for strings that aren't valid CX reactions
TEST(CXReactionTest, InvalidReactionStrings) {
  Options base_options{
      {"test", {{"data_ids", "H.2_3.1.8"}}},
      {"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}}};

  // Invalid CX reaction strings
  std::string too_few_reactants = "h -> h+ + h";
  std::string too_many_reactants = "h + d+ + t -> d + h+";
  std::string too_few_products = "h + h+ -> h+";
  std::string too_many_products = "h + h+ -> h+ + h + d";
  std::string not_cx = "h + e -> h+ + 2e";

  // Test that constructor throws for each invalid reaction string
  for (const auto& invalid_reaction_str : {too_few_reactants, too_many_reactants,
                                           too_few_products, too_many_products, not_cx}) {
    ReactionBase::reset_instance_counter();
    Options options = base_options.copy();
    options["test"]["type"] = invalid_reaction_str;
    ASSERT_THROW(CXReaction("test", options), BoutException);
  }
}

/// @brief CXReaction should accept strings with reactants and products in either order
TEST(CXReactionTest, OrderIndependentReactionStrs) {
  Options base_options{
      {"test", {{"data_ids", "H.2_3.1.8"}}},
      {"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}}};

  std::string DT_cx_order1 = "d + t+ -> t + d+";
  std::string DT_cx_order2 = "d + t+ -> d+ + t";
  std::string DT_cx_order3 = "t+ + d -> t + d+";
  std::string DT_cx_order4 = "t+ + d -> d+ + t";

  // Test that all of the strings allow construction without throwing
  for (const auto& valid_reaction_str :
       {DT_cx_order1, DT_cx_order2, DT_cx_order3, DT_cx_order4}) {
    ReactionBase::reset_instance_counter();
    Options options = base_options.copy();
    options["test"]["type"] = valid_reaction_str;
    ASSERT_NO_THROW(CXReaction("test", options));
  }
}

//====================== Reaction source regression tests =====================

// H isotopes ionization
TEST_F(HIznTest, SourcesRegression) { sources_regression_test(); }
TEST_F(DIznTest, SourcesRegression) { sources_regression_test(); }
TEST_F(TIznTest, SourcesRegression) { sources_regression_test(); }

// H isotopes recombination
TEST_F(HRecTest, SourcesRegression) { sources_regression_test(); }
TEST_F(DRecTest, SourcesRegression) { sources_regression_test(); }
TEST_F(TRecTest, SourcesRegression) { sources_regression_test(); }

// H isotopes CX (non-exhaustive)
TEST_F(HHpCXTest, SourcesRegression) { sources_regression_test(); }
TEST_F(DDpCXTest, SourcesRegression) { sources_regression_test(); }
TEST_F(TTpCXTest, SourcesRegression) { sources_regression_test(); }
TEST_F(HDpCXTest, SourcesRegression) { sources_regression_test(); }
TEST_F(THpCXTest, SourcesRegression) { sources_regression_test(); }
TEST_F(DTpCXTest, SourcesRegression) { sources_regression_test(); }

// H/H+ CX with neutral momentum gain turned off
TEST_F(HHpCXTest_noNeutralMomGain, SourcesRegression) { sources_regression_test(); }

// He ionization
TEST_F(HeIzn01Test, SourcesRegression) { sources_regression_test(); }

// He recombination
TEST_F(HeRec10Test, SourcesRegression) { sources_regression_test(); }

// C ions ionization (non-exhaustive)
TEST_F(ADASCIznTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASCpIznTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASC5pIznTest, SourcesRegression) { sources_regression_test(); }

// C ions recombination (non-exhaustive)
TEST_F(ADASCRecTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASCpRecTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASC5pRecTest, SourcesRegression) { sources_regression_test(); }

// C ions CX  (non-exhaustive)
TEST_F(ADASCpHCXTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASCp3DCXTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASCp5TCXTest, SourcesRegression) { sources_regression_test(); }

// Li ions ionization
TEST_F(ADASLiIznTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASLipIznTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASLi2pIznTest, SourcesRegression) { sources_regression_test(); }

// Li ions recombination
TEST_F(ADASLiRecTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASLipRecTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASLi2pRecTest, SourcesRegression) { sources_regression_test(); }

// Li ions CX (non-exhaustive)
TEST_F(ADASLipHCXTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASLip2DCXTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASLip3TCXTest, SourcesRegression) { sources_regression_test(); }

// Ne ions ionization (non-exhaustive)
TEST_F(ADASNeIznTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASNepIznTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASNe9pIznTest, SourcesRegression) { sources_regression_test(); }

// Ne ions recombination (non-exhaustive)
TEST_F(ADASNeRecTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASNepRecTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASNe9pRecTest, SourcesRegression) { sources_regression_test(); }

// Ne ions CX  (non-exhaustive)
TEST_F(ADASNepHCXTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASNep5DCXTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASNep9TCXTest, SourcesRegression) { sources_regression_test(); }
