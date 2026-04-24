#include "test_adas_reactions.hxx"
#include "test_cx_reactions.hxx"
#include "test_izn_rec_reactions.hxx"
#include "test_molecule_reactions.hxx"

namespace hermes {

// ======================= General reaction class tests =======================
/// @brief Test parsing of various input optionsReactionBase constructor should throw if
/// the reaction type string is not valid
TEST(ReactionTest, InputOptions) {
  const std::string comp_name = "test";

  // Base input with two reaction strings
  Options base_input{
      {comp_name, {{"type", "(h + h+ -> h+ + h, d + d+ -> d+ + d)"}}},
      {"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}}};

  // Setting one or two data_srcs should work
  Options valid_input1 = base_input.copy();
  valid_input1[comp_name]["data_srcs"] = "(Amjuel)";
  Options valid_input2 = base_input.copy();
  valid_input2[comp_name]["data_srcs"] = "(Amjuel,Amjuel)";
  ReactionBase::reset_instance_counter();
  ASSERT_NO_THROW(CXReaction(comp_name, valid_input1));
  ReactionBase::reset_instance_counter();
  ASSERT_NO_THROW(CXReaction(comp_name, valid_input2));

  // Setting num_data_srcs != (1 || num_reactions) should throw
  Options invalid_input1 = base_input.copy();
  invalid_input1[comp_name]["data_srcs"] = "(Amjuel,Amjuel,Amjuel)";
  ReactionBase::reset_instance_counter();
  ASSERT_THROW(CXReaction(comp_name, invalid_input1), BoutException);

  // Setting one or two data IDs should work
  Options valid_input3 = base_input.copy();
  valid_input3[comp_name]["data_ids"] = "H.2_3.1.8";
  Options valid_input4 = base_input.copy();
  valid_input4[comp_name]["data_ids"] = "H.2_3.1.8,H.2_3.1.8";
  ReactionBase::reset_instance_counter();
  ASSERT_NO_THROW(CXReaction(comp_name, valid_input3));
  ReactionBase::reset_instance_counter();
  ASSERT_NO_THROW(CXReaction(comp_name, valid_input4));

  // Setting num_data_ids != (1 || num_reactions) should throw
  Options invalid_input2 = base_input.copy();
  invalid_input2[comp_name]["data_ids"] = "H.2_3.1.8,H.2_3.1.8,H.2_3.1.8";
  ReactionBase::reset_instance_counter();
  ASSERT_THROW(CXReaction(comp_name, invalid_input2), BoutException);
}

// ========================== CX reaction class tests =========================

/// @brief CXReaction constructor should throw for strings that aren't valid CX
/// reactions
TEST(CXReactionTest, InvalidReactionStrings) {
  Options base_options{
      {"test", {{"data_ids", "H.2_3.1.8"}}},
      {"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}}};

  // Invalid CX reaction strings
  std::string too_few_reactants = "h -> h+ + h";
  std::string too_many_reactants = "h + d+ + t -> d + h+";
  std::string too_few_products = "h + h+ -> h+";
  std::string too_many_products = "h + h+ -> h+ + h + d";
  std::string invalid_cx1 = "h + d+ -> d + t+";
  std::string invalid_cx2 = "h + d+ -> t + h+";

  // Test that constructor throws for each invalid reaction string
  for (const auto& invalid_reaction_str :
       {too_few_reactants, too_many_reactants, too_few_products, too_many_products,
        invalid_cx1, invalid_cx2}) {
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

// =================== ElasticCollision reaction class tests ==================
/// @brief Test validity of various reaction strings for elastic collisions
TEST(ElasticCollisionTest, ReactionStrings) {
  Options base_options{
      {"test", {{"data_ids", "H.2_0.3T"}}},
      {"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}}};

  // Swapped order should work
  std::string valid_reaction_str = "h2 + d+ -> d+ + h2 ";
  ReactionBase::reset_instance_counter();
  Options options = base_options.copy();
  options["test"]["type"] = valid_reaction_str;
  ASSERT_NO_THROW(ElasticCollision("test", options));

  // Invalid elastic collision reaction strings
  std::string too_few_species = "h -> h";
  std::string too_many_species = "h + d+ + t -> h + d+ + t";
  std::string reactants_products_mismatch = "h2 + d+ -> h2 +d";

  // Test that constructor throws for each invalid reaction string
  for (const auto& invalid_reaction_str :
       {too_few_species, too_many_species, reactants_products_mismatch}) {
    ReactionBase::reset_instance_counter();
    Options options = base_options.copy();
    options["test"]["type"] = invalid_reaction_str;
    ASSERT_THROW(ElasticCollision("test", options), BoutException);
  }
}

// ===================== Reaction source regression tests =====================

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

/* Various dissociation reactions involving neutral/ionic molecules of H
 * isotopes(non-exhaustive)
 */
TEST_F(H2DissTest, SourcesRegression) { sources_regression_test(); }
TEST_F(T2pDissExcTest, SourcesRegression) { sources_regression_test(); }
TEST_F(D2DissIznTest, SourcesRegression) { sources_regression_test(); }
TEST_F(H2NonDissIznTest, SourcesRegression) { sources_regression_test(); }
TEST_F(T2pDissRecTest, SourcesRegression) { sources_regression_test(); }

// CX involving molecules of H isotopes (non-exhaustive)
TEST_F(D2DpCXTest, SourcesRegression) { sources_regression_test(); }

// Elastic collision between H2 and H+
TEST_F(H2HpElasticCollisionTest, SourcesRegression) { sources_regression_test(); }

// He ionization
TEST_F(HeIzn01Test, SourcesRegression) { sources_regression_test(); }

// He recombination
TEST_F(HeRec10Test, SourcesRegression) { sources_regression_test(); }

} // namespace hermes

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
