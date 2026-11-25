#include "test_adas_reactions.hxx"
#include "test_cx_reactions.hxx"
#include "test_izn_rec_reactions.hxx"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <vector>

namespace hermes {

//======================== General reaction class tests =======================
/// @brief Test parsing of various input optionsReactionBase constructor should throw if
/// the reaction type string is not

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
  ASSERT_NO_THROW(ConcreteCXReaction(comp_name, valid_input1));
  ReactionBase::reset_instance_counter();
  ASSERT_NO_THROW(ConcreteCXReaction(comp_name, valid_input2));

  // Setting num_data_srcs != (1 || num_reactions) should throw
  Options invalid_input1 = base_input.copy();
  invalid_input1[comp_name]["data_srcs"] = "(Amjuel,Amjuel,Amjuel)";
  ReactionBase::reset_instance_counter();
  ASSERT_THROW(ConcreteCXReaction(comp_name, invalid_input1), BoutException);

  // Setting one or two data IDs should work
  Options valid_input3 = base_input.copy();
  valid_input3[comp_name]["data_ids"] = "H.2_3.1.8";
  Options valid_input4 = base_input.copy();
  valid_input4[comp_name]["data_ids"] = "H.2_3.1.8,H.2_3.1.8";
  ReactionBase::reset_instance_counter();
  ASSERT_NO_THROW(ConcreteCXReaction(comp_name, valid_input3));
  ReactionBase::reset_instance_counter();
  ASSERT_NO_THROW(ConcreteCXReaction(comp_name, valid_input4));

  // Setting num_data_ids != (1 || num_reactions) should throw
  Options invalid_input2 = base_input.copy();
  invalid_input2[comp_name]["data_ids"] = "H.2_3.1.8,H.2_3.1.8,H.2_3.1.8";
  ReactionBase::reset_instance_counter();
  ASSERT_THROW(ConcreteCXReaction(comp_name, invalid_input2), BoutException);
}

//======================== CX reaction class tests =======================

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
    ASSERT_THROW(ConcreteCXReaction("test", options), BoutException);
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
    ASSERT_NO_THROW(ConcreteCXReaction("test", options));
  }
}

//====================== Reaction source regression tests =====================

// H isotopes ionization
TEST_F(HIznTest, SourcesRegression) {
  sources_regression_test<bout::build::use_metric_3d>();
}
TEST_F(DIznTest, SourcesRegression) {
  sources_regression_test<bout::build::use_metric_3d>();
}
TEST_F(TIznTest, SourcesRegression) {
  sources_regression_test<bout::build::use_metric_3d>();
}

// H isotopes recombination
TEST_F(HRecTest, SourcesRegression) {
  sources_regression_test<bout::build::use_metric_3d>();
}
TEST_F(DRecTest, SourcesRegression) {
  sources_regression_test<bout::build::use_metric_3d>();
}
TEST_F(TRecTest, SourcesRegression) {
  sources_regression_test<bout::build::use_metric_3d>();
}

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
TEST_F(HeIzn01Test, SourcesRegression) {
  sources_regression_test<bout::build::use_metric_3d>();
}

// He recombination
TEST_F(HeRec10Test, SourcesRegression) {
  sources_regression_test<bout::build::use_metric_3d>();
}

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

// Ne ionisation: radiated power in SI units at a node of the PLT table (W m^3)
TEST_F(ADASNeIznTest, RadiatedPowerSI) {
  // Non-trivial normalisations, so that the unit conversion is tested
  const BoutReal Tnorm = 5.0;  // eV
  const BoutReal Nnorm = 1e18; // m^-3
  const BoutReal tnorm = 1e-6; // s

  Options state = generate_state();
  state["units"]["eV"].force(Tnorm);
  state["units"]["inv_meters_cubed"].force(Nnorm);
  state["units"]["seconds"].force(tnorm);

  // PLT table node closest to Te = 10 eV, ne = 1e19 m^-3
  OpenADASRateCoefficient plt(get_json_db_dir(state) / "plt96_ne.json", 0);
  const auto closest = [](const std::vector<BoutReal>& v, BoutReal x) {
    return static_cast<std::size_t>(std::distance(
        v.begin(), std::min_element(v.begin(), v.end(), [x](BoutReal a, BoutReal b) {
          return std::abs(a - x) < std::abs(b - x);
        })));
  };
  const std::size_t iT = closest(plt.log_temperature, 1.0);
  const std::size_t in = closest(plt.log_density, 19.0);
  const BoutReal Te_SI = std::pow(10., plt.log_temperature[iT]);  // eV
  const BoutReal ne_SI = std::pow(10., plt.log_density[in]);      // m^-3
  const BoutReal n1_SI = 1e17;                                    // m^-3
  const BoutReal plt_node = std::pow(10., plt.log_coeff[iT][in]); // W m^3

  state["species"]["e"]["density"].force(Field3D(ne_SI / Nnorm));
  state["species"]["e"]["temperature"].force(Field3D(Te_SI / Tnorm));
  state["species"]["ne"]["density"].force(Field3D(n1_SI / Nnorm));

  hermes::ReactionBase::reset_instance_counter();
  ADASNeonIonisation<0> component("test" + lbl, state, nullptr);
  component.transform(state);

  const Field3D rate = get<Field3D>(state["species"]["ne+"]["density_source"]);
  const Field3D energy = get<Field3D>(state["species"]["e"]["energy_source"]);
  const BoutReal rate_to_SI = Nnorm / tnorm;                   // m^-3 s^-1
  const BoutReal power_to_SI = SI::qe * Tnorm * Nnorm / tnorm; // W m^-3

  // Electron energy loss = radiated power + ionisation energy * ionisation rate
  const BoutReal expected = ne_SI * n1_SI * plt_node; // W m^-3
  BOUT_FOR_SERIAL(i, energy.getRegion("RGN_NOBNDRY")) {
    const BoutReal radiated = -energy[i] * power_to_SI
                              - SI::qe * neon_ionisation_energy[0] * rate[i] * rate_to_SI;
    ASSERT_NEAR(radiated, expected, 1e-8 * expected);
  }
}

// Ne ions CX  (non-exhaustive)
TEST_F(ADASNepHCXTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASNep5DCXTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASNep9TCXTest, SourcesRegression) { sources_regression_test(); }
