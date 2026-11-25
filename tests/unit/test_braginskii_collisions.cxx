#include <cmath>
#include <string>

#include <bout/bout_types.hxx>
#include <bout/constants.hxx>
#include <bout/region.hxx>
#include "gtest/gtest.h"

#include "../../include/braginskii_collisions.hxx"
#include "../../include/component.hxx"
#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using BraginskiiCollisionsTest = FakeMeshFixture;

TEST_F(BraginskiiCollisionsTest, OnlyElectrons) {
  Options options;

  options["units"]["eV"] = 1.0;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["inv_meters_cubed"] = 1.0;

  BraginskiiCollisions component("test", options, nullptr);

  Options state;
  state["species"]["e"]["density"] = 1e19;
  state["species"]["e"]["temperature"] = 10.;

  component.declareAllSpecies(SpeciesInformation({"e"}, {}, {}, {}));
  component.transform(state);

  ASSERT_TRUE(state["species"]["e"].isSet("collision_frequency"));
}

TEST_F(BraginskiiCollisionsTest, OneOrTwoSpeciesCharged) {
  Options options;

  options["units"]["eV"] = 1.0;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["inv_meters_cubed"] = 1.0;

  BraginskiiCollisions component("test", options, nullptr);

  Options state1;
  state1["species"]["s1"]["density"] = 1e19;
  state1["species"]["s1"]["temperature"] = 10;
  state1["species"]["s1"]["charge"] = 1;
  state1["species"]["s1"]["AA"] = 2;

  Options state2;
  state2["species"]["s1"]["density"] = 5e18; // Half density
  state2["species"]["s1"]["temperature"] = 10;
  state2["species"]["s1"]["charge"] = 1;
  state2["species"]["s1"]["AA"] = 2;

  state2["species"]["s2"] = state2["species"]["s1"].copy();

  // Run calculations
  component.declareAllSpecies(SpeciesInformation({}, {}, {"s1", "s2"}, {}));
  component.transform(state1);
  component.transform(state2);

  Field3D nu1 = get<Field3D>(state1["species"]["s1"]["collision_frequency"]);
  Field3D nu21 = get<Field3D>(state2["species"]["s1"]["collision_frequency"]);
  Field3D nu22 = get<Field3D>(state2["species"]["s2"]["collision_frequency"]);

  BOUT_FOR_SERIAL(i, nu1.getRegion("RGN_ALL")) {
    // Collision rates for species1 and species2 should be equal
    ASSERT_DOUBLE_EQ(nu21[i], nu22[i]);

    // Whether one or two species, the collision rates should be similar
    // Note: Not exactly the same, because coulomb logarithm is slightly different
    ASSERT_TRUE(abs(nu1[i] - nu21[i]) / (nu1[i] + nu21[i]) < 0.05);
  }
}

TEST_F(BraginskiiCollisionsTest, TnormDependence) {
  // Calculate rates with normalisation factors 1
  Options options{
      {"units",
       {{"eV", 1.0}, {"meters", 1.0}, {"seconds", 1.0}, {"inv_meters_cubed", 1.0}}},
      {"test",
       {{"electron_neutral", true},
        {"ion_neutral", true},
        {"electron_electron", true},
        {"ion_ion", true},
        {"neutral_neutral", true}}}};

  BraginskiiCollisions component("test", options, nullptr);

  Options state{
      {"species",
       {{"e",
         {{"density", 1e19}, {"temperature", 10}, {"charge", -1}, {"AA", 1. / 1836}}},
        {"d+", {{"density", 2e19}, {"temperature", 20}, {"charge", 1}, {"AA", 2}}},
        {"d", {{"density", 1e18}, {"temperature", 3}, {"AA", 2}}}}}};

  component.declareAllSpecies(SpeciesInformation({"e"}, {"d"}, {"d+"}, {}));
  component.transform(state);

  ASSERT_TRUE(state["species"]["e"].isSet("collision_frequency"));
  ASSERT_TRUE(state["species"]["d"].isSet("collision_frequency"));
  ASSERT_TRUE(state["species"]["d+"].isSet("collision_frequency"));

  // Collision frequencies should be positive, non-zero
  ASSERT_GT(get<Field3D>(state["species"]["e"]["collision_frequency"])(0, 0, 0), 0.0);
  ASSERT_GT(get<Field3D>(state["species"]["d"]["collision_frequency"])(0, 0, 0), 0.0);
  ASSERT_GT(get<Field3D>(state["species"]["d+"]["collision_frequency"])(0, 0, 0), 0.0);

  // Re-calculate with Tnorm != 1
  // To keep frequency normalisation fixed, rho_s0 scales like sqrt(Tnorm)
  const BoutReal Tnorm = 100;
  Options options2{{"units",
                    {{"eV", Tnorm},
                     {"meters", sqrt(Tnorm)},
                     {"seconds", 1.0},
                     {"inv_meters_cubed", 1.0}}},
                   {"test",
                    {{"electron_neutral", true},
                     {"ion_neutral", true},
                     {"electron_electron", true},
                     {"ion_ion", true},
                     {"neutral_neutral", true}}}};

  BraginskiiCollisions component2("test", options2, nullptr);

  Options state2{
      {"species",
       {{"e",
         {{"density", 1e19},
          {"temperature", 10 / Tnorm},
          {"charge", -1},
          {"AA", 1. / 1836}}},
        {"d+",
         {{"density", 2e19}, {"temperature", 20 / Tnorm}, {"charge", 1}, {"AA", 2}}},
        {"d", {{"density", 1e18}, {"temperature", 3 / Tnorm}, {"AA", 2}}}}}};

  component2.declareAllSpecies(SpeciesInformation({"e"}, {"d"}, {"d+"}, {}));
  component2.transform(state2);

  // Normalised frequencies should be unchanged
  ASSERT_DOUBLE_EQ(get<Field3D>(state["species"]["e"]["collision_frequency"])(0, 0, 0),
                   get<Field3D>(state2["species"]["e"]["collision_frequency"])(0, 0, 0));
  ASSERT_DOUBLE_EQ(get<Field3D>(state["species"]["d"]["collision_frequency"])(0, 0, 0),
                   get<Field3D>(state2["species"]["d"]["collision_frequency"])(0, 0, 0));
  ASSERT_DOUBLE_EQ(get<Field3D>(state["species"]["d+"]["collision_frequency"])(0, 0, 0),
                   get<Field3D>(state2["species"]["d+"]["collision_frequency"])(0, 0, 0));
}

TEST_F(BraginskiiCollisionsTest, IonNeutral) {
  Options options{
      {"units",
       {{"eV", 1.0}, {"meters", 1.0}, {"seconds", 1.0}, {"inv_meters_cubed", 1.0}}},
      {"test",
       {{"electron_neutral", false},
        {"ion_neutral", true},
        {"electron_electron", false},
        {"ion_ion", false},
        {"neutral_neutral", false}}}};

  BraginskiiCollisions component("test", options, nullptr);

  // Put ion in between neutral species to test both orders (ion-neutral and neutral-ion)
  Options state{
      {"species",
       {{"d", {{"density", 1e18}, {"temperature", 3}, {"AA", 2}}},
        {"d+", {{"density", 2e19}, {"temperature", 20}, {"charge", 1}, {"AA", 2}}},
        {"d2", {{"density", 1e17}, {"temperature", 1}, {"AA", 4}}}}}};

  component.declareAllSpecies(SpeciesInformation({}, {"d", "d2"}, {"d+"}, {}));
  component.transform(state);

  ASSERT_TRUE(state["species"]["d"].isSet("collision_frequency"));
  ASSERT_TRUE(state["species"]["d"]["collision_frequencies"].isSet("d_d+_coll"));

  ASSERT_TRUE(state["species"]["d2"].isSet("collision_frequency"));
  ASSERT_TRUE(state["species"]["d2"]["collision_frequencies"].isSet("d2_d+_coll"));

  ASSERT_TRUE(state["species"]["d+"].isSet("collision_frequency"));
  ASSERT_TRUE(state["species"]["d+"]["collision_frequencies"].isSet("d+_d_coll"));
  ASSERT_TRUE(state["species"]["d+"]["collision_frequencies"].isSet("d+_d2_coll"));

  Options options2{
      {"units",
       {{"eV", 1.0}, {"meters", 1.0}, {"seconds", 1.0}, {"inv_meters_cubed", 1.0}}},
      {"test",
       {{"electron_neutral", false},
        {"ion_neutral", false}, // Turn off ion_neutral collisions
        {"electron_electron", false},
        {"ion_ion", false},
        {"neutral_neutral", false}}}};

  BraginskiiCollisions component2("test", options2, nullptr);

  Options state2{
      {"species",
       {{"d", {{"density", 1e18}, {"temperature", 3}, {"AA", 2}}},
        {"d+", {{"density", 2e19}, {"temperature", 20}, {"charge", 1}, {"AA", 2}}},
        {"d2", {{"density", 1e17}, {"temperature", 1}, {"AA", 4}}}}}};

  component2.declareAllSpecies(SpeciesInformation({}, {"d", "d2"}, {"d+"}, {}));
  component2.transform(state2);

  // All ion-neutral collision frequencies should be missing
  ASSERT_FALSE(state2["species"]["d"].isSet("collision_frequency"));
  ASSERT_FALSE(state2["species"]["d"]["collision_frequencies"].isSet("d_d+_coll"));

  ASSERT_FALSE(state2["species"]["d2"].isSet("collision_frequency"));
  ASSERT_FALSE(state2["species"]["d2"]["collision_frequencies"].isSet("d2_d+_coll"));

  ASSERT_FALSE(state2["species"]["d+"].isSet("collision_frequency"));
  ASSERT_FALSE(state2["species"]["d+"]["collision_frequencies"].isSet("d+_d_coll"));
  ASSERT_FALSE(state2["species"]["d+"]["collision_frequencies"].isSet("d+_d2_coll"));
}

// Options ion_neutral_sigma and ion_neutral_vrel_reduced.
// The ion-neutral rate nu = v_rel * n2 * sigma * rho_s0 (rho_s0 = 1 here) is compared
// with the hand formula in both branches of the calculation: d_d+_coll (neutral species
// first) and d+_d2_coll (charged species first), with the default options and with
// sigma = 1e-18 m^2 and the mean relative speed of two Maxwellians.
TEST_F(BraginskiiCollisionsTest, IonNeutralSigmaAndReducedVrelBothBranches) {
  auto make_options = [](bool set_new) {
    Options options{
        {"units",
         {{"eV", 1.0}, {"meters", 1.0}, {"seconds", 1.0}, {"inv_meters_cubed", 1.0}}},
        {"test",
         {{"electron_neutral", false},
          {"ion_neutral", true},
          {"electron_electron", false},
          {"ion_ion", false},
          {"neutral_neutral", false}}}};
    if (set_new) {
      options["test"]["ion_neutral_sigma"] = 1e-18;
      options["test"]["ion_neutral_vrel_reduced"] = true;
    }
    return options;
  };
  auto make_state = []() {
    return Options{
        {"species",
         {{"d", {{"density", 1e18}, {"temperature", 3}, {"AA", 2}}},
          {"d+", {{"density", 2e19}, {"temperature", 20}, {"charge", 1}, {"AA", 2}}},
          {"d2", {{"density", 1e17}, {"temperature", 1}, {"AA", 4}}}}}};
  };
  auto rate = [](Options& state, const std::string& species, const std::string& coll) {
    return get<Field3D>(state["species"][species]["collision_frequencies"][coll])(1, 1,
                                                                                  0);
  };

  // Default options: previous hard-coded rate
  Options opt_default = make_options(false);
  BraginskiiCollisions comp_default("test", opt_default, nullptr);
  Options state_default = make_state();
  comp_default.declareAllSpecies(SpeciesInformation({}, {"d", "d2"}, {"d+"}, {}));
  comp_default.transform(state_default);

  EXPECT_DOUBLE_EQ(rate(state_default, "d", "d_d+_coll"),
                   sqrt(3. / 2. + 20. / 2.) * 2e19 * 5e-19);
  EXPECT_DOUBLE_EQ(rate(state_default, "d+", "d+_d2_coll"),
                   sqrt(20. / 2. + 1. / 4.) * 1e17 * 5e-19);

  // sigma = 1e-18 and v_rel = sqrt(8/pi (T1/m1 + T2/m2))
  Options opt_new = make_options(true);
  BraginskiiCollisions comp_new("test", opt_new, nullptr);
  Options state_new = make_state();
  comp_new.declareAllSpecies(SpeciesInformation({}, {"d", "d2"}, {"d+"}, {}));
  comp_new.transform(state_new);

  EXPECT_DOUBLE_EQ(rate(state_new, "d", "d_d+_coll"),
                   sqrt((8. / PI) * (3. / 2. + 20. / 2.)) * 2e19 * 1e-18);
  EXPECT_DOUBLE_EQ(rate(state_new, "d+", "d+_d2_coll"),
                   sqrt((8. / PI) * (20. / 2. + 1. / 4.)) * 1e17 * 1e-18);

  // Different masses and temperatures (d+: T = 20 eV, AA = 2; d2: T = 1 eV, AA = 4),
  // value calculated by hand: sqrt(8/pi * (20/2 + 1/4)) * 1e17 * 1e-18 = 0.51089539700
  EXPECT_NEAR(rate(state_new, "d+", "d+_d2_coll"), 0.51089539700, 1e-10);

  // Equal masses (d, d+): the same as sqrt(8 T / (pi mu)) with T = (T1 + T2) / 2 and
  // 1/mu = 1/m1 + 1/m2
  EXPECT_DOUBLE_EQ(rate(state_new, "d", "d_d+_coll"),
                   sqrt((8. / PI) * 0.5 * (3. + 20.) * (1. / 2. + 1. / 2.)) * 2e19
                       * 1e-18);

  // The rates of the other species of each pair change by the same factor
  EXPECT_NEAR(rate(state_new, "d+", "d+_d_coll") / rate(state_default, "d+", "d+_d_coll"),
              rate(state_new, "d", "d_d+_coll") / rate(state_default, "d", "d_d+_coll"),
              1e-12);
  EXPECT_NEAR(
      rate(state_new, "d2", "d2_d+_coll") / rate(state_default, "d2", "d2_d+_coll"),
      rate(state_new, "d+", "d+_d2_coll") / rate(state_default, "d+", "d+_d2_coll"),
      1e-12);
}

// Equal masses and temperatures: sigma = 1e-18 and the mean relative speed increase
// the rate by 2 sqrt(8 / pi) in both directions
TEST_F(BraginskiiCollisionsTest, IonNeutralSigmaAndReducedVrelEqualMasses) {
  auto make_options = [](bool set_new) {
    Options options{
        {"units",
         {{"eV", 1.0}, {"meters", 1.0}, {"seconds", 1.0}, {"inv_meters_cubed", 1.0}}},
        {"test",
         {{"electron_neutral", false},
          {"ion_neutral", true},
          {"electron_electron", false},
          {"ion_ion", false},
          {"neutral_neutral", false}}}};
    if (set_new) {
      options["test"]["ion_neutral_sigma"] = 1e-18;
      options["test"]["ion_neutral_vrel_reduced"] = true;
    }
    return options;
  };
  auto make_state = []() {
    return Options{
        {"species",
         {{"ar", {{"density", 1e19}, {"temperature", 0.3}, {"AA", 40}}},
          {"ar+",
           {{"density", 1e18}, {"temperature", 0.3}, {"charge", 1}, {"AA", 40}}}}}};
  };

  Options opt_old = make_options(false);
  BraginskiiCollisions comp_old("test", opt_old, nullptr);
  Options state_old = make_state();
  comp_old.declareAllSpecies(SpeciesInformation({}, {"ar"}, {"ar+"}, {}));
  comp_old.transform(state_old);
  Field3D nu_old =
      get<Field3D>(state_old["species"]["ar+"]["collision_frequencies"]["ar+_ar_coll"]);
  const BoutReal expected_old =
      sqrt(0.3 / (1.0 * 40.0) + 0.3 / (1.0 * 40.0)) * 1e19 * 5e-19 * 1.0;
  ASSERT_DOUBLE_EQ(nu_old(1, 1, 0), expected_old);

  Options opt_new = make_options(true);
  BraginskiiCollisions comp_new("test", opt_new, nullptr);
  Options state_new = make_state();
  comp_new.declareAllSpecies(SpeciesInformation({}, {"ar"}, {"ar+"}, {}));
  comp_new.transform(state_new);
  Field3D nu_new =
      get<Field3D>(state_new["species"]["ar+"]["collision_frequencies"]["ar+_ar_coll"]);
  const BoutReal factor = nu_new(1, 1, 0) / nu_old(1, 1, 0);
  EXPECT_NEAR(factor, 2.0 * sqrt(8.0 / PI), 1e-9); // = 3.192

  Field3D nu_ni_old =
      get<Field3D>(state_old["species"]["ar"]["collision_frequencies"]["ar_ar+_coll"]);
  Field3D nu_ni_new =
      get<Field3D>(state_new["species"]["ar"]["collision_frequencies"]["ar_ar+_coll"]);
  EXPECT_NEAR(nu_ni_new(1, 1, 0) / nu_ni_old(1, 1, 0), 2.0 * sqrt(8.0 / PI), 1e-9);
}
