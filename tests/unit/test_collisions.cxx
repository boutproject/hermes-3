
#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh

#include "../../include/collisions.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using CollisionsTest = FakeMeshFixture;

TEST_F(CollisionsTest, CreateComponent) {
  Options options;
  
  options["units"]["eV"] = 1.0;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["inv_meters_cubed"] = 1e19;
  Collisions component("test", options, nullptr);
}

TEST_F(CollisionsTest, OnlyElectrons) {
  Options options;

  options["units"]["eV"] = 1.0;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["inv_meters_cubed"] = 1.0;
  
  Collisions component("test", options, nullptr);

  Options state;
  state["species"]["e"]["density"] = 1e19;
  state["species"]["e"]["temperature"] = 10.;

  component.transform(state);

  ASSERT_TRUE(state["species"]["e"].isSet("collision_frequency"));
}

TEST_F(CollisionsTest, OneOrTwoSpeciesCharged) {
  Options options;

  options["units"]["eV"] = 1.0;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["inv_meters_cubed"] = 1.0;
  
  Collisions component("test", options, nullptr);

  Options state1;
  state1["species"]["s1"]["density"] = 1e19;
  state1["species"]["s1"]["temperature"] = 10;
  state1["species"]["s1"]["charge"] = 1;
  state1["species"]["s1"]["AA"] = 2;

  // State with two species, both the same but half the density
  Options state2;
  state2["species"]["s1"]["density"] = 5e18; // Half density
  state2["species"]["s1"]["temperature"] = 10;
  state2["species"]["s1"]["charge"] = 1;
  state2["species"]["s1"]["AA"] = 2;

  state2["species"]["s2"] = state2["species"]["s1"];

  // Run calculations
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
    ASSERT_TRUE(abs(nu1[i] - nu21[i]) / (nu1[i] + nu21[i]) < 0.05 );
  }
}

TEST_F(CollisionsTest, TnormDependence) {
  // Calculate rates with normalisation factors 1
  Options options {{"units", {{"eV", 1.0},
                              {"meters", 1.0},
                              {"seconds", 1.0},
                              {"inv_meters_cubed", 1.0}}},
                   {"test", {{"electron_neutral", true},
                             {"ion_neutral", true},
                             {"electron_electron", true},
                             {"ion_ion", true},
                             {"neutral_neutral", true}}}};

  Collisions component("test", options, nullptr);

  Options state {{"species", {{"e", {{"density", 1e19},
                                     {"temperature", 10},
                                     {"charge", -1},
                                     {"AA", 1./1836}}},
                              {"d+", {{"density", 2e19},
                                      {"temperature", 20},
                                      {"charge", 1},
                                      {"AA", 2}}},
                              {"d", {{"density", 1e18},
                                     {"temperature", 3},
                                     {"AA", 2}}}}}};

  component.transform(state);

  ASSERT_TRUE(state["species"]["e"].isSet("collision_frequency"));
  ASSERT_TRUE(state["species"]["d"].isSet("collision_frequency"));
  ASSERT_TRUE(state["species"]["d+"].isSet("collision_frequency"));

  // Collision frequencies should be positive, non-zero
  ASSERT_GT(get<Field3D>(state["species"]["e"]["collision_frequency"])(0,0,0), 0.0);
  ASSERT_GT(get<Field3D>(state["species"]["d"]["collision_frequency"])(0,0,0), 0.0);
  ASSERT_GT(get<Field3D>(state["species"]["d+"]["collision_frequency"])(0,0,0), 0.0);

  // Re-calculate with Tnorm != 1
  // To keep frequency normalisation fixed, rho_s0 scales like sqrt(Tnorm)
  const BoutReal Tnorm = 100;
  Options options2 {{"units", {{"eV", Tnorm},
                               {"meters", sqrt(Tnorm)},
                               {"seconds", 1.0},
                               {"inv_meters_cubed", 1.0}}},
                    {"test", {{"electron_neutral", true},
                              {"ion_neutral", true},
                              {"electron_electron", true},
                              {"ion_ion", true},
                              {"neutral_neutral", true}}}};

  Collisions component2("test", options2, nullptr);

  Options state2 {{"species", {{"e", {{"density", 1e19},
                                      {"temperature", 10 / Tnorm},
                                      {"charge", -1},
                                      {"AA", 1./1836}}},
                               {"d+", {{"density", 2e19},
                                       {"temperature", 20 / Tnorm},
                                       {"charge", 1},
                                       {"AA", 2}}},
                               {"d", {{"density", 1e18},
                                      {"temperature", 3 / Tnorm},
                                      {"AA", 2}}}}}};

  component2.transform(state2);

  // Normalised frequencies should be unchanged
  ASSERT_FLOAT_EQ(get<Field3D>(state["species"]["e"]["collision_frequency"])(0,0,0),
                  get<Field3D>(state2["species"]["e"]["collision_frequency"])(0,0,0));
  ASSERT_FLOAT_EQ(get<Field3D>(state["species"]["d"]["collision_frequency"])(0,0,0),
                  get<Field3D>(state2["species"]["d"]["collision_frequency"])(0,0,0));
  ASSERT_FLOAT_EQ(get<Field3D>(state["species"]["d+"]["collision_frequency"])(0,0,0),
                  get<Field3D>(state2["species"]["d+"]["collision_frequency"])(0,0,0));
}
