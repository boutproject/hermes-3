
#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh
#include "fake_mesh_fixture.hxx"

#include "../../include/braginskii_heat_exchange.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using BraginskiiHeatExchangeTest = FakeMeshFixture;


TEST_F(BraginskiiHeatExchangeTest, OnlyElectrons) {
  Options options;

  options["units"]["eV"] = 1.0;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["inv_meters_cubed"] = 1.0;
  
  BraginskiiHeatExchange component("test", options, nullptr);

  Options state;
  state["species"]["e"]["density"] = 1e19;
  state["species"]["e"]["temperature"] = 10.;
  state["species"]["e"]["velocity"] = 1.;
  state["species"]["e"]["AA"] = 1./1836;
  state["species"]["e"]["collision_frequencies"]["e_e_coll"] = 1.;

  component.transform(state);

  // A species can't exchange heat with itself
  if (state["species"]["e"]["energy_source"].isSet()) {
    ASSERT_FLOAT_EQ(get<Field3D>(state["species"]["e"]["energy_source"])(0,0,0), 0.);
  }
}

TEST_F(BraginskiiHeatExchangeTest, TwoEqualTempSpeciesCharged) {
  Options options;

  options["units"]["eV"] = 1.0;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["inv_meters_cubed"] = 1.0;
  
  BraginskiiHeatExchange component("test", options, nullptr);

  Options state;
  state["species"]["s1"]["density"] = 5e18; // Half density
  state["species"]["s1"]["temperature"] = 10;
  state["species"]["s1"]["charge"] = 1;
  state["species"]["s1"]["AA"] = 2;
  state["species"]["s1"]["velocity"] = 1;

  state["species"]["s2"] = state["species"]["s1"].copy();
  state["species"]["s1"]["collision_frequencies"]["s1_s1_coll"] = 1.;
  state["species"]["s1"]["collision_frequencies"]["s1_s2_coll"] = 0.5;
  state["species"]["s2"]["collision_frequencies"]["s2_s2_coll"] = 0.5;
  state["species"]["s2"]["collision_frequencies"]["s2_s1_coll"] = 0.25;

  // Run calculations
  component.transform(state);
  ASSERT_TRUE(state["species"]["s1"].isSet("energy_source"));
  ASSERT_TRUE(state["species"]["s2"].isSet("energy_source"));

  Field3D es1 = get<Field3D>(state["species"]["s1"]["energy_source"]);
  Field3D es2 = get<Field3D>(state["species"]["s2"]["energy_source"]);
  
  BOUT_FOR_SERIAL(i, es1.getRegion("RGN_ALL")) {
    // If the species have the same temperature, there should be no heat exchange
    ASSERT_DOUBLE_EQ(es1[i], 0.);
    ASSERT_DOUBLE_EQ(es2[i], 0.);
  }
}


TEST_F(BraginskiiHeatExchangeTest, TwoSpeciesCharged) {
  Options options;

  options["units"]["eV"] = 1.0;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["inv_meters_cubed"] = 1.0;
  
  BraginskiiHeatExchange component("test", options, nullptr);

  Options state;
  state["species"]["s1"]["density"] = 5e18; // Half density
  state["species"]["s1"]["charge"] = 1;
  state["species"]["s1"]["AA"] = 2;
  state["species"]["s1"]["velocity"] = 1;

  state["species"]["s2"] = state["species"]["s1"].copy();
  state["species"]["s1"]["collision_frequencies"]["s1_s1_coll"] = 1.;
  state["species"]["s1"]["collision_frequencies"]["s1_s2_coll"] = 0.5;
  state["species"]["s2"]["collision_frequencies"]["s2_s2_coll"] = 0.5;
  state["species"]["s2"]["collision_frequencies"]["s2_s1_coll"] = 0.25;

  state["species"]["s1"]["temperature"] = 10;
  state["species"]["s2"]["temperature"] = 20;

  // Run calculations
  component.transform(state);

  Field3D es1 = get<Field3D>(state["species"]["s1"]["energy_source"]);
  Field3D es2 = get<Field3D>(state["species"]["s2"]["energy_source"]);

  BOUT_FOR_SERIAL(i, es1.getRegion("RGN_ALL")) {
    ASSERT_NE(es1[i], 0.);
    // The species should have opposite heat exchange
    ASSERT_DOUBLE_EQ(es1[i], -es2[i]);
  }
}

TEST_F(BraginskiiHeatExchangeTest, DoubleCollisionRates) {
  Options options;

  options["units"]["eV"] = 1.0;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["inv_meters_cubed"] = 1.0;
  
  BraginskiiHeatExchange component("test", options, nullptr);

  Options state1, state2;
  state1["species"]["s1"]["density"] = 5e18; // Half density
  state1["species"]["s1"]["charge"] = 0;
  state1["species"]["s1"]["AA"] = 2;
  state1["species"]["s1"]["velocity"] = 1;
  state1["species"]["s2"] = state1["species"]["s1"].copy();
  state1["species"]["s1"]["temperature"] = 10;
  state1["species"]["s2"]["temperature"] = 20;
  state1["species"]["s1"]["collision_frequencies"]["s1_s1_coll"] = 1.;
  state1["species"]["s2"]["collision_frequencies"]["s2_s2_coll"] = 0.25;

  state2 = state1.copy();
  state1["species"]["s1"]["collision_frequencies"]["s1_s2_coll"] = 0.5;
  state1["species"]["s2"]["collision_frequencies"]["s2_s1_coll"] = 0.5;
  state2["species"]["s1"]["collision_frequencies"]["s1_s2_coll"] = 1.0;
  state2["species"]["s2"]["collision_frequencies"]["s2_s1_coll"] = 1.0;

  // Run calculations
  component.transform(state1);
  component.transform(state2);

  Field3D es11 = get<Field3D>(state1["species"]["s1"]["energy_source"]);
  Field3D es21 = get<Field3D>(state2["species"]["s1"]["energy_source"]);

  const BoutReal factor = 2. * 2. / sqrt(3.);
  BOUT_FOR_SERIAL(i, es11.getRegion("RGN_ALL")) {
    ASSERT_NE(es11[i], 0.);
    ASSERT_DOUBLE_EQ(2 * es11[i], es21[i]);
  }
}
