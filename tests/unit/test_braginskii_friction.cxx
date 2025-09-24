
#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh
#include "fake_mesh_fixture.hxx"

#include "../../include/braginskii_collisions.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using BraginskiiFrictionTest = FakeMeshFixture;


TEST_F(BraginskiiFrictionTest, OnlyElectrons) {
  Options options;

  options["units"]["eV"] = 1.0;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["inv_meters_cubed"] = 1.0;
  
  BraginskiiCollisions component("test", options, nullptr);

  Options state;
  state["species"]["e"]["density"] = 1e19;
  state["species"]["e"]["temperature"] = 10.;
  state["species"]["e"]["velocity"] = 1.;

  component.transform(state);

  // A species can't exert friction on itself, so momentum and energy transfer won't be set or will be 0.
  if (state["species"]["e"]["momentum_source"].isSet()) {
    ASSERT_FLOAT_EQ(get<Field3D>(state["species"]["e"]["momentum_source"])(0, 0, 0), 0.);
    ASSERT_FLOAT_EQ(get<Field3D>(state["species"]["e"]["energy_source"])(0, 0, 0), 0.);
  } else {
    ASSERT_FALSE(state["species"]["e"]["energy_source"].isSet());
  }
}

TEST_F(BraginskiiFrictionTest, TwoComovingSpeciesCharged) {
  Options options;

  options["units"]["eV"] = 1.0;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["inv_meters_cubed"] = 1.0;
  
  BraginskiiCollisions component("test", options, nullptr);

  // State with two species, both the same but half the density
  Options state;
  state["species"]["s1"]["density"] = 5e18; // Half density
  state["species"]["s1"]["temperature"] = 10;
  state["species"]["s1"]["charge"] = 1;
  state["species"]["s1"]["AA"] = 2;
  state["species"]["s1"]["velocity"] = 1;

  state["species"]["s2"] = state["species"]["s1"].copy();

  // Run calculations
  component.transform(state);
  ASSERT_TRUE(state["species"]["s1"].isSet("momentum_source"));
  ASSERT_TRUE(state["species"]["s2"].isSet("momentum_source"));
  ASSERT_TRUE(state["species"]["s1"].isSet("energy_source"));
  ASSERT_TRUE(state["species"]["s2"].isSet("energy_source"));

  Field3D ms1 = get<Field3D>(state["species"]["s1"]["momentum_source"]);
  Field3D ms2 = get<Field3D>(state["species"]["s2"]["momentum_source"]);
  Field3D es1 = get<Field3D>(state["species"]["s1"]["energy_source"]);
  Field3D es2 = get<Field3D>(state["species"]["s2"]["energy_source"]);
  
  BOUT_FOR_SERIAL(i, ms1.getRegion("RGN_ALL")) {
    // If the species have the same velocities, there should be no friction
    ASSERT_DOUBLE_EQ(ms1[i], 0.);
    ASSERT_DOUBLE_EQ(ms2[i], 0.);
    ASSERT_DOUBLE_EQ(es1[i], 0.);
    ASSERT_DOUBLE_EQ(es2[i], 0.);
  }
}


TEST_F(BraginskiiFrictionTest, TwoSpeciesCharged) {
  Options options;

  options["units"]["eV"] = 1.0;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["inv_meters_cubed"] = 1.0;
  
  BraginskiiCollisions component("test", options, nullptr);

  // State with two species, both the same but half the density
  Options state;
  state["species"]["s1"]["density"] = 5e18; // Half density
  state["species"]["s1"]["temperature"] = 10;
  state["species"]["s1"]["charge"] = 1;
  state["species"]["s1"]["AA"] = 2;

  state["species"]["s2"] = state["species"]["s1"].copy();

  state["species"]["s1"]["velocity"] = 1;
  state["species"]["s2"]["velocity"] = 2;

  // Run calculations
  component.transform(state);

  Field3D ms1 = get<Field3D>(state["species"]["s1"]["momentum_source"]);
  Field3D ms2 = get<Field3D>(state["species"]["s2"]["momentum_source"]);
  Field3D es1 = get<Field3D>(state["species"]["s1"]["energy_source"]);
  Field3D es2 = get<Field3D>(state["species"]["s2"]["energy_source"]);

  BOUT_FOR_SERIAL(i, ms1.getRegion("RGN_ALL")) {
    // The species should have opposite momentum sources
    ASSERT_NE(ms1[i], 0.);
    ASSERT_DOUBLE_EQ(ms1[i], -ms2[i]);
    // The species should have equal frictional heating
    ASSERT_NE(es1[i], 0.);
    ASSERT_DOUBLE_EQ(es1[i], es2[i]);
  }
}


TEST_F(BraginskiiFrictionTest, DoubleRelativeVelocities) {
  Options options;

  options["units"]["eV"] = 1.0;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["inv_meters_cubed"] = 1.0;
  
  BraginskiiCollisions component("test", options, nullptr);

  // State with two species, both the same but half the density
  Options state1, state2;
  state1["species"]["s1"]["density"] = 5e18; // Half density
  state1["species"]["s1"]["temperature"] = 10;
  state1["species"]["s1"]["charge"] = 1;
  state1["species"]["s1"]["AA"] = 2;

  state1["species"]["s2"] = state1["species"]["s1"].copy();
  state2 = state1.copy();

  state1["species"]["s1"]["velocity"] = 1;
  state1["species"]["s2"]["velocity"] = 2;
  state2["species"]["s1"]["velocity"] = 1;
  state2["species"]["s2"]["velocity"] = 3;

  // Run calculations
  component.transform(state1);
  component.transform(state2);

  Field3D ms11 = get<Field3D>(state1["species"]["s1"]["momentum_source"]);
  Field3D ms21 = get<Field3D>(state2["species"]["s1"]["momentum_source"]);
  Field3D es11 = get<Field3D>(state1["species"]["s1"]["energy_source"]);
  Field3D es21 = get<Field3D>(state2["species"]["s1"]["energy_source"]);
  
  BOUT_FOR_SERIAL(i, ms11.getRegion("RGN_ALL")) {
    // The relative velocity in the second state is double that in the
    // first, so the friction should be double too.
    ASSERT_NE(ms11[i], 0.);
    ASSERT_DOUBLE_EQ(2. * ms11[i], ms21[i]);
    // The the frictional heating depends on the relative velocity
    // squared, so that of the second species should be quadruple that
    // of the first.
    ASSERT_NE(es11[i], 0.);
    ASSERT_DOUBLE_EQ(4. * es11[i], es21[i]);
  }
}

TEST_F(BraginskiiFrictionTest, TwoSpeciesNoHeating) {
  Options options;

  options["units"]["eV"] = 1.0;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["inv_meters_cubed"] = 1.0;
  options["test"]["frictional_heating"] = false;
  
  BraginskiiCollisions component("test", options, nullptr);

  // State with two species, both the same but half the density
  Options state;
  state["species"]["s1"]["density"] = 5e18; // Half density
  state["species"]["s1"]["temperature"] = 10;
  state["species"]["s1"]["charge"] = 1;
  state["species"]["s1"]["AA"] = 2;

  state["species"]["s2"] = state["species"]["s1"].copy();

  state["species"]["s1"]["velocity"] = 1;
  state["species"]["s2"]["velocity"] = 2;

  // Run calculations
  component.transform(state);

  Field3D es1 = get<Field3D>(state["species"]["s1"]["energy_source"]);
  Field3D es2 = get<Field3D>(state["species"]["s2"]["energy_source"]);
  
  BOUT_FOR_SERIAL(i, es1.getRegion("RGN_ALL")) {
    // Frictional heating is turned off, so this should be zero
    ASSERT_DOUBLE_EQ(es1[i], 0.);
    ASSERT_DOUBLE_EQ(es2[i], 0.);
  }
}
