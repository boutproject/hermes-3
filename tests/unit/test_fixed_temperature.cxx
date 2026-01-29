
#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/fixed_temperature.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using FixedTemperatureTest = FakeMeshFixture;

TEST_F(FixedTemperatureTest, CreateComponent) {
  Options options = {{"units", {{"eV", 1.0}}}, {"test", {{"temperature", 1.0}}}};

  FixedTemperature component("test", options, nullptr);
}

TEST_F(FixedTemperatureTest, MustSetTemperature) {
  Options options = {{"units", {{"eV", 1.0}}}, {"test", {{"density", 1.0}}}};

  ASSERT_THROW(FixedTemperature component("test", options, nullptr), BoutException);
}

TEST_F(FixedTemperatureTest, SetTemperature) {
  Options options = {{"units", {{"eV", 5.0}}}, {"test", {{"temperature", 1.4}}}};

  FixedTemperature component("test", options, nullptr);

  Options state;
  component.transform(state);

  ASSERT_TRUE(state["species"]["test"].isSet("temperature"));
  // Note: Sets normalised value
  ASSERT_TRUE(
      IsFieldEqual(get<Field3D>(state["species"]["test"]["temperature"]), 1.4 / 5));
  // Pressure not set because density wasn't set
  ASSERT_FALSE(state["species"]["test"].isSet("pressure"));
}

TEST_F(FixedTemperatureTest, SetPressure) {
  Options options = {{"units", {{"eV", 5.0}}}, {"test", {{"temperature", 1.4}}}};

  FixedTemperature component("test", options, nullptr);

  Options state{{"species", {{"test", {{"density", 2.0}}}}}};
  component.transform(state);

  EXPECT_TRUE(state["species"]["test"].isSet("temperature"));
  EXPECT_TRUE(
      IsFieldEqual(get<Field3D>(state["species"]["test"]["temperature"]), 1.4 / 5));
  ASSERT_TRUE(state["species"]["test"].isSet("pressure"));
  ASSERT_TRUE(
      IsFieldEqual(get<Field3D>(state["species"]["test"]["pressure"]), 2.0 * 1.4 / 5));
}

TEST_F(FixedTemperatureTest, outputVarsDiagnose) {
  Options options = {{"units", {{"eV", 5.0}}},
                     {"test", {{"temperature", 1.4}, {"diagnose", true}}}};

  FixedTemperature component("test", options, nullptr);

  Options state{{"species", {{"test", {{"density", 2.0}}}}}};
  component.transform(state);

  Options outputs{{"Tnorm", 5.0}, {"Nnorm", 1e19}};
  component.outputVars(outputs);

  ASSERT_TRUE(outputs.isSet("Ttest"));
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(outputs["Ttest"]), 1.4 / 5));
  ASSERT_TRUE(outputs.isSet("Ptest"));
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(outputs["Ptest"]), 2.0 * 1.4 / 5));
}
