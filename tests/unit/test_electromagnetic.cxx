#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/electromagnetic.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

#include <bout/field_factory.hxx> // For generating functions

// Reuse the "standard" fixture for FakeMesh
using ElectromagneticTest = FakeMeshFixture;

#if not BOUT_USE_METRIC_3D

TEST_F(ElectromagneticTest, CreateComponent) {
  Options options = {
      {"units", {{"Tesla", 1.0}, {"eV", 1.0}, {"inv_meters_cubed", 1e19}}}};

  Electromagnetic component("test", options, nullptr);
}

TEST_F(ElectromagneticTest, TransformNoChargedSpecies) {
  Options options = {
      {"units", {{"Tesla", 1.0}, {"eV", 1.0}, {"inv_meters_cubed", 1e19}}}};

  Electromagnetic component("test", options, nullptr);

  // Transform with no charged species
  Options state = {{"species", {{"d", {{"AA", 2.0}}}}}};

  component.transform(state);

  // Apar should be set
  EXPECT_TRUE(state["fields"]["Apar"].isSet());

  // Apar should be zero
  auto Apar = get<Field3D>(state["fields"]["Apar"]);
  BOUT_FOR_SERIAL(i, Apar.getRegion("RGN_NOBNDRY")) { ASSERT_DOUBLE_EQ(Apar[i], 0.0); }

  // Apar_flutter should not be set
  EXPECT_FALSE(state["fields"]["Apar_flutter"].isSet());
}

TEST_F(ElectromagneticTest, FlutterSetsField) {
  Options options = {{"units", {{"Tesla", 1.0}, {"eV", 1.0}, {"inv_meters_cubed", 1e19}}},
                     {"test", {{"magnetic_flutter", true}}}};

  Electromagnetic component("test", options, nullptr);

  // Transform with no charged species
  Options state = {{"species", {{"d", {{"AA", 2.0}}}}}};

  component.transform(state);

  // Apar should be set
  EXPECT_TRUE(state["fields"]["Apar"].isSet());

  // Apar_flutter should be set
  EXPECT_TRUE(state["fields"]["Apar_flutter"].isSet());

  // Apar_flutter should be zero
  auto Apar_flutter = get<Field3D>(state["fields"]["Apar_flutter"]);
  BOUT_FOR_SERIAL(i, Apar_flutter.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(Apar_flutter[i], 0.0);
  }
}
#else

TEST_F(ElectromagneticTest, CreateComponent) {
  Options options = {
      {"units", {{"Tesla", 1.0}, {"eV", 1.0}, {"inv_meters_cubed", 1e19}}}};

  ASSERT_THROW(Electromagnetic component("test", options, nullptr),
	       BoutException);
}

#endif
