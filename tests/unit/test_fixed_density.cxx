
#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh

#include "../../include/fixed_density.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using FixedDensityTest = FakeMeshFixture;

TEST_F(FixedDensityTest, CreateComponent) {
  Options options = {{"units", {{"inv_meters_cubed", 1.0}}},
                     {"test", {{"density", 1.0}}}};

  FixedDensity component("test", options, nullptr);
}
TEST_F(FixedDensityTest, MustSetDensity) {
  Options options = {{"units", {{"inv_meters_cubed", 1.0}}},
                     {"test", {{"charge", 1.0}}}};

  // Density isn't set, so this should throw
  ASSERT_THROW(FixedDensity component("test", options, nullptr), BoutException);
}

TEST_F(FixedDensityTest, SetValues) {
  Options options = {{"units", {{"inv_meters_cubed", 1e10}}},
                     {"e", {{"density", 2.4e11},
                            {"charge", 2},
                            {"AA", 3}}}};

  FixedDensity component("e", options, nullptr);

  Options state;
  component.transform(state);
  ASSERT_FLOAT_EQ(2, state["species"]["e"]["charge"].as<BoutReal>());
  ASSERT_FLOAT_EQ(3, state["species"]["e"]["AA"].as<BoutReal>());
  // Check that density has been normalised
  ASSERT_TRUE(
      IsFieldEqual(state["species"]["e"]["density"].as<Field3D>(), 24.0, "RGN_ALL"));
}
