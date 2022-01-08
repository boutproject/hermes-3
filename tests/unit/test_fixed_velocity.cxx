
#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh

#include "../../include/fixed_velocity.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using FixedVelocityTest = FakeMeshFixture;

TEST_F(FixedVelocityTest, CreateComponent) {
  Options options = {{"units", {{"meters", 1.0},
                                {"seconds", 1.0}}},
                     {"test", {{"velocity", 1.0}}}};

  FixedVelocity component("test", options, nullptr);
}

TEST_F(FixedVelocityTest, MustSetVelocity) {
  Options options = {{"units", {{"meters", 1.0},
                                {"seconds", 1.0}}},
                     {"test", {{"charge", 1.0}}}};

  // Density isn't set, so this should throw
  ASSERT_THROW(FixedVelocity component("test", options, nullptr), BoutException);
}

TEST_F(FixedVelocityTest, SetValues) {
  Options options = {{"units", {{"meters", 10.0},
                                {"seconds", 2.0}}},
                     {"e", {{"velocity", 3}}}};

  FixedVelocity component("e", options, nullptr);

  Options state = {{"species", {{"e", {{"AA", 2},
                                       {"density", 7}}}}}};

  component.transform(state);
  // Check that velocity has been normalised
  ASSERT_TRUE(IsFieldEqual(state["species"]["e"]["velocity"].as<Field3D>(), 3 / (10. / 2),
                           "RGN_ALL"));
  // Check that momentum is calculated
  ASSERT_TRUE(IsFieldEqual(state["species"]["e"]["momentum"].as<Field3D>(),
                           2 * 7 * 3 / (10. / 2), "RGN_ALL",1e-8));
}
