
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
  options["units"]["inv_meters_cubed"] = 1e19;
  Collisions component("test", options, nullptr);
}
