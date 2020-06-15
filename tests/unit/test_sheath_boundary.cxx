
#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh

#include "../../include/sheath_boundary.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using SheathBoundaryTest = FakeMeshFixture;

TEST_F(SheathBoundaryTest, CreateComponent) {
  Options options;
  options["units"]["meters"] = 1.0;
  options["test"]["connection_length"] = 10;
  
  SheathBoundary component("test", options, nullptr);
}
