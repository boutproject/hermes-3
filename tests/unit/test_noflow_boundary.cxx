
#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh

#include "../../include/noflow_boundary.hxx"


/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

#include <bout/constants.hxx>
#include <field_factory.hxx>  // For generating functions

// Reuse the "standard" fixture for FakeMesh
using NoFlowBoundaryTest = FakeMeshFixture;

TEST_F(NoFlowBoundaryTest, CreateComponent) {
  Options options;

  NoFlowBoundary component("test", options, nullptr);
}
