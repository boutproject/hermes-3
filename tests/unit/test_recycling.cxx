
#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh

#include "../../include/recycling.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

#include <field_factory.hxx>  // For generating functions

// Reuse the "standard" fixture for FakeMesh
using RecyclingTest = FakeMeshFixture;

TEST_F(RecyclingTest, CreateComponent) {
  Options options;
  options["recycling"]["species"] = "";  // No species to recycle

  Recycling component("recycling", options, nullptr);
}

