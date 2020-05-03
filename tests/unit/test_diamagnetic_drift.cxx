
#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh

#include "../../include/diamagnetic_drift.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using DiamagneticDriftTest = FakeMeshFixture;

TEST_F(DiamagneticDriftTest, CreateComponent) {
  Options options;
  
  options["units"]["Tesla"] = 1.0;
  options["units"]["meters"] = 1.0;
  DiamagneticDrift component("test", options, nullptr);
}
