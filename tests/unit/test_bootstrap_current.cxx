#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/bootstrap_current.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using BootstrapCurrentTest = FakeMeshFixture;

TEST_F(BootstrapCurrentTest, CreateComponent) {
  Options options;

  options["units"]["eV"] = 1.0;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["inv_meters_cubed"] = 1e19;
  options["units"]["Tesla"] = 1.0;

  BootstrapCurrent component("test", options, nullptr);
}
