#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/penalisation.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using PenalisationTest = FakeMeshFixture;

TEST_F(PenalisationTest, CreateComponent) {
  FakeSolver solver;

  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";
  Options options{{"units", {{"seconds", 1.0}, {"Tesla", 1.0}, {"meters", 1.0}}}};

  Penalisation component("test", options, &solver);

  Options state = solver.getState();
}

TEST_F(PenalisationTest, Transform) {
  FakeSolver solver;

  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";

  Options options{{"units", {{"seconds", 1.0}, {"Tesla", 1.0}, {"meters", 1.0}}}};

  Penalisation component("test", options, &solver);

  Options state = solver.getState();

  component.transform(state);

  ASSERT_TRUE(state["fields"].isSet("penalisation_wall"));
}
