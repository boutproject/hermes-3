#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/kmodel.hxx"

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
using KmodelTest = FakeMeshFixture;

TEST_F(KmodelTest, CreateComponent) {
  FakeSolver solver;

  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";
  Options options{{"units", {{"seconds", 1.0}, {"Tesla", 1.0}, {"meters", 1.0}}}};

  Kmodel component("test", options, &solver);

  Options state = solver.getState();

  ASSERT_TRUE(state.isSet("k"));
}

TEST_F(KmodelTest, Transform) {
  FakeSolver solver;

  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";
  Options options{{"units", {{"seconds", 1.0}, {"Tesla", 1.0}, {"meters", 1.0}}}};

  Permissions permissions{{readOnly("fields"), readWrite("fields:D_k")}};

  Kmodel component("test", options, &solver);

  Options state = solver.getState();

  Field3D A = FieldFactory::get()->create3D("1 + y * (x - 0.5)", &options, mesh);
  mesh->communicate(A);
  mesh->getCoordinates()->Bxy = 1.0;
  state["sound_speed"] = A;
  state["fastest_wave"] = A;

  component.transform(state);
  ASSERT_TRUE(state.isSet("k"));
  ASSERT_TRUE(state["fields"].isSet("D_k"));
}
