#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/relax_potential.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using RelaxPotentialTest = FakeMeshFixture;

TEST_F(RelaxPotentialTest, CreateComponent) {
  FakeSolver solver;

  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";
  Options options{{"units",
                   {{"seconds", 1.0},
                    {"Tesla", 1.0},
                    {"meters", 1.0},
                    {"eV", 100.},
                    {"inv_meters_cubed", 1e19}}}};

  RelaxPotential component("test", options, &solver);

  Options state = solver.getState();

  ASSERT_TRUE(state.isSet("Vort"));
  ASSERT_TRUE(state.isSet("phi1"));
}

TEST_F(RelaxPotentialTest, Transform) {
  FakeSolver solver;

  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";
  Options options{{"units",
                   {{"seconds", 1.0},
                    {"Tesla", 1.0},
                    {"meters", 1.0},
                    {"eV", 100.},
                    {"inv_meters_cubed", 1e19}}}};

  RelaxPotential component("test", options, &solver);

  Options state;
  component.transform(state);

  ASSERT_TRUE(state["fields"].isSet("vorticity"));
  ASSERT_TRUE(state["fields"].isSet("phi"));
}
