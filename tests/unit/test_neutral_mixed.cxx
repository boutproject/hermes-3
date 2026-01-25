
#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/neutral_mixed.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using NeutralMixedTest = FakeMeshFixture;

TEST_F(NeutralMixedTest, CreateComponent) {
  FakeSolver solver;

  Options options{
      {"units",
       {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}, {"meters", 1.0}}},
      {"n", {{"AA", 2.0}}}};
  NeutralMixed component("n", options, &solver);

  Options state = solver.getState();

  EXPECT_TRUE(state.isSet("Nn"));
  EXPECT_TRUE(state.isSet("Pn"));
  EXPECT_TRUE(state.isSet("NVn"));
}

TEST_F(NeutralMixedTest, Transform) {
  FakeSolver solver;

  Options options{
      {"units",
       {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}, {"meters", 1.0}}},
      {"n", {{"AA", 2.0}}}};
  NeutralMixed component("n", options, &solver);

  Options state;
  component.transform(state);

  auto& species = state["species"]["n"];
  EXPECT_TRUE(species.isSet("density"));
  EXPECT_TRUE(species.isSet("AA"));
  EXPECT_TRUE(species.isSet("pressure"));
  EXPECT_TRUE(species.isSet("momentum"));
  EXPECT_TRUE(species.isSet("velocity"));
  EXPECT_TRUE(species.isSet("temperature"));
}
