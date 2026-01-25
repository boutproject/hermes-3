#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh
#include "fake_solver.hxx"
#include "fake_mesh_fixture.hxx"

#include "../../include/vorticity.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using VorticityTest = FakeMeshFixture;

TEST_F(VorticityTest, CreateComponent) {
  FakeSolver solver;

  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";
  Options options {{"units",
                    {{"seconds", 1.0},
                     {"Tesla", 1.0},
                     {"meters", 1.0}}}};

  Vorticity component("test", options, &solver);

  Options state = solver.getState();

  ASSERT_TRUE(state.isSet("Vort"));
}

TEST_F(VorticityTest, Transform) {
  FakeSolver solver;

  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";
  Options options {{"units",
                    {{"seconds", 1.0},
                     {"Tesla", 1.0},
                     {"meters", 1.0}}}};

  Vorticity component("test", options, &solver);

  Options state;
  component.transform(state);

  ASSERT_TRUE(state["fields"].isSet("vorticity"));
  ASSERT_TRUE(state["fields"].isSet("phi"));
}

TEST_F(VorticityTest, TransformNoDiamagnetic) {
  FakeSolver solver;
  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";

  Options options {{"units",
                    {{"seconds", 1.0},
                     {"Tesla", 1.0},
                     {"meters", 1.0}}},
                   {"test",
                    {{"diamagnetic", false}}}};

  Vorticity component("test", options, &solver);
  component.declareAllSpecies({{"d+"},
                               {}});

  Options state{{"species",
                 {{"d+",
                   {{"AA", 2.0},
                    {"charge", 1.0},
                    {"pressure", 1.0}}}}}};

  component.transform(state);

  // No diamagnetic terms so values are not set
  ASSERT_FALSE(state["species"]["d+"].isSet("energy_source"));
  ASSERT_FALSE(state["fields"].isSet("DivJdia"));
}

TEST_F(VorticityTest, TransformWithCurlb_B) {
  FakeSolver solver;
  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";

  // Add curvature data
  static_cast<FakeMesh*>(bout::globals::mesh)->setGridDataSource(
                                                                 new FakeGridDataSource {{{"bxcvx", 1.0},
                                                                                          {"bxcvy", 0.0},
                                                                                          {"bxcvz", 0.0}}});

  Options options {{"units",
                    {{"seconds", 1.0},
                     {"Tesla", 1.0},
                     {"meters", 1.0}}},
                   {"test",
                    {{"diamagnetic", true}}}};

  Vorticity component("test", options, &solver);
  component.declareAllSpecies({{"d+"},
                               {}});

  Options state{{"species",
                 {{"d+",
                   {{"AA", 2.0},
                    {"charge", 1.0},
                    {"pressure", 1.0}}}}}};

  component.transform(state);

  // Diamagnetic terms set both the divergence of the current,
  // and energy exchange with charged species
  ASSERT_TRUE(state["species"]["d+"].isSet("energy_source"));
  ASSERT_TRUE(state["fields"].isSet("DivJdia"));
}

TEST_F(VorticityTest, KinematicViscosity) {
  FakeSolver solver;

  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";
  Options options {{"units",
                    {{"seconds", 1.0},
                     {"Tesla", 1.0},
                     {"meters", 1.0}}},
                   {"test",
                    {{"viscosity", 1.0},
                     {"diamagnetic", false}}}};

  Vorticity component("test", options, &solver);
  component.declareAllSpecies({{"e"},
                               {"d"},
                               {"d+"},
                               {}});

  Options state {{"species",
                  {{"d+",
                    {{"AA", 2.0},
                     {"charge", 1.0},
                     {"density", 1.0}}},
                   {"d",
                    {{"AA", 2.0},
                     {"density", 1.0}}}}}};
  component.transform(state);

  ASSERT_TRUE(state["species"]["d+"].isSet("energy_source"));
  ASSERT_FALSE(state["species"]["d"].isSet("energy_source"));
}

TEST_F(VorticityTest, calculatePihat) {
  FakeSolver solver;
  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";

  Options options {{"units",
                    {{"seconds", 1.0},
                     {"Tesla", 1.0},
                     {"meters", 1.0}}},
                   {"test",
                    {{"average_atomic_mass", 2.5}}}};

  Vorticity component("test", options, &solver);

  Options state {{"species",
                  {{"d+",
                    {{"pressure", 1.2},
                     {"AA", 2.0},
                     {"charge", 1.4}}}}}};
  Permissions permissions {{readOnly("species"),
                              readWrite("species:d+:energy_source")}};
  GuardedOptions guarded_state {&state, &permissions};

  Field3D Pi_hat = component.calculatePihat(guarded_state["species"]);

  ASSERT_TRUE(IsFieldEqual(Pi_hat,
                           1.2 * 2.0 / 2.5 / 1.4, "RGN_NOBNDRY"));
}

TEST_F(VorticityTest, calculateDivJdia) {
  FakeSolver solver;
  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";

  Options options {{"units",
                    {{"seconds", 1.0},
                     {"Tesla", 1.0},
                     {"meters", 1.0}}},
                   {"test",
                    {{"average_atomic_mass", 2.5}}}};

  Vorticity component("test", options, &solver);

  Options state {{"species",
                  {{"d+",
                    {{"pressure", 1.2},
                     {"AA", 2.0},
                     {"charge", 1.4}}}}}};
  Permissions permissions {{readOnly("species"),
                              readWrite("species:d+:energy_source")}};
  GuardedOptions guarded_state {&state, &permissions};

  Field3D DivJdia = component.calculateDivJdia(0.0, guarded_state["species"]);

  ASSERT_TRUE(IsFieldEqual(DivJdia,
                           0.0, "RGN_NOBNDRY"));
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["d+"]["energy_source"]),
                           0.0, "RGN_NOBNDRY"));
}
