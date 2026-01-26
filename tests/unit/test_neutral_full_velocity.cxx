#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/neutral_full_velocity.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using NeutralFullVelocityTest = FakeMeshFixture;

TEST_F(NeutralFullVelocityTest, CreateComponentRequiresRxy) {
  FakeSolver solver;

  static_cast<FakeMesh*>(bout::globals::mesh)
      ->setGridDataSource(new FakeGridDataSource{{// Missing Rxy
                                                  {"Zxy", 0.0},
                                                  {"hthe", 1.0},
                                                  {"Bpxy", 1.0}}});

  Options options{{"units",
                   {{"seconds", 1.0},
                    {"Tesla", 1.0},
                    {"meters", 1.0},
                    {"inv_meters_cubed", 1e19},
                    {"eV", 100}}},
                  {"d", {{"AA", 2.0}}}};

  // Need Rxy
  EXPECT_THROW(NeutralFullVelocity component("d", options, &solver), BoutException);
}

TEST_F(NeutralFullVelocityTest, CreateComponentRequiresZxy) {
  FakeSolver solver;

  static_cast<FakeMesh*>(bout::globals::mesh)
      ->setGridDataSource(new FakeGridDataSource{{{"Rxy", 0.0},
                                                  // Missing Zxy
                                                  {"hthe", 1.0},
                                                  {"Bpxy", 1.0}}});

  Options options{{"units",
                   {{"seconds", 1.0},
                    {"Tesla", 1.0},
                    {"meters", 1.0},
                    {"inv_meters_cubed", 1e19},
                    {"eV", 100}}},
                  {"d", {{"AA", 2.0}}}};

  // Need Zxy
  EXPECT_THROW(NeutralFullVelocity component("d", options, &solver), BoutException);
}

TEST_F(NeutralFullVelocityTest, CreateComponentRequiresHthe) {
  FakeSolver solver;

  static_cast<FakeMesh*>(bout::globals::mesh)
      ->setGridDataSource(new FakeGridDataSource{{{"Rxy", 0.0},
                                                  {"Zxy", 1.0},
                                                  // Missing hthe
                                                  {"Bpxy", 1.0}}});

  Options options{{"units",
                   {{"seconds", 1.0},
                    {"Tesla", 1.0},
                    {"meters", 1.0},
                    {"inv_meters_cubed", 1e19},
                    {"eV", 100}}},
                  {"d", {{"AA", 2.0}}}};

  // Need hthe
  EXPECT_THROW(NeutralFullVelocity component("d", options, &solver), BoutException);
}

TEST_F(NeutralFullVelocityTest, CreateComponentRequiresBpxy) {
  FakeSolver solver;

  static_cast<FakeMesh*>(bout::globals::mesh)
      ->setGridDataSource(new FakeGridDataSource{{
          {"Rxy", 0.0}, {"Zxy", 1.0}, {"hthe", 1.0},
          // Missing Bpxy
      }});

  Options options{{"units",
                   {{"seconds", 1.0},
                    {"Tesla", 1.0},
                    {"meters", 1.0},
                    {"inv_meters_cubed", 1e19},
                    {"eV", 100}}},
                  {"d", {{"AA", 2.0}}}};

  // Need Bpxy
  EXPECT_THROW(NeutralFullVelocity component("d", options, &solver), BoutException);
}

TEST_F(NeutralFullVelocityTest, CreateComponent) {
  FakeSolver solver;

  static_cast<FakeMesh*>(bout::globals::mesh)
      ->setGridDataSource(new FakeGridDataSource{
          {{"Rxy", 1.0}, {"Zxy", 0.0}, {"hthe", 1.0}, {"Bpxy", 1.0}}});

  Options options{{"units",
                   {{"seconds", 1.0},
                    {"Tesla", 1.0},
                    {"meters", 1.0},
                    {"inv_meters_cubed", 1e19},
                    {"eV", 100}}},
                  {"d", {{"AA", 2.0}}}};

  NeutralFullVelocity component("d", options, &solver);

  Options state = solver.getState();

  ASSERT_TRUE(state.isSet("Nd"));
  ASSERT_TRUE(state.isSet("Pd"));
  ASSERT_TRUE(state.isSet("Vdx"));
  ASSERT_TRUE(state.isSet("Vdy"));
  ASSERT_TRUE(state.isSet("Vdz"));
}

TEST_F(NeutralFullVelocityTest, Transform) {
  FakeSolver solver;

  static_cast<FakeMesh*>(bout::globals::mesh)
      ->setGridDataSource(new FakeGridDataSource{
          {{"Rxy", 1.0}, {"Zxy", 0.0}, {"hthe", 1.0}, {"Bpxy", 1.0}}});

  Options options{{"units",
                   {{"seconds", 1.0},
                    {"Tesla", 1.0},
                    {"meters", 1.0},
                    {"inv_meters_cubed", 1e19},
                    {"eV", 100}}},
                  {"neutral", {{"AA", 2.0}}}};

  NeutralFullVelocity component("neutral", options, &solver);

  Options state;
  component.transform(state);

  ASSERT_TRUE(state["species"]["neutral"].isSet("density"));
  ASSERT_DOUBLE_EQ(get<BoutReal>(state["species"]["neutral"]["AA"]), 2.0);
  ASSERT_TRUE(state["species"]["neutral"].isSet("pressure"));
  ASSERT_TRUE(state["species"]["neutral"].isSet("momentum"));
  ASSERT_TRUE(state["species"]["neutral"].isSet("velocity"));
  ASSERT_TRUE(state["species"]["neutral"].isSet("temperature"));
}

TEST_F(NeutralFullVelocityTest, OutputVars) {
  FakeSolver solver;

  static_cast<FakeMesh*>(bout::globals::mesh)
      ->setGridDataSource(new FakeGridDataSource{
          {{"Rxy", 1.0}, {"Zxy", 0.0}, {"hthe", 1.0}, {"Bpxy", 1.0}}});

  Options options{{"units",
                   {{"seconds", 1.0},
                    {"Tesla", 1.0},
                    {"meters", 1.0},
                    {"inv_meters_cubed", 1e19},
                    {"eV", 100}}},
                  {"neutral", {{"AA", 2.0}}}};

  NeutralFullVelocity component("neutral", options, &solver);

  Options state;
  component.transform(state);

  Options outputs {
    {"Nnorm", get<BoutReal>(options["units"]["inv_meters_cubed"])},
    {"Tnorm", get<BoutReal>(options["units"]["eV"])},
    {"Omega_ci", 1./get<BoutReal>(options["units"]["seconds"])},
    {"Cs0", get<BoutReal>(options["units"]["meters"]) /
     get<BoutReal>(options["units"]["seconds"])}};
  component.outputVars(outputs);

  ASSERT_TRUE(outputs.isSet("Urx"));
  ASSERT_TRUE(outputs.isSet("Tyr"));
}
