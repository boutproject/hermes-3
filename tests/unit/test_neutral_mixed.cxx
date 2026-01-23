
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

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}, {"meters", 1.0}}},
                  {"d", {{"type","neutral_mixed"}, {"AA", 2.0}, {"evolve_momentum", true}}}};
  NeutralMixed component("d", options, &solver);

  Options state = solver.getState();

  ASSERT_TRUE(state.isSet("Nd"));
  ASSERT_TRUE(state.isSet("Pd"));
  ASSERT_TRUE(state.isSet("NVd"));
}

TEST_F(NeutralMixedTest, CreateComponentEvolveMomentumFalse) {
  FakeSolver solver;

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}, {"meters", 1.0}}},
                  {"d", {{"type","neutral_mixed"}, {"AA", 2.0}, {"evolve_momentum", false}}}};
  NeutralMixed component("d", options, &solver);

  Options state = solver.getState();

  ASSERT_TRUE(state.isSet("Nd"));
  ASSERT_TRUE(state.isSet("Pd"));
  ASSERT_FALSE(state.isSet("NVd"));
}

TEST_F(NeutralMixedTest, Transform) {
  FakeSolver solver;

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}, {"meters", 1.0}}},
                  {"d", {{"type","neutral_mixed"}, {"AA", 2.0}, {"evolve_momentum", true}}}};
  NeutralMixed component("d", options, &solver);

  Options state;
  component.transform(state);

  Options& species = state["species"]["d"];
  ASSERT_TRUE(species.isSet("density"));
  ASSERT_TRUE(species.isSet("AA"));
  ASSERT_TRUE(species.isSet("pressure"));
  ASSERT_TRUE(species.isSet("momentum"));
  ASSERT_TRUE(species.isSet("velocity"));
  ASSERT_TRUE(species.isSet("temperature"));
}

TEST_F(NeutralMixedTest, Finally) {
  FakeSolver solver;

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}, {"meters", 1.0}}},
                  {"d", {{"type","neutral_mixed"}, {"AA", 2.0}, {"evolve_momentum", true}}}};
  NeutralMixed component("d", options, &solver);

  // Call the finally() method with a density, energy, and momentum source
  const Options state = {
      {"species", {{"d", {{"density", 1.0}, {"density_source", 0.5},
                         {"pressure", 1.0}, {"energy_source", 1.5},
                         {"momentum", 1.0}, {"momentum_source", 0.75},
                         {"temperature", 1.0}, {"velocity", 1.0}}}}}};
  component.finally(state);

  Options ddt = solver.getTimeDerivs();

  ASSERT_TRUE(ddt.isSet("Nd"));
  Field3D ddt_Nd = ddt["Nd"].as<Field3D>();

  BOUT_FOR_SERIAL(i, ddt_Nd.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(0.5, ddt_Nd[i]);
  }

  ASSERT_TRUE(ddt.isSet("Pd"));
  Field3D ddt_Pd = ddt["Pd"].as<Field3D>();

  BOUT_FOR_SERIAL(i, ddt_Pd.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(1.0, ddt_Pd[i]);
  }

  ASSERT_TRUE(ddt.isSet("NVd"));
  Field3D ddt_NVd = ddt["NVd"].as<Field3D>();

  BOUT_FOR_SERIAL(i, ddt_NVd.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(0.75, ddt_NVd[i]);
  }
}

TEST_F(NeutralMixedTest, FinallyCollisionalityOverride) {
  FakeSolver solver;

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}, {"meters", 1.0}}},
                  {"d", {{"type","neutral_mixed"}, {"AA", 2.0}, {"evolve_momentum", true},
                         {"collisionality_override", 1.0}}}};
  NeutralMixed component("d", options, &solver);

  // Call the finally() method with a density, energy, and momentum source
  const Options state = {
      {"species", {{"d", {{"density", 1.0}, {"density_source", 0.5},
                         {"pressure", 1.0}, {"energy_source", 1.5},
                         {"momentum", 1.0}, {"momentum_source", 0.75},
                         {"temperature", 1.0}, {"velocity", 1.0}}}}}};
  component.finally(state);

  Options ddt = solver.getTimeDerivs();

  ASSERT_TRUE(ddt.isSet("Nd"));
  Field3D ddt_Nd = ddt["Nd"].as<Field3D>();

  BOUT_FOR_SERIAL(i, ddt_Nd.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(0.5, ddt_Nd[i]);
  }

  ASSERT_TRUE(ddt.isSet("Pd"));
  Field3D ddt_Pd = ddt["Pd"].as<Field3D>();

  BOUT_FOR_SERIAL(i, ddt_Pd.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(1.0, ddt_Pd[i]);
  }

  ASSERT_TRUE(ddt.isSet("NVd"));
  Field3D ddt_NVd = ddt["NVd"].as<Field3D>();

  BOUT_FOR_SERIAL(i, ddt_NVd.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(0.75, ddt_NVd[i]);
  }
}

TEST_F(NeutralMixedTest, FinallyEvolveMomentumFalse) {
  FakeSolver solver;

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}, {"meters", 1.0}}},
                  {"d", {{"type","neutral_mixed"}, {"AA", 2.0}, {"evolve_momentum", false}}}};
  NeutralMixed component("d", options, &solver);

  // Call the finally() method with a density, energy, and momentum source
  const Options state = {
      {"species", {{"d", {{"density", 1.0}, {"density_source", 0.5},
                         {"pressure", 1.0}, {"energy_source", 1.5},
                         {"momentum", 1.0}, {"momentum_source", 0.75},
                         {"temperature", 1.0}, {"velocity", 1.0}}}}}};
  component.finally(state);

  Options ddt = solver.getTimeDerivs();

  ASSERT_TRUE(ddt.isSet("Nd"));
  Field3D ddt_Nd = ddt["Nd"].as<Field3D>();

  BOUT_FOR_SERIAL(i, ddt_Nd.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(0.5, ddt_Nd[i]);
  }

  ASSERT_TRUE(ddt.isSet("Pd"));
  Field3D ddt_Pd = ddt["Pd"].as<Field3D>();

  BOUT_FOR_SERIAL(i, ddt_Pd.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(1.0, ddt_Pd[i]);
  }

  ASSERT_FALSE(ddt.isSet("NVd"));
}
