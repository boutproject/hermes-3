
#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/evolve_energy.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using EvolveEnergyTest = FakeMeshFixture;

TEST_F(EvolveEnergyTest, CreateComponent) {
  FakeSolver solver;

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"i", {{"AA", 2.0}, {"charge", 1.0}}}};
  EvolveEnergy component("i", options, &solver);

  Options state = solver.getState();

  ASSERT_TRUE(state.isSet("Ei"));
}

TEST_F(EvolveEnergyTest, Transform) {
  FakeSolver solver;

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"i", {{"AA", 2.0}, {"charge", 1.0}}}};
  EvolveEnergy component("i", options, &solver);

  Options state;

  // Exception if density, velocity and AA are not set
  EXPECT_THROW(component.transform(state), BoutException);

  state = {{"species", {{"i", {{"density", 1.0}, {"velocity", 0.0}, {"AA", 1.0}}}}}};
  component.transform(state);

  Options& species = state["species"]["i"];
  ASSERT_TRUE(species.isSet("pressure"));
  ASSERT_TRUE(species.isSet("temperature"));
}

TEST_F(EvolveEnergyTest, Finally) {
  FakeSolver solver;

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"i", {{"AA", 2.0}, {"charge", 1.0}}}};
  EvolveEnergy component("i", options, &solver);

  Options state = {
      {"species", {{"i", {{"density", 1.0}, {"velocity", 0.0}, {"AA", 1.0}}}}}};
  component.transform(state);

  component.finally(state);

  Options ddt = solver.getTimeDerivs();

  ASSERT_TRUE(ddt.isSet("Ei"));
  ASSERT_TRUE(IsFieldEqual(ddt["Ei"].as<Field3D>(), 0.0));
}
