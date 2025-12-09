
#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/evolve_density.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using EvolveDensityTest = FakeMeshFixture;

TEST_F(EvolveDensityTest, CreateComponent) {
  FakeSolver solver;

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"i", {{"AA", 2.0}, {"charge", 1.0}}}};
  EvolveDensity component("i", options, &solver);

  Options state = solver.getState();

  ASSERT_TRUE(state.isSet("Ni"));
}

TEST_F(EvolveDensityTest, Transform) {
  FakeSolver solver;

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"i", {{"AA", 2.0}, {"charge", 1.0}}}};
  EvolveDensity component("i", options, &solver);

  Options state;
  component.performAllSubstitutions({"i"});
  component.transform(state);

  Options& species = state["species"]["i"];
  ASSERT_TRUE(species.isSet("density"));
  ASSERT_TRUE(species.isSet("AA"));
  ASSERT_TRUE(species.isSet("charge"));
}

TEST_F(EvolveDensityTest, Finally) {
  FakeSolver solver;

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"i", {{"AA", 2.0}, {"charge", 1.0}}}};
  EvolveDensity component("i", options, &solver);

  // Call the finally() method with a density source
  const Options state = {
      {"species", {{"i", {{"density", 1.0}, {"density_source", 0.5}}}}}};
  component.finally(state);

  Options ddt = solver.getTimeDerivs();

  ASSERT_TRUE(ddt.isSet("Ni"));
  Field3D ddt_Ni = ddt["Ni"].as<Field3D>();

  BOUT_FOR_SERIAL(i, ddt_Ni.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(0.5, ddt_Ni[i]);
  }
}
