#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/evolve_pressure.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using EvolvePressureTest = FakeMeshFixture;

TEST_F(EvolvePressureTest, CreateComponent) {
  FakeSolver solver;

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"i", {{"AA", 2.0}, {"charge", 1.0}}}};
  EvolvePressure component("i", options, &solver);

  Options state = solver.getState();

  ASSERT_TRUE(state.isSet("Pi"));
}

TEST_F(EvolvePressureTest, Transform) {
  FakeSolver solver;

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"i", {{"AA", 2.0}, {"charge", 1.0}}}};
  EvolvePressure component("i", options, &solver);

  Options state;
  // Exception if density is not set
  EXPECT_THROW(component.transform(state), BoutException);

  state = {{"species", {{"i", {{"density", 1.0}}}}}};
  component.transform(state);

  Options& species = state["species"]["i"];
  ASSERT_TRUE(species.isSet("pressure"));
  ASSERT_TRUE(species.isSet("temperature"));
}

TEST_F(EvolvePressureTest, Finally) {
  FakeSolver solver;

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"i", {{"AA", 2.0}, {"charge", 1.0}, {"thermal_conduction", false}}}};
  EvolvePressure component("i", options, &solver);

  const Options state = {{"species",
                          {{"i",
                            {{"density", 1.0},
                             {"pressure", 1.0},
                             {"temperature", 1.0},
                             {"pressure_source", 0.5}}}}}};
  // Throws exception due to pressure_source
  EXPECT_THROW(component.finally(state), BoutException);

  const Options state2 = {{"species",
                           {{"i",
                             {{"density", 1.0},
                              {"pressure", 1.0},
                              {"temperature", 1.0},
                              {"energy_source", 0.5}}}}}};
  component.finally(state2);

  Options ddt = solver.getTimeDerivs();

  ASSERT_TRUE(ddt.isSet("Pi"));
  Field3D ddt_Pi = ddt["Pi"].as<Field3D>();

  BOUT_FOR_SERIAL(i, ddt_Pi.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ((2. / 3) * 0.5, ddt_Pi[i]);
  }
}
