#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/evolve_momentum.hxx"
#include "../../include/viscosity_precon.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using EvolveMomentumTest = FakeMeshFixture;

TEST_F(EvolveMomentumTest, CreateComponent) {
  FakeSolver solver;

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"i", {{"AA", 2.0}, {"charge", 1.0}}}};
  EvolveMomentum component("i", options, &solver);

  Options state = solver.getState();

  ASSERT_TRUE(state.isSet("NVi"));
}

TEST_F(EvolveMomentumTest, Transform) {
  FakeSolver solver;

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"i", {{"AA", 2.0}, {"charge", 1.0}}}};
  EvolveMomentum component("i", options, &solver);

  Options state;
  // Exception if density is not set
  EXPECT_THROW(component.transform(state), BoutException);
  state = {{"species", {{"i", {{"density", 1.0}, {"AA", 1.0}}}}}};
  component.transform(state);

  Options& species = state["species"]["i"];
  ASSERT_TRUE(species.isSet("velocity"));
  ASSERT_TRUE(species.isSet("momentum"));
}

TEST_F(EvolveMomentumTest, Finally) {
  FakeSolver solver;

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"i", {{"AA", 2.0}, {"charge", 1.0}}}};
  EvolveMomentum component("i", options, &solver);

  const Options state = {{"fastest_wave", 0.0},
                         {"species",
                          {{"i",
                            {{"AA", 1.0},
                             {"density", 1.0},
                             {"momentum", 0.0},
                             {"velocity", 0.0},
                             {"momentum_source", 0.5}}}}}};
  component.finally(state);

  Options ddt = solver.getTimeDerivs();

  ASSERT_TRUE(ddt.isSet("NVi"));
  Field3D ddt_NVi = ddt["NVi"].as<Field3D>();

  BOUT_FOR_SERIAL(i, ddt_NVi.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(0.5, ddt_NVi[i]);
  }
}

TEST_F(EvolveMomentumTest, ViscosityPrecon) {
  namespace pc = hermes::viscosity_precon;
  FakeSolver solver;

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"i", {{"AA", 2.0}, {"charge", 1.0}}}};
  EvolveMomentum component("i", options, &solver);

  Options state = {{"fastest_wave", 0.0},
                   {"species",
                    {{"i",
                      {{"AA", 1.0},
                       {"density", 1.0},
                       {"momentum", 0.0},
                       {"velocity", 0.0},
                       {"momentum_source", 0.5}}}}}};
  component.finally(state);

  // Nothing published: precon leaves ddt unchanged
  pc::frozen.clear();
  component.precon(state, 1.0);
  Field3D ddt_NVi = solver.getTimeDerivs()["NVi"].as<Field3D>();
  BOUT_FOR_SERIAL(i, ddt_NVi.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(0.5, ddt_NVi[i]);
  }

  // Zero viscosity: the inversion reduces to the identity
  pc::frozen["i"] = {Field3D(0.0)};
  component.precon(state, 1.0);
  ddt_NVi = solver.getTimeDerivs()["NVi"].as<Field3D>();
  BOUT_FOR_SERIAL(i, ddt_NVi.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(0.5, ddt_NVi[i]);
  }

  // Unsupported option
  state["scale_timederivs"] = 1.0;
  EXPECT_THROW(component.precon(state, 1.0), BoutException);
  pc::frozen.clear();
}
