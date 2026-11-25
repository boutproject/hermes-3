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
#if CHECKLEVEL >= 1
  // Throws exception due to pressure_source
  EXPECT_THROW(component.finally(state), BoutException);
#endif

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

// keep_solver_pressure: same RHS as the default, but the evolving field keeps the solver
// value in cells with N < density_floor (by default it is replaced by N * T).
TEST_F(EvolvePressureTest, KeepSolverPressure) {
  // Boundary regions, so that setBoundaryTo() acts (the fixture creates none)
  static_cast<FakeMesh*>(bout::globals::mesh)->createBoundaries();

  auto density = []() {
    return makeField<Field3D>(
        [](Field3D::ind_type& i) {
          return 0.01 + 0.07 * i.x() + 0.03 * i.y() + 0.002 * i.z();
        },
        bout::globals::mesh);
  };
  auto pressure = []() {
    return makeField<Field3D>(
        [](Field3D::ind_type& i) {
          return 0.05 + 0.02 * i.x() + 0.011 * i.y() * i.y() + 0.003 * i.z() * i.z();
        },
        bout::globals::mesh);
  };
  auto run = [&](bool keep, Field3D& ddt_P, Field3D& P_after) {
    FakeSolver solver;
    Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                    {"i",
                     {{"AA", 2.0},
                      {"charge", 1.0},
                      {"thermal_conduction", false},
                      {"density_floor", 0.1},
                      {"hyper_z", 0.1},
                      {"keep_solver_pressure", keep}}}};
    EvolvePressure component("i", options, &solver);
    solver.setState({{"Pi", pressure()}}); // as the solver does before each RHS call
    Options state{{"species", {{"i", {{"density", density()}}}}}};
    component.transform(state);
    component.finally(state);
    ddt_P = copy(solver.getTimeDerivs()["Pi"].as<Field3D>());
    P_after = copy(solver.getState()["Pi"].as<Field3D>());
  };
  Field3D ddt_default, P_default, ddt_keep, P_keep;
  run(false, ddt_default, P_default);
  run(true, ddt_keep, P_keep);

  const Field3D P0 = pressure();
  const Field3D N = density();
  int dilute = 0;
  BOUT_FOR_SERIAL(i, P0.getRegion("RGN_NOBNDRY")) {
    ASSERT_EQ(ddt_default[i], ddt_keep[i]) << i; // same RHS, bit for bit
    ASSERT_EQ(P0[i], P_keep[i]) << i;            // solver value kept
    const BoutReal Nlim = softFloor(N[i], 0.1);
    EXPECT_NEAR(P0[i] * N[i] / Nlim, P_default[i], 1e-14) << i; // default: N * T
    if (N[i] / Nlim < 0.99) {
      ++dilute;
    }
  }
  EXPECT_GT(dilute, 0);
}
