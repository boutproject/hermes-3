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

using PressureHeatingTest = FakeMeshFixture_tmpl<5, 9, 3, false, 2, 2>;

TEST_F(PressureHeatingTest, DiagnosticHeatingMatchesAppliedEnergyCorrection) {
  getMesh().ix_separatrix = 0;
  Field3D velocity{0.0};
  for (const auto& i : velocity.getRegion("RGN_ALL")) {
    velocity[i] = 0.1 * i.y() * i.y();
  }
  const Options state{{"fastest_wave", 3.0},
                      {"species",
                       {{"e",
                         {{"AA", 1.0},
                          {"density", 1.0},
                          {"pressure", 1.0},
                          {"temperature", 1.0},
                          {"velocity", velocity}}}}}};
  for (bool fixflux : {false, true}) {
    BoutReal previous_rhs = 0.0;
    BoutReal previous_available = 0.0;
    for (bool heating : {false, true}) {
      FakeSolver solver;
      Options options{
          {"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
          {"e",
           {{"thermal_conduction", false},
            {"diagnose", true},
            {"numerical_viscous_heating", heating},
            {"fix_momentum_boundary_flux", fixflux}}}};
      EvolvePressure component("e", options, &solver);
      component.finally(state);
      Options output{{"Nnorm", 1.0}, {"Tnorm", 1.0}, {"rho_s0", 1.0}, {"Omega_ci", 1.0}};
      component.outputVars(output);
      const BoutReal available =
          output["balance_pressure_e_thermal_e_numerical_heat_available"];
      const BoutReal rhs = output["balance_pressure_e_thermal_e_rhs"];
      const BoutReal applied = output["balance_pressure_e_thermal_e_numerical_heat"];
      if (heating) {
        EXPECT_NEAR(available, previous_available, 1e-12);
        EXPECT_NEAR(applied, available, 1e-12);
        EXPECT_NEAR(rhs - previous_rhs, available, 1e-12);
      } else {
        EXPECT_GT(std::abs(available), 1e-8);
        EXPECT_EQ(applied, 0.0);
        previous_rhs = rhs;
        previous_available = available;
      }
    }
  }
}

TEST_F(EvolvePressureTest, MmsForcingIsNotAddedAsPhysicalPressureSource) {
  FakeSolver solver;
  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"solver", {{"mms", true}}},
                  {"Pi", {{"source", 1.0}}},
                  {"i", {{"AA", 2.0}, {"charge", 1.0}, {"thermal_conduction", false}}}};
  EvolvePressure component("i", options, &solver);
  EXPECT_FALSE(options["Pi"]["source"].valueUsed());
  Options state{{"species", {{"i", {{"density", 1.0}, {"energy_source", 0.75}}}}}};
  component.transform(state);
  component.finally(state);
  Field3D derivative = solver.getTimeDerivs()["Pi"].as<Field3D>();
  BOUT_FOR_SERIAL(i, derivative.getRegion("RGN_NOBNDRY")) {
    EXPECT_DOUBLE_EQ(0.5, derivative[i]);
  }
}

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
