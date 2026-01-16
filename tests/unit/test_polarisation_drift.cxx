#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh
#include "fake_mesh_fixture.hxx"

#include "../../include/polarisation_drift.hxx"
#include <bout/invert_laplace.hxx>

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using PolarisationDriftTest = FakeMeshFixture;

TEST_F(PolarisationDriftTest, CreateComponent) {
  Options options;

  PolarisationDrift component("test", options, nullptr);
}

TEST_F(PolarisationDriftTest, calcDivJZero) {
  Options options;

  PolarisationDrift component("polarisation_drift", options, nullptr);

  Options state {};
  Permissions permissions {};
  GuardedOptions guarded_state {&state, &permissions};
  Field3D DivJ = component.calcDivJ(guarded_state);
  ASSERT_TRUE(IsFieldEqual(DivJ, 0.0, "RGN_NOBNDRY"));
}

TEST_F(PolarisationDriftTest, calcDivJdia) {
  Options options;

  PolarisationDrift component("polarisation_drift", options, nullptr);

  Options state {{"fields",
                  {{"DivJdia", 1.0}}}};
  Permissions permissions {{readOnly("fields")}};
  GuardedOptions guarded_state {&state, &permissions};
  Field3D DivJ = component.calcDivJ(guarded_state);
  ASSERT_TRUE(IsFieldEqual(DivJ, 1.0, "RGN_NOBNDRY"));
}

TEST_F(PolarisationDriftTest, diamagneticCompression) {
  Options options {{"polarisation_drift",
                    {{"average_atomic_mass", 1.5}}}};

  PolarisationDrift component("polarisation_drift", options, nullptr);

  Options state {{"species",
                  {{"d+",
                    {{"pressure", 1.2},
                     {"AA", 2.0},
                     {"charge", 1.4}}}}}};
  Permissions permissions {{readOnly("species"),
                              readWrite("species:d+:energy_source")}};
  GuardedOptions guarded_state {&state, &permissions};
  component.diamagneticCompression(guarded_state, 3.0);
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["d+"]["energy_source"]),
                           1.2 * 2.0 * 3.0 / 1.5 / 1.4, "RGN_NOBNDRY"));
}

TEST_F(PolarisationDriftTest, calcMassDensityBoussinesq) {
  Options options {{"polarisation_drift",
                    {{"boussinesq", true},
                     {"average_atomic_mass", 1.5}}}};

  PolarisationDrift component("polarisation_drift", options, nullptr);

  Options state {};
  Permissions permissions {};
  GuardedOptions guarded_state {&state, &permissions};
  Field3D mass_density = component.calcMassDensity(guarded_state);
  ASSERT_TRUE(IsFieldEqual(mass_density, 1.5, "RGN_NOBNDRY"));
}

TEST_F(PolarisationDriftTest, calcMassDensityNonBoussinesq) {
  Options options {{"polarisation_drift",
                    {{"boussinesq", false}}}};

  PolarisationDrift component("polarisation_drift", options, nullptr);

  Options state {{"species",
                  {{"d+",
                    {{"AA", 2.0},
                     {"charge", 1.0},
                     {"density", 3.0}}}}}};
  Permissions permissions {{readOnly("species")}};
  GuardedOptions guarded_state {&state, &permissions};
  Field3D mass_density = component.calcMassDensity(guarded_state);
  ASSERT_TRUE(IsFieldEqual(mass_density, 6.0, "RGN_NOBNDRY"));
}

TEST_F(PolarisationDriftTest, polarisationAdvection) {
  Options options;

  PolarisationDrift component("polarisation_drift", options, nullptr);

  Options state {{"species",
                  {{"d+",
                    {{"charge", 1.0},
                     {"AA", 2.0},
                     {"density", 1.0},
                     {"pressure", 1.0},
                     {"momentum", 1.0}}},
                   {"i",
                    {{"charge", 1.0},
                     {"AA", 2.0},
                     {"density", 1.0},
                     {"momentum", 1.0}}},
                   {"d",
                    {{"AA", 2.0}}},
                   {"d2",
                    {{"charge", 0.0},
                     {"AA", 4.0}}},
                   {"massless",
                    {{"charge", 1.0}}}}}};

  Permissions permissions {{
      readOnly("species"),
      readWrite("species:d+:density_source"),
      readWrite("species:d+:energy_source"),
      readWrite("species:d+:momentum_source"),
      readWrite("species:i:density_source"),
      readWrite("species:i:momentum_source"),
    }};
  GuardedOptions guarded_state {&state, &permissions};
  component.polarisationAdvection(guarded_state, 0.0);
  ASSERT_TRUE(state["species"]["d+"].isSet("density_source"));
  ASSERT_TRUE(state["species"]["d+"].isSet("energy_source"));
  ASSERT_TRUE(state["species"]["d+"].isSet("momentum_source"));

  ASSERT_TRUE(state["species"]["i"].isSet("density_source"));
  ASSERT_TRUE(state["species"]["i"].isSet("momentum_source"));

  ASSERT_FALSE(state["species"]["d"].isSet("density_source"));
  ASSERT_FALSE(state["species"]["d"].isSet("energy_source"));
  ASSERT_FALSE(state["species"]["d"].isSet("momentum_source"));
}
