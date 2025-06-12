#include "gtest/gtest.h"
#include <cmath>

#include "hermes_utils.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/sheath_boundary.hxx"
#include "../../include/sheath_boundary_simple.hxx"
#include <bout/bout_types.hxx>
#include <bout/output.hxx>
#include <bout/region.hxx>

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

#include <bout/constants.hxx>
#include <bout/field_factory.hxx>  // For generating functions

// Reuse the "standard" fixture for FakeMesh
using SheathBoundaryTest = FakeMeshFixture;

TEST_F(SheathBoundaryTest, CreateComponent) {
  Options options;
  options["units"]["eV"] = 1.0; // Voltage normalisation

  SheathBoundary component("test", options, nullptr);
}

TEST_F(SheathBoundaryTest, DontSetPotential) {
  Options options;
  options["units"]["eV"] = 1.0; // Voltage normalisation

  SheathBoundary component("test", options, nullptr);

  Field3D N = FieldFactory::get()->create3D("1 + y", &options, mesh);
  BoutReal Te = 2.0;
  BoutReal Ti = 3.0;
  BoutReal Zi = 1.1;
  BoutReal si = 0.5;

  Options state {{"species",
                  {// Electrons
                   {"e", {{"density", N},
                          {"temperature", Te},
                          {"velocity", 0.0}}},
                   // Ion species
                   {"h", {{"density", si*N},
                          {"temperature", Ti},
                          {"AA", 1.0},
                          {"charge", Zi},
                          {"velocity", 0.0}}}}}};

  component.transform(state);

  // Should have calculated, but not set potential
  ASSERT_FALSE(state["fields"].isSet("phi"));
}

TEST_F(SheathBoundaryTest, CalculatePotential) {
  Options options{{"test", {{"always_set_phi", true}}}};
  options["units"]["eV"] = 1.0; // Voltage normalisation

  SheathBoundary component("test", options, nullptr);

  Field3D N = FieldFactory::get()->create3D("1 + y", &options, mesh);
  BoutReal Te = 2.0;
  BoutReal Ti = 3.0;
  BoutReal Zi = 1.1;
  BoutReal si = 0.5;

  Options state{{"species",
                 {// Electrons
                  {"e", {{"density", N}, {"temperature", Te}, {"velocity", 0.0}}},
                  // Ion species
                  {"h",
                   {{"density", si * N},
                    {"temperature", Ti},
                    {"AA", 1.0},
                    {"charge", Zi},
                    {"velocity", 0.0}}}}}};

  component.transform(state);

  // Should have calculated, but not set potential
  ASSERT_TRUE(state["fields"].isSet("phi"));

  // Calculate the expected value of phi
  const BoutReal adiabatic = 5./3;
  BoutReal Vzi = sqrt(adiabatic * Ti + Zi * Te);
  BoutReal phi_ref = Te * log(sqrt(Te * SI::Mp / SI::Me / TWOPI) / (si * Zi * Vzi));

  output.write("TEST: {:e} {:e} {:e}\n", Te, si * Zi * Vzi, phi_ref);
  output.write("ION: {:e} {:e} {:e}\n", adiabatic * Ti, Zi * Te * si / (si + 1), Vzi);

  Field3D phi = state["fields"]["phi"];

  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
    for (int jz = 0; jz < mesh->LocalNz; jz++) {
      ASSERT_DOUBLE_EQ(phi_ref, phi(r.ind, mesh->yend, jz));
    }
  }
}

// Reuse the "standard" fixture for FakeMesh
using SheathBoundarySimpleTest = FakeMeshFixture;

TEST_F(SheathBoundarySimpleTest, CreateComponent) {
  WithQuietOutput quiet{output_info};
  Options options;
  options["units"]["eV"] = 1.0; // Voltage normalisation

  SheathBoundarySimple component("test", options, nullptr);
}

TEST_F(SheathBoundarySimpleTest, DontSetPotential) {
  WithQuietOutput quiet{output_info};
  Options options;
  options["units"]["eV"] = 1.0; // Voltage normalisation
  options["test"]["lower_y"] = false;
  options["test"]["upper_y"] = true;

  SheathBoundarySimple component("test", options, nullptr);

  Field3D N_in = FieldFactory::get()->create3D("1 - (1e-10 * exp(y))", &options, mesh);

  const auto lower = indexAt(N_in, 1, mesh->ystart, 0);
  const auto upper = indexAt(N_in, 1, mesh->yend, 0);

  // Make sure boundary is not set
  N_in[upper.yp()] = BoutNaN;

  constexpr BoutReal Te_in = 2.0;
  constexpr BoutReal Ti_in = 3.0;
  constexpr BoutReal Zi_in = 1.1;
  constexpr BoutReal si_in = 0.5;

  Options state {{"species",
                  {// Electrons
                   {"e", {{"density", N_in},
                          {"temperature", Te_in},
                          {"velocity", 0.0}}},
                   // Ion species
                   {"h", {{"density", si_in*N_in},
                          {"temperature", Ti_in},
                          {"AA", 1.0},
                          {"charge", Zi_in},
                          {"velocity", 0.0}}}}}};

  component.transform(state);

  // Should have calculated, but not set potential
  ASSERT_FALSE(state["fields"].isSet("phi"));

  auto& electrons = state["species"]["e"];
  auto& ions = state["species"]["h"];

  ASSERT_TRUE(electrons.isSet("velocity"));
  ASSERT_TRUE(electrons.isSet("density"));
  ASSERT_TRUE(electrons.isSet("temperature"));

  ASSERT_TRUE(ions.isSet("velocity"));
  ASSERT_TRUE(ions.isSet("density"));
  ASSERT_TRUE(ions.isSet("temperature"));
  
  const auto& Ve = electrons["velocity"].as<Field3D>();
  const auto& Ne = electrons["density"].as<Field3D>();
  const auto& Te = electrons["temperature"].as<Field3D>();

  const auto& Vi = ions["velocity"].as<Field3D>();
  const auto& Ni = ions["density"].as<Field3D>();
  const auto& Ti = ions["temperature"].as<Field3D>();

  // Velocity at sheath boundary end should have risen
  EXPECT_GT(Ve[upper.yp()], 2.5);
  EXPECT_GT(Vi[upper.yp()], 2.5);
  // But not at other end
  EXPECT_DOUBLE_EQ(Ve[lower.ym()], 0.0);
  EXPECT_DOUBLE_EQ(Vi[lower.ym()], 0.0);

  // Density and temperature should be set in boundary
  EXPECT_FALSE(std::isnan(Ne[upper.yp()]));
  EXPECT_FALSE(std::isnan(Te[upper.yp()]));
  EXPECT_FALSE(std::isnan(Ni[upper.yp()]));
  EXPECT_FALSE(std::isnan(Ti[upper.yp()]));

  const auto& Ee = electrons["energy_source"].as<Field3D>();
  const auto& Ei = ions["energy_source"].as<Field3D>();

  // Energy should be set at sheath (and be negative)
  EXPECT_LT(Ee[upper], 0.0);
  EXPECT_LT(Ei[upper], 0.0);
  // but not at other end
  EXPECT_DOUBLE_EQ(Ee[lower], 0.0);
  EXPECT_DOUBLE_EQ(Ei[lower], 0.0);

  const auto& Qe = electrons["energy_flow_ylow"].as<Field3D>();
  const auto& Qi = ions["energy_flow_ylow"].as<Field3D>();

  // Energy should be set at sheath (and be positive)
  EXPECT_GT(Qe[upper.yp()], 0.0);
  EXPECT_GT(Qi[upper.yp()], 0.0);
  // but not at other end
  EXPECT_DOUBLE_EQ(Qe[lower], 0.0);
  EXPECT_DOUBLE_EQ(Qi[lower], 0.0);
}

TEST_F(SheathBoundarySimpleTest, CalculatePotential) {
  WithQuietOutput quiet{output_info};
  Options options{{"test", {{"always_set_phi", true}}}};
  options["units"]["eV"] = 1.0; // Voltage normalisation

  SheathBoundarySimple component("test", options, nullptr);

  Field3D N = FieldFactory::get()->create3D("1 + y", &options, mesh);
  BoutReal Te = 2.0;
  BoutReal Ti = 3.0;
  BoutReal Zi = 1.1;
  BoutReal si = 0.5;

  Options state{{"species",
                 {// Electrons
                  {"e", {{"density", N}, {"temperature", Te}, {"velocity", 0.0}}},
                  // Ion species
                  {"h",
                   {{"density", si * N},
                    {"temperature", Ti},
                    {"AA", 1.0},
                    {"charge", Zi},
                    {"velocity", 0.0}}}}}};

  component.transform(state);

  ASSERT_TRUE(state["fields"].isSet("phi"));

  // Golden answer
  BoutReal phi_ref = 5.9177;

  Field3D phi = state["fields"]["phi"];

  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
    for (int jz = 0; jz < mesh->LocalNz; jz++) {
      ASSERT_NEAR(phi_ref, phi(r.ind, mesh->yend, jz), 1e-3);
    }
  }
}
