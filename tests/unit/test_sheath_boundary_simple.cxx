
#include "gtest/gtest.h"
#include <cmath>

#include "fake_mesh_fixture.hxx"
#include "hermes_utils.hxx"
#include "test_extras.hxx" // FakeMesh

#include <bout/bout_types.hxx>
#include <bout/constants.hxx>
#include <bout/field_factory.hxx> // For generating functions
#include <bout/output.hxx>
#include <bout/sys/range.hxx>

#include "../../include/sheath_boundary_simple.hxx"

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using SheathBoundarySimpleTest = FakeMeshFixture;

// Reuse the tests for SheathBoundary that will still apply to SheathBoundarySimple

TEST_F(SheathBoundarySimpleTest, DontSetPotential) {
  const WithQuietOutput quiet{output_info};
  Options options;
  options["units"]["eV"] = 1.0; // Voltage normalisation
  options["test"]["lower_y"] = false;
  options["test"]["upper_y"] = true;
  options["test"]["diagnose"] = true;

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

  Options state{{"species",
                 {// Electrons
                  {"e", {{"density", N_in}, {"temperature", Te_in}, {"velocity", 0.0}}},
                  // Ion species
                  {"h+",
                   {{"density", si_in * N_in},
                    {"temperature", Ti_in},
                    {"AA", 1.0},
                    {"charge", Zi_in},
                    {"velocity", 0.0}}}}}};

  component.declareAllSpecies({"e", "h+"});
  component.transform(state);

  // Should have calculated, but not set potential
  ASSERT_FALSE(state["fields"].isSet("phi"));

  auto& electrons = state["species"]["e"];
  auto& ions = state["species"]["h+"];

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

  state["Nnorm"] = 1.0;
  state["rho_s0"] = 1.0;
  state["Omega_ci"] = 1.0;
  state["Tnorm"] = 1.0;

  component.outputVars(state);
  ASSERT_TRUE(state.isSet("Ee_sheath"));
  ASSERT_TRUE(state.isSet("Eh+_sheath"));
  ASSERT_TRUE(state.isSet("Sh+_sheath"));

  ASSERT_FALSE(state.isSet("Se_sheath"));

  const Field3D nothing = Field3D{}.setDirectionY(YDirectionType::Aligned);
  EXPECT_NEAR(state["Ee_sheath"].as<Field3D>(nothing)[upper], -8.577, 1e-2);
  EXPECT_NEAR(state["Eh+_sheath"].as<Field3D>(nothing)[upper], -11.696, 1e-2);
  EXPECT_NEAR(state["Sh+_sheath"].as<Field3D>(nothing)[upper], -1.1139, 1e-2);
}

TEST_F(SheathBoundarySimpleTest, PotentialChangeSymmetricOnYBoundaries) {
  // Use a y-symmetric density (no y dependence) and vary Te/Ti over x and z.
  // With symmetric input, the ion-SEE-induced change in floating potential should be
  // identical at both y boundaries.

  Options field_options;
  Field3D N = FieldFactory::get()->create3D("1 + x", &field_options, mesh);
  Field3D Te =
      FieldFactory::get()->create3D("2.0 * (1 + 0.5*sin(z)*x)", &field_options, mesh);
  Field3D Ti =
      FieldFactory::get()->create3D("3.0 * (1 + 0.25*cos(z)*x)", &field_options, mesh);

  const BoutReal Zi = 1.1;
  const BoutReal si = 0.5;

  Options state{{"species",
                 {// Electrons
                  {"e", {{"density", N}, {"temperature", Te}, {"velocity", 0.0}}},
                  // Ion species
                  {"h+",
                   {{"density", si * N},
                    {"temperature", Ti},
                    {"AA", 1.0},
                    {"charge", Zi},
                    {"velocity", 0.0}}}}}};

  Options state_iee = state.copy();

  // Baseline: no ion-induced SEE
  {
    Options options{{"test", {{"always_set_phi", true}}}};
    options["units"]["eV"] = 1.0; // Voltage normalisation

    SheathBoundarySimple component("test", options, nullptr);
    component.declareAllSpecies({"e", "h+"});
    component.transform(state);
  }

  // With ion-induced SEE: use low thresholds so the yield is non-zero.
  {
    Options options{{"test",
                     {{"always_set_phi", true},
                      {"ion_ee_gamma_max", 0.5},
                      {"ion_ee_E_th", 0.1},
                      {"ion_ee_E_max", 10.0},
                      {"ion_ee_p", 1.0}}}};
    options["units"]["eV"] = 1.0; // Voltage normalisation

    SheathBoundarySimple component("test", options, nullptr);
    component.declareAllSpecies({"e", "h+"});
    component.transform(state_iee);
  }

  const Field3D phi0 = state["fields"]["phi"];
  const Field3D phi1 = state_iee["fields"]["phi"];

  for (int ix = mesh->xstart; ix <= mesh->xend; ++ix) {
    for (int kz = mesh->zstart; kz <= mesh->zend; ++kz) {
      const BoutReal dphi_lower = phi1(ix, mesh->ystart, kz) - phi0(ix, mesh->ystart, kz);
      const BoutReal dphi_upper = phi1(ix, mesh->yend, kz) - phi0(ix, mesh->yend, kz);
      ASSERT_NEAR(dphi_lower, dphi_upper, 1e-10);
    }
  }
}

TEST_F(SheathBoundarySimpleTest, CalculatePotential) {
  const WithQuietOutput quiet{output_info};
  Options options{{"test", {{"always_set_phi", true}}}};
  options["units"]["eV"] = 1.0; // Voltage normalisation

  SheathBoundarySimple component("test", options, nullptr);
  component.declareAllSpecies({"e", "h+"});

  Field3D N = FieldFactory::get()->create3D("1 + y", &options, mesh);
  BoutReal Te = 2.0;
  BoutReal Ti = 3.0;
  BoutReal Zi = 1.1;
  const BoutReal si = 0.5;

  Options state{{"species",
                 {// Electrons
                  {"e", {{"density", N}, {"temperature", Te}, {"velocity", 0.0}}},
                  // Ion species
                  {"h+",
                   {{"density", si * N},
                    {"temperature", Ti},
                    {"AA", 1.0},
                    {"charge", Zi},
                    {"velocity", 0.0}}}}}};

  component.transform(state);

  ASSERT_TRUE(state["fields"].isSet("phi"));

  // Golden answer
  const BoutReal phi_ref = 5.9177;

  Field3D phi = state["fields"]["phi"];

  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
    for (int jz = 0; jz < mesh->LocalNz; jz++) {
      ASSERT_NEAR(phi_ref, phi(r.ind, mesh->yend, jz), 1e-3);
    }
  }
}
