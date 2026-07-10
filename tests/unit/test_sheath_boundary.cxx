
#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

#include <bout/bout_types.hxx>
#include <bout/output.hxx>

#include "../../include/sheath_boundary.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

#include <bout/constants.hxx>
#include <bout/field_factory.hxx> // For generating functions

// Reuse the "standard" fixture for FakeMesh
using SheathBoundaryTest = FakeMeshFixture;

TEST_F(SheathBoundaryTest, DontSetPotential) {
  Options options;
  options["units"]["eV"] = 1.0; // Voltage normalisation

  SheathBoundary component("test", options, nullptr);

  Field3D N = FieldFactory::get()->create3D("1 + y", &options, mesh);
  const BoutReal Te = 2.0;
  const BoutReal Ti = 3.0;
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

  component.declareAllSpecies({"e", "h+"});
  component.transform(state);

  // Should have calculated, but not set potential
  ASSERT_FALSE(state["fields"].isSet("phi"));
}

TEST_F(SheathBoundaryTest, CalculatePotential) {
  Options options{{"test", {{"always_set_phi", true}}}};
  options["units"]["eV"] = 1.0; // Voltage normalisation

  SheathBoundary component("test", options, nullptr);

  Field3D N = FieldFactory::get()->create3D("1 + y", &options, mesh);
  const BoutReal Te = 2.0;
  const BoutReal Ti = 3.0;
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

  component.declareAllSpecies({"e", "h+"});
  component.transform(state);

  // Should have calculated, but not set potential
  ASSERT_TRUE(state["fields"].isSet("phi"));

  // Calculate the expected value of phi
  const BoutReal adiabatic = 5. / 3;
  BoutReal Vzi = sqrt((adiabatic * Ti) + (Zi * Te));
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

TEST_F(SheathBoundaryTest, CalculatePotentialIonSEE) {
  // Same setup as CalculatePotential, but enable ion-induced secondary electron emission
  // and check that the floating potential decreases.

  Field3D N = FieldFactory::get()->create3D("1 + y", &Options::root(), mesh);
  const BoutReal Te = 2.0;
  const BoutReal Ti = 3.0;
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

    SheathBoundary component("test", options, nullptr);
    component.declareAllSpecies({"e", "h+"});
    component.transform(state);
  }

  // With ion-induced SEE: use low thresholds so the yield is non-zero in this unit test
  {
    Options options{{"test",
                     {{"always_set_phi", true},
                      {"ion_ee_gamma_max", 0.5},
                      {"ion_ee_E_th", 0.1},
                      {"ion_ee_E_max", 10.0},
                      {"ion_ee_p", 1.0}}}};
    options["units"]["eV"] = 1.0; // Voltage normalisation

    SheathBoundary component("test", options, nullptr);
    component.declareAllSpecies({"e", "h+"});
    component.transform(state_iee);
  }

  const Field3D phi0 = state["fields"]["phi"];
  const Field3D phi1 = state_iee["fields"]["phi"];

  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
    for (int jz = 0; jz < mesh->LocalNz; jz++) {
      ASSERT_LT(phi1(r.ind, mesh->yend, jz), phi0(r.ind, mesh->yend, jz));
    }
  }
}

TEST_F(SheathBoundaryTest, PotentialChangeSymmetricOnYBoundaries) {
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

    SheathBoundary component("test", options, nullptr);
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

    SheathBoundary component("test", options, nullptr);
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

TEST_F(SheathBoundaryTest, IonSEE) {
  // Test ion-induced secondary electron emission

  const BoutReal Ne = 0.4;
  const BoutReal Te = 2.0;
  const BoutReal Ti = 3.0;
  const BoutReal Zi = 1.1;
  const BoutReal si = 0.5;
  const BoutReal phi = 4.0;

  // Prepare a state and a copy
  Options state{{"fields", {{"phi", phi}}},
                {"species",
                 {// Electrons
                  {"e", {{"density", Ne}, {"temperature", Te}, {"velocity", 0.0}}},
                  // Ion species
                  {"h+",
                   {{"density", si * Ne},
                    {"temperature", Ti},
                    {"AA", 1.0},
                    {"charge", Zi},
                    {"velocity", 0.0}}}}}};
  Options state2 = state.copy();

  {
    // Run without ion-induced secondary electron emission
    Options options;
    options["units"]["eV"] = 1.0; // Voltage normalisation

    SheathBoundary component("test", options, nullptr);

    component.declareAllSpecies({"e", "h+"});
    component.transform(state);

    // Check that the velocity of both species is antisymmetric and non-zero
    const auto Ve = state["species"]["e"]["velocity"].as<Field3D>();
    const auto Vi = state["species"]["h+"]["velocity"].as<Field3D>();

    for (int ix = mesh->xstart; ix <= mesh->xend; ++ix) {
      for (int kz = mesh->zstart; kz <= mesh->zend; ++kz) {
        ASSERT_DOUBLE_EQ(Ve(ix, mesh->ystart - 1, kz), -Ve(ix, mesh->yend + 1, kz));
        ASSERT_DOUBLE_EQ(Vi(ix, mesh->ystart - 1, kz), -Vi(ix, mesh->yend + 1, kz));
      }
    }
  }

  {
    // Enable ion-induced secondary electron emission
    Options options;
    options["units"]["eV"] = 1.0; // Voltage normalisation
    options["test"] = {{"ion_ee_gamma_max", 0.5}};

    SheathBoundary component("test", options, nullptr);

    component.declareAllSpecies({"e", "h+"});
    component.transform(state2);

    // Check that the velocity of both species is antisymmetric and non-zero
    const auto Ve = state2["species"]["e"]["velocity"].as<Field3D>();
    const auto Vi = state2["species"]["h+"]["velocity"].as<Field3D>();

    for (int ix = mesh->xstart; ix <= mesh->xend; ++ix) {
      for (int kz = mesh->zstart; kz <= mesh->zend; ++kz) {
        ASSERT_DOUBLE_EQ(Ve(ix, mesh->ystart - 1, kz), -Ve(ix, mesh->yend + 1, kz));
        ASSERT_DOUBLE_EQ(Vi(ix, mesh->ystart - 1, kz), -Vi(ix, mesh->yend + 1, kz));
      }
    }
  }

  // Check that the ion velocity is unchanged by IEE
  const auto Vi1 = state["species"]["h+"]["velocity"].as<Field3D>();
  const auto Vi2 = state2["species"]["h+"]["velocity"].as<Field3D>();
  for (int ix = mesh->xstart; ix <= mesh->xend; ++ix) {
    for (int kz = mesh->zstart; kz <= mesh->zend; ++kz) {
      ASSERT_DOUBLE_EQ(Vi1(ix, mesh->ystart - 1, kz), Vi2(ix, mesh->ystart - 1, kz));
      ASSERT_DOUBLE_EQ(Vi1(ix, mesh->yend + 1, kz), Vi2(ix, mesh->yend + 1, kz));
    }
  }

  // Check that electron velocity changes so that some electrons are coming into the domain
  const auto Ve1 = state["species"]["e"]["velocity"].as<Field3D>();
  const auto Ve2 = state2["species"]["e"]["velocity"].as<Field3D>();
  for (int ix = mesh->xstart; ix <= mesh->xend; ++ix) {
    for (int kz = mesh->zstart; kz <= mesh->zend; ++kz) {
      ASSERT_GT(Ve2(ix, mesh->ystart - 1, kz), Ve1(ix, mesh->ystart - 1, kz));
      ASSERT_LT(Ve2(ix, mesh->yend + 1, kz), Ve1(ix, mesh->yend + 1, kz));
    }
  }
}

TEST_F(SheathBoundaryTest, IonSEEYieldFunction) {
  Options options;
  options["units"]["eV"] = 1.0; // Voltage normalisation

  // Disabled by default (ion_ee_gamma_max < 0)
  {
    SheathBoundary component("test", options, nullptr);
    ASSERT_DOUBLE_EQ(component.ionSecondaryElectronEmissionGamma(100.0), 0.0);
  }

  // Enabled: check value at E_i = E_max where the model reduces to smooth_turn_on * gamma_max
  {
    options["test"] = {{"ion_ee_gamma_max", 0.5},
                       {"ion_ee_E_th", 50.0},
                       {"ion_ee_E_max", 5e3},
                       {"ion_ee_p", 1.0}};

    SheathBoundary component("test", options, nullptr);

    const BoutReal Ei = 5e3;
    const BoutReal smooth = 0.5 * (1. + std::tanh((Ei - 50.0) / (0.3 * 50.0)));
    const BoutReal expected = smooth * 0.5;

    ASSERT_DOUBLE_EQ(component.ionSecondaryElectronEmissionGamma(Ei), expected);
  }
}
