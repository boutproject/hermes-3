
#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

#include <bout/bout_types.hxx>
#include <bout/output.hxx>

#include "../../include/sheath_boundary_simple.hxx"

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
using SheathBoundarySimpleTest = FakeMeshFixture;

// Reuse the tests for SheathBoundary that will still apply to SheathBoundarySimple

TEST_F(SheathBoundarySimpleTest, DontSetPotential) {
  Options options;
  options["units"]["eV"] = 1.0; // Voltage normalisation

  SheathBoundarySimple component("test", options, nullptr);

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
