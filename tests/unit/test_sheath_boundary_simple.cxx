
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

namespace {
// Electrons with an energy source already set (as if from an earlier component),
// ions without one, so both ways into the sheath are covered
Options sheathEnergyState() {
  Options field_options;
  Field3D N = FieldFactory::get()->create3D("1 + y", &field_options, mesh);
  Field3D Te = FieldFactory::get()->create3D("2 + 0.5*x", &field_options, mesh);
  Field3D Ti = FieldFactory::get()->create3D("3 + 0.2*y", &field_options, mesh);
  Field3D Se = FieldFactory::get()->create3D("0.1 * (1 + x + y)", &field_options, mesh);

  return Options{{"species",
                  {{"e",
                    {{"density", N},
                     {"temperature", Te},
                     {"pressure", N * Te},
                     {"velocity", 0.0},
                     {"energy_source", Se}}},
                   {"h+",
                    {{"density", N},
                     {"temperature", Ti},
                     {"pressure", N * Ti},
                     {"AA", 2.0},
                     {"charge", 1.0},
                     {"velocity", 0.0},
                     {"momentum", 0.0}}}}}};
}

Options sheathOutputs() {
  return Options{{"Nnorm", 1.0}, {"Omega_ci", 1.0}, {"Tnorm", 1.0}};
}
} // namespace

TEST_F(SheathBoundarySimpleTest, SheathExtraIsEnergySourceChange) {
  Options options{{"test", {{"diagnose", true}, {"gamma_e", 4.5}, {"gamma_i", 2.5}}}};
  options["units"]["eV"] = 1.0;

  SheathBoundarySimple component("test", options, nullptr);
  component.declareAllSpecies({"e", "h+"});

  Options state = sheathEnergyState();
  // from a fresh state, a copy would share data and the sheath adds in place
  const Field3D Se_before =
      sheathEnergyState()["species"]["e"]["energy_source"].as<Field3D>();
  component.transform(state);

  Options outputs = sheathOutputs();
  component.outputVars(outputs);

  // electrons had a source already, so compare the change
  const Field3D Se_after = get<Field3D>(state["species"]["e"]["energy_source"]);
  const Field3D Ee_extra = get<Field3D>(outputs["Ee_sheath_extra"]);
  const Field3D dSe = Se_after - Se_before;
  EXPECT_TRUE(IsFieldEqual(Ee_extra, dSe, "RGN_NOBNDRY", 1e-12));

  // ions had none, so it's the whole energy_source
  const Field3D Si_after = get<Field3D>(state["species"]["h+"]["energy_source"]);
  const Field3D Ei_extra = get<Field3D>(outputs["Eh+_sheath_extra"]);
  EXPECT_TRUE(IsFieldEqual(Ei_extra, Si_after, "RGN_NOBNDRY", 0.0));

  ASSERT_EQ("h+", outputs["Eh+_sheath_extra"].attributes["species"].as<std::string>());

  // only the cells next to the sheath, and not the same as the total
  const Field3D Ei_total = get<Field3D>(outputs["Eh+_sheath"]);
  for (int ix = mesh->xstart; ix <= mesh->xend; ++ix) {
    for (int kz = mesh->zstart; kz <= mesh->zend; ++kz) {
      for (int jy = mesh->ystart; jy <= mesh->yend; ++jy) {
        const bool at_sheath = (jy == mesh->ystart) or (jy == mesh->yend);
        if (at_sheath) {
          EXPECT_NE(Ei_extra(ix, jy, kz), 0.0);
          EXPECT_NE(Ei_extra(ix, jy, kz), Ei_total(ix, jy, kz));
          EXPECT_NE(Ee_extra(ix, jy, kz), 0.0);
        } else {
          EXPECT_EQ(Ei_extra(ix, jy, kz), 0.0);
          EXPECT_EQ(Ee_extra(ix, jy, kz), 0.0);
        }
      }
    }
  }
}

TEST_F(SheathBoundarySimpleTest, DiagnoseDoesNotChangeState) {
  Options state_off = sheathEnergyState();
  Options state_on = sheathEnergyState();

  Options outputs_off = sheathOutputs();
  {
    Options options{{"test", {{"diagnose", false}}}};
    options["units"]["eV"] = 1.0;
    SheathBoundarySimple component("test", options, nullptr);
    component.declareAllSpecies({"e", "h+"});
    component.transform(state_off);
    component.outputVars(outputs_off);
  }
  {
    Options options{{"test", {{"diagnose", true}}}};
    options["units"]["eV"] = 1.0;
    SheathBoundarySimple component("test", options, nullptr);
    component.declareAllSpecies({"e", "h+"});
    component.transform(state_on);
  }

  // nothing extra saved when it's off
  EXPECT_FALSE(outputs_off.isSet("Ee_sheath_extra"));
  EXPECT_FALSE(outputs_off.isSet("Eh+_sheath_extra"));

  // and the state is the same either way, guard cells included
  for (const std::string sp : {"e", "h+"}) {
    for (const std::string var : {"density", "temperature", "pressure", "velocity",
                                  "energy_source", "energy_flow_ylow"}) {
      const Field3D on = get<Field3D>(state_on["species"][sp][var]);
      const Field3D off = get<Field3D>(state_off["species"][sp][var]);
      EXPECT_TRUE(IsFieldEqual(on, off, "RGN_ALL", 0.0)) << sp << ":" << var;
    }
  }
  const Field3D NVi_on = get<Field3D>(state_on["species"]["h+"]["momentum"]);
  const Field3D NVi_off = get<Field3D>(state_off["species"]["h+"]["momentum"]);
  EXPECT_TRUE(IsFieldEqual(NVi_on, NVi_off, "RGN_ALL", 0.0));
}
