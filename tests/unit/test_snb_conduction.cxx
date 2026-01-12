#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/snb_conduction.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

#include <bout/field_factory.hxx> // For generating functions

// Reuse the "standard" fixture for FakeMesh
using SNBConductionTest = FakeMeshFixture;

TEST_F(SNBConductionTest, CreateComponent) {
  Options options{
      {"units",
       {{"meters", 1.0}, {"eV", 1.0}, {"inv_meters_cubed", 1e19}, {"seconds", 1e-6}}}};
  SNBConduction component("test", options, nullptr);
}

TEST_F(SNBConductionTest, Transform) {
  Options options{
      {"units",
       {{"meters", 1.0}, {"eV", 1.0}, {"inv_meters_cubed", 1e19}, {"seconds", 1e-6}}}};
  SNBConduction component("test", options, nullptr);

  Options state{
      {"units",
       {{"meters", 1.0}, {"eV", 1.0}, {"inv_meters_cubed", 1e19}, {"seconds", 1e-6}}},
      {"species", {{"e", {{"temperature", 1.0}, {"density", 1.0}}}}}};
  component.transform(state);

  ASSERT_TRUE(state["species"]["e"].isSet("energy_source"));
  // Zero temperature gradient everywhere -> No divergence of heat flux
  auto source = get<Field3D>(state["species"]["e"]["energy_source"]);
  BOUT_FOR_SERIAL(i, source.getRegion("RGN_NOBNDRY")) {
    ASSERT_LT(abs(source[i]), 1e-20);
  }
}

TEST_F(SNBConductionTest, OutputDiagnose) {
  Options options{
      {"units",
       {{"meters", 1.0}, {"eV", 1.0}, {"inv_meters_cubed", 1e19}, {"seconds", 1e-6}}}};
  options["test"]["diagnose"] = true;
  SNBConduction component("test", options, nullptr);

  Options state{
      {"units",
       {{"meters", 1.0}, {"eV", 1.0}, {"inv_meters_cubed", 1e19}, {"seconds", 1e-6}}},
      {"species", {{"e", {{"temperature", 1.0}, {"density", 1.0}}}}}};
  component.transform(state);

  Options outputs = {{"Tnorm", 1.0}, {"Nnorm", 1.0}, {"Omega_ci", 1.0}};
  component.outputVars(outputs);

  ASSERT_TRUE(outputs.isSet("Div_Q_SH"));
  ASSERT_TRUE(outputs.isSet("Div_Q_SNB"));

  // Zero temperature gradient everywhere -> No divergence of heat flux
  auto Div_Q_SH = get<Field3D>(outputs["Div_Q_SH"]);
  auto Div_Q_SNB = get<Field3D>(outputs["Div_Q_SNB"]);
  BOUT_FOR_SERIAL(i, Div_Q_SH.getRegion("RGN_NOBNDRY")) {
    ASSERT_LT(abs(Div_Q_SH[i]), 1e-20);
    ASSERT_LT(abs(Div_Q_SNB[i]), 1e-20);
  }
}
