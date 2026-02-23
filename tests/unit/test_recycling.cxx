
#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/recycling.hxx"

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
using RecyclingTest = FakeMeshFixture;

TEST_F(RecyclingTest, CreateComponent) {
  Options options;
  options["units"]["eV"] = 5;           // Normalisation temperature
  options["recycling"]["species"] = ""; // No species to recycle

  Recycling component("recycling", options, nullptr);
}

// Increase velocity, pressure and temperature
// and check that recycled density and energy sources increase
TEST_F(RecyclingTest, RecycleSourceChange) {
  Options options;
  options["units"]["eV"] = 1.0;
  options["test"]["species"] = "d+";
  options["d+"]["recycle_as"] = "d";
  options["d+"]["target_recycle"] = true;

  Options state1;
  Options state2;

  state1["species"]["d+"]["AA"] = 2.0;
  state1["species"]["d+"]["density"] = 1.0;

  state1["species"]["d"]["AA"] = 2.0;
  state1["species"]["d"]["density"] = 1.0;
  state1["species"]["d"]["velocity"] = 1.0;
  state1["species"]["d"]["temperature"] = 1.0;
  state1["species"]["d"]["pressure"] = 1.0;

  state2 = state1.copy();
  state1["species"]["d+"]["velocity"] = 1.0;
  state1["species"]["d+"]["temperature"] = 1.0;
  state1["species"]["d+"]["pressure"] = 1.0;

  state2["species"]["d+"]["velocity"] = 2.0;
  state2["species"]["d+"]["pressure"] = 2.0;
  state2["species"]["d+"]["temperature"] = 2.0;

  Recycling component("test", options, nullptr);
  component.declareAllSpecies({"d+", "d"});
  component.transform(state1);
  component.transform(state2);

  Field3D density_source1 = state1["species"]["d"]["density_source"];
  Field3D density_source2 = state2["species"]["d"]["density_source"];

  Field3D energy_source1 = state1["species"]["d"]["energy_source"];
  Field3D energy_source2 = state2["species"]["d"]["energy_source"];

  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
    for (int jz = 0; jz < mesh->LocalNz; jz++) {

      auto i = indexAt(density_source1, r.ind, mesh->yend, jz);

      ASSERT_GT(density_source2[i], density_source1[i]);
      ASSERT_GT(energy_source2[i], energy_source1[i]);
    }
  }
}

// Test that enabling the pump leads to a non-zero pump neutral sink
// on targets
TEST_F(RecyclingTest, TargetPumpSource) {
  Options options;
  options["units"]["eV"] = 1.0;
  options["test"]["species"] = "d+";
  options["d+"]["recycle_as"] = "d";
  options["d+"]["target_recycle"] = true;
  options["d+"]["neutral_pump"] = true;
  options["d+"]["pump_recycle_multiplier"] = 0.5;
  options["d+"]["diagnose"] = true;
  options["d"]["diagnose"] = true;

  Options state;

  state["species"]["d+"]["AA"] = 2.0;
  state["species"]["d+"]["density"] = 1.0;
  state["species"]["d+"]["velocity"] = 1.0;
  state["species"]["d+"]["temperature"] = 1.0;
  state["species"]["d+"]["pressure"] = 1.0;

  state["species"]["d"]["AA"] = 2.0;
  state["species"]["d"]["density"] = 1.0;
  state["species"]["d"]["velocity"] = 1.0;
  state["species"]["d"]["temperature"] = 1.0;
  state["species"]["d"]["pressure"] = 1.0;

  static_cast<FakeMesh*>(mesh)->setGridDataSource(
      new FakeGridDataSource({{"is_pump", 1.0}}));

  Recycling component("test", options, nullptr);
  component.declareAllSpecies({"d+", "d"});
  component.transform(state);

  Options outputs = {
      {"Tnorm", 1.0},
      {"Nnorm", 1.0},
      {"Omega_ci", 1.0},
  };

  component.outputVars(outputs);

  ASSERT_TRUE(outputs.isSet("Sd_pump"));
  ASSERT_TRUE(outputs.isSet("Ed_pump"));

  auto pump_density_source = get<Field3D>(outputs["Sd_pump"]);
  auto pump_energy_source = get<Field3D>(outputs["Ed_pump"]);
  auto is_pump = get<Field3D>(outputs["is_pump"]);

  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
    for (int jz = 0; jz < mesh->LocalNz; jz++) {

      auto i = indexAt(pump_density_source, r.ind, mesh->yend, jz);

      ASSERT_LT(pump_density_source[i], 0.0);
      ASSERT_LT(pump_energy_source[i], 0.0);
    }
  }
}
