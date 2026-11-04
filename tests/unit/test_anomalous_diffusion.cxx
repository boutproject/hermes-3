
#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/anomalous_diffusion.hxx"

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
using AnomalousDiffusionTest = FakeMeshFixture;

TEST_F(AnomalousDiffusionTest, NoDiffusion) {
  Options options;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;

  AnomalousDiffusion component("h", options, nullptr);

  Field3D N = FieldFactory::get()->create3D("1 + y * (x - 0.5)", &options, mesh);
  mesh->communicate(N);

  Options state;
  state["species"]["h"]["density"] = N;

  // If D is not set, then the diffusion should not be calculated
  component.transform(state);

  ASSERT_FALSE(state["species"]["h"].isSet("density_source"));
  ASSERT_FALSE(state["species"]["h"].isSet("momentum_source"));
  ASSERT_FALSE(state["species"]["h"].isSet("energy_source"));
}

TEST_F(AnomalousDiffusionTest, ParticleDiffusion) {

  Coordinates* coords = mesh->getCoordinates();
  coords->setBxy(1.0); // Note: This is non-finite or zero?

  Options options;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;

  options["h"]["anomalous_D"] = 1.0; // Set particle diffusion for "h" species

  AnomalousDiffusion component("h", options, nullptr);

  Options state;
  state["species"]["h"]["density"] =
      FieldFactory::get()->create3D("1 + y * (x - 0.5)", &options, mesh);
  state["species"]["h"]["AA"] = 1.0; // Atomic mass number

  component.transform(state);

  // Expect all sources to be set
  ASSERT_TRUE(state["species"]["h"].isSet("density_source"));
  ASSERT_TRUE(state["species"]["h"].isSet("momentum_source"));
  ASSERT_TRUE(state["species"]["h"].isSet("energy_source"));

  // Expect momentum and energy sources to be zero
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["h"]["momentum_source"]), 0.0,
                           "RGN_NOBNDRY"));
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["h"]["energy_source"]), 0.0,
                           "RGN_NOBNDRY"));

  // Expect the sum over all cells of density source to be zero
  auto dV = coords->cell_volume(); // Cell volume

  Field3D source = get<Field3D>(state["species"]["h"]["density_source"]);
  BoutReal integral = 0.0;
  BOUT_FOR_SERIAL(i, source.getRegion("RGN_NOBNDRY")) {
    ASSERT_TRUE(std::isfinite(dV[i]))
        << "Volume element not finite at " << i.x() << ", " << i.y() << ", " << i.z();
    ASSERT_TRUE(std::isfinite(source[i]))
        << "Density source not finite at " << i.x() << ", " << i.y() << ", " << i.z();
    integral += source[i] * dV[i];
  }
  ASSERT_LT(abs(integral), 1e-3) << "Integral of density source should be close to zero";
}

namespace {
Options anomalousOptions(bool diagnose, bool D, bool chi) {
  Options options;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  if (D) {
    options["h"]["anomalous_D"] = 1.0;
  }
  if (chi) {
    options["h"]["anomalous_chi"] = 0.5;
  }
  options["h"]["diagnose"] = diagnose;
  return options;
}

Options anomalousState() {
  Options field_options;
  Options state;
  state["species"]["h"]["density"] =
      FieldFactory::get()->create3D("1 + y * (x - 0.5)", &field_options, mesh);
  state["species"]["h"]["temperature"] =
      FieldFactory::get()->create3D("2 + 0.3 * x * x + 0.1 * y", &field_options, mesh);
  state["species"]["h"]["velocity"] =
      FieldFactory::get()->create3D("0.2 * y", &field_options, mesh);
  state["species"]["h"]["AA"] = 1.0;
  return state;
}

Options anomalousOutputs() {
  return Options{{"Nnorm", 1.0}, {"Tnorm", 1.0}, {"Omega_ci", 1.0}, {"rho_s0", 1.0}};
}
} // namespace

TEST_F(AnomalousDiffusionTest, EnergyChannelsAreEnergySource) {
  mesh->getCoordinates()->setBxy(1.0);

  // each term on its own is the whole energy_source, both together sum to it
  for (const auto& [D, chi] : std::vector<std::pair<bool, bool>>{
           {true, false}, {false, true}, {true, true}}) {
    Options options = anomalousOptions(true, D, chi);
    AnomalousDiffusion component("h", options, nullptr);
    Options state = anomalousState();
    component.transform(state);
    Options outputs = anomalousOutputs();
    component.outputVars(outputs);

    const Field3D energy_source = get<Field3D>(state["species"]["h"]["energy_source"]);
    Field3D sum = zeroFrom(energy_source);
    EXPECT_EQ(D, outputs.isSet("Eh_anomalous_D"));
    EXPECT_EQ(chi, outputs.isSet("Eh_anomalous_chi"));
    for (const std::string ch : {"Eh_anomalous_D", "Eh_anomalous_chi"}) {
      if (outputs.isSet(ch)) {
        const Field3D E = get<Field3D>(outputs[ch]);
        EXPECT_GT(max(abs(E), false, "RGN_NOBNDRY"), 0.0) << ch;
        EXPECT_EQ("h", outputs[ch].attributes["species"].as<std::string>());
        sum += E;
      }
    }
    EXPECT_TRUE(IsFieldEqual(sum, energy_source, "RGN_NOBNDRY", 0.0))
        << "D = " << D << ", chi = " << chi;
  }
}

TEST_F(AnomalousDiffusionTest, EnergyChannelsAreCopies) {
  mesh->getCoordinates()->setBxy(1.0);

  // only D on, so energy_source is set from the D field and would share its data
  Options options = anomalousOptions(true, true, false);
  AnomalousDiffusion component("h", options, nullptr);
  Options state = anomalousState();
  component.transform(state);
  const Field3D expected = copy(get<Field3D>(state["species"]["h"]["energy_source"]));

  // a later component adding in place (like the sheath) mustn't change the channel
  Field3D energy_source = state["species"]["h"]["energy_source"].as<Field3D>();
  BOUT_FOR_SERIAL(i, energy_source.getRegion("RGN_ALL")) { energy_source[i] += 1.0; }

  Options outputs = anomalousOutputs();
  component.outputVars(outputs);
  EXPECT_TRUE(
      IsFieldEqual(get<Field3D>(outputs["Eh_anomalous_D"]), expected, "RGN_NOBNDRY", 0.0));
}

TEST_F(AnomalousDiffusionTest, DiagnoseDoesNotChangeState) {
  mesh->getCoordinates()->setBxy(1.0);

  Options options_off = anomalousOptions(false, true, true);
  Options options_on = anomalousOptions(true, true, true);
  AnomalousDiffusion component_off("h", options_off, nullptr);
  AnomalousDiffusion component_on("h", options_on, nullptr);
  Options state_off = anomalousState();
  Options state_on = anomalousState();
  component_off.transform(state_off);
  component_on.transform(state_on);

  Options outputs_off = anomalousOutputs();
  component_off.outputVars(outputs_off);
  EXPECT_FALSE(outputs_off.isSet("Eh_anomalous_D"));
  EXPECT_FALSE(outputs_off.isSet("Eh_anomalous_chi"));

  for (const std::string var :
       {"density_source", "momentum_source", "energy_source", "particle_flow_xlow",
        "particle_flow_ylow", "momentum_flow_xlow", "momentum_flow_ylow",
        "energy_flow_xlow", "energy_flow_ylow"}) {
    const Field3D on = get<Field3D>(state_on["species"]["h"][var]);
    const Field3D off = get<Field3D>(state_off["species"]["h"][var]);
    EXPECT_TRUE(IsFieldEqual(on, off, "RGN_NOBNDRY", 0.0)) << var;
  }
}
