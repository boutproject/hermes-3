#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx"

#include "../../include/classical_diffusion.hxx"

#include <bout/field_factory.hxx>

using bout::globals::mesh;

class ClassicalDiffusionTest : public FakeMeshFixture_tmpl<7, 5, 8, false, 2> {
protected:
  void SetUp() override {
    FakeMeshFixture_tmpl::SetUp();
    mesh->getCoordinates()->setBxy(2.0);
  }

  Options diagnosticState() {
    return Options{{"rho_s0", 2.0}, {"Omega_ci", 2.0}, {"Nnorm", 1.0}, {"Tnorm", 1.0}};
  }
};

TEST_F(ClassicalDiffusionTest, SpatialCustomDiffusionUsesSIUnits) {
  Options options;
  options["units"]["meters"] = 2.0;
  options["units"]["seconds"] = 0.5;
  options["classical_diffusion"]["custom_D"] = "8*(1+4*(x-0.5)^2)*(1+0.2*cos(z))";
  options["classical_diffusion"]["diagnose"] = true;
  const auto expected_D =
      FieldFactory::get()->create3D("(1+4*(x-0.5)^2)*(1+0.2*cos(z))", &options, mesh);
  const auto density = FieldFactory::get()->create3D("2+cos(z)", &options, mesh);

  ClassicalDiffusion component("classical_diffusion", options, nullptr);
  component.declareAllSpecies({"n"});
  Options initial_state;
  auto& species = initial_state["species"]["n"];
  species["density"] = density;
  species["pressure"] = 2.0 * density;
  species["velocity"] = Field3D{3.0};
  species["AA"] = 2.0;

  // Custom diffusion works for neutrals without charge or collision frequencies,
  // and the supplied coefficient survives repeated transforms.
  for (int repeat = 0; repeat < 2; ++repeat) {
    Options state = initial_state.copy();
    component.transform(state);
    auto output = diagnosticState();
    component.outputVars(output);
    const auto D = get<Field3D>(output["D_classical"]);
    EXPECT_TRUE(IsFieldEqual(D, expected_D, "RGN_NOBNDRY"));
    EXPECT_DOUBLE_EQ(output["D_classical"].attributes["conversion"].as<BoutReal>(), 8.0);

    const auto source = get<Field3D>(state["species"]["n"]["density_source"]);
    EXPECT_GT(max(abs(source), false, "RGN_NOBNDRY"), 0.0);
    EXPECT_NEAR(DomainBalance::integral(source), 0.0, 1e-12);
    EXPECT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["n"]["momentum_source"]),
                             Field3D{6.0 * source}, "RGN_NOBNDRY", 1e-12));
    EXPECT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["n"]["energy_source"]),
                             Field3D{3.0 * source}, "RGN_NOBNDRY", 1e-12));
  }
}

TEST_F(ClassicalDiffusionTest, ZeroCustomDiffusionNeedsNoCollisions) {
  Options options;
  options["units"]["meters"] = 2.0;
  options["units"]["seconds"] = 0.5;
  options["classical_diffusion"]["custom_D"] = 0.0;
  ClassicalDiffusion component("classical_diffusion", options, nullptr);
  component.declareAllSpecies({"n"});
  Options state;
  auto& species = state["species"]["n"];
  species["density"] = FieldFactory::get()->create3D("2+cos(z)", &options, mesh);
  species["temperature"] = FieldFactory::get()->create3D("2+sin(z)", &options, mesh);
  component.transform(state);

  EXPECT_TRUE(IsFieldEqual(get<Field3D>(species["density_source"]), 0.0, "RGN_NOBNDRY"));
  EXPECT_TRUE(IsFieldEqual(get<Field3D>(species["energy_source"]), 0.0, "RGN_NOBNDRY"));
}

TEST_F(ClassicalDiffusionTest, UnsetCustomDiffusionUsesCollisionsAndConduction) {
  Options options;
  options["classical_diffusion"]["diagnose"] = true;
  ClassicalDiffusion component("classical_diffusion", options, nullptr);
  component.declareAllSpecies({"e", "n"});
  const auto temperature = FieldFactory::get()->create3D("2+cos(z)", &options, mesh);
  Options state;
  auto& electrons = state["species"]["e"];
  electrons["charge"] = -1.0;
  electrons["AA"] = 0.5;
  electrons["density"] = Field3D{2.0};
  electrons["pressure"] = 2.0 * temperature;
  electrons["temperature"] = temperature;
  electrons["collision_frequency"] = Field3D{4.0};
  state["species"]["n"]["density"] = temperature;
  component.transform(state);
  auto output = diagnosticState();
  component.outputVars(output);

  EXPECT_TRUE(IsFieldEqual(get<Field3D>(output["D_classical"]),
                           Field3D{0.5 * temperature}, "RGN_NOBNDRY"));
  EXPECT_TRUE(
      IsFieldEqual(get<Field3D>(electrons["density_source"]), 0.0, "RGN_NOBNDRY"));
  const auto energy = get<Field3D>(electrons["energy_source"]);
  EXPECT_GT(max(abs(energy), false, "RGN_NOBNDRY"), 0.0);
  EXPECT_NEAR(DomainBalance::integral(energy), 0.0, 1e-12);
  EXPECT_FALSE(state["species"]["n"].isSet("density_source"));
}
