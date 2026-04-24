#include <cmath>

#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/field3d.hxx>
#include <bout/options.hxx>
#include <bout/region.hxx>
#include "gtest/gtest.h"

#include "../../include/component.hxx"
#include "../../include/neutral_parallel_diffusion.hxx"
#include "fake_mesh_fixture.hxx"

#include <algorithm>

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using NeutralParallelDiffusionTest = FakeMeshFixture;

TEST_F(NeutralParallelDiffusionTest, ComponentRequiresDneut) {
  Options options; // No dneut setting
  EXPECT_THROW(const NeutralParallelDiffusion component("test", options, nullptr),
               BoutException);
}

TEST_F(NeutralParallelDiffusionTest, CreateComponent) {
  Options options{{"test", {{"dneut", 1.0}}}};
  const NeutralParallelDiffusion component("test", options, nullptr);
}

TEST_F(NeutralParallelDiffusionTest, NoNeutrals) {
  Options options{{"test", {{"dneut", 1.0}}}};
  NeutralParallelDiffusion component("test", options, nullptr);

  // State without any neutrals
  Options state{{"species",
                 {{"e", {{"AA", 1. / 1836}, {"charge", -1.0}}},
                  {"d+", {{"AA", 2.0}, {"charge", 1.0}}},
                  {"d2-", {{"AA", 2.0}, {"charge", -1.0}}}}}};

  component.declareAllSpecies(SpeciesInformation({"e"}, {}, {"d+"}, {"d2-"}));
  component.transform(state);

  ASSERT_FALSE(state["species"]["e"].isSet("density_source"));
  ASSERT_FALSE(state["species"]["d+"].isSet("density_source"));
  ASSERT_FALSE(state["species"]["d2-"].isSet("density_source"));
}

TEST_F(NeutralParallelDiffusionTest, NeutralNoCollisions) {
  Options options{{"test", {{"dneut", 1.0}}}};
  NeutralParallelDiffusion component("test", options, nullptr);

  // State without any collision frequencies
  Options state{
      {"species", {{"d", {{"AA", 2.0}, {"density", 1e19}, {"temperature", 1.0}}}}}};

  component.declareAllSpecies(SpeciesInformation({}, {"d"}, {}, {}));
  EXPECT_THROW(component.transform(state), BoutException);
}

TEST_F(NeutralParallelDiffusionTest, NeutralNoVelocity) {
  Options options{{"test", {{"dneut", 1.0}}}};
  NeutralParallelDiffusion component("test", options, nullptr);

  // State without any collision frequencies
  Options state{{"species",
                 {{"d",
                   {{"AA", 2.0},
                    {"density", 1e19},
                    {"temperature", 1.0},
                    {"collision_frequencies", {{"d_d+_iz", 1.0}}}}}}}};

  component.declareAllSpecies(SpeciesInformation({}, {"d"}, {}, {}));
  component.transform(state);

  ASSERT_TRUE(state["species"]["d"].isSet("density_source"));
  ASSERT_TRUE(state["species"]["d"].isSet("energy_source"));
  ASSERT_FALSE(state["species"]["d"].isSet("momentum_source"));
}

TEST_F(NeutralParallelDiffusionTest, NeutralWithVelocity) {
  Options options{{"test", {{"dneut", 1.0}}}};
  NeutralParallelDiffusion component("test", options, nullptr);

  // State without any collision frequencies
  Options state{{"species",
                 {{"d",
                   {{"AA", 2.0},
                    {"density", 1e19},
                    {"temperature", 1.0},
                    {"velocity", 0.0},
                    {"momentum", 0.0},
                    {"collision_frequencies", {{"d_d+_iz", 1.0}}}}}}}};

  component.declareAllSpecies(SpeciesInformation({}, {"d"}, {}, {}));
  component.transform(state);

  ASSERT_TRUE(state["species"]["d"].isSet("density_source"));
  ASSERT_TRUE(state["species"]["d"].isSet("energy_source"));
  ASSERT_TRUE(state["species"]["d"].isSet("momentum_source"));
}

TEST_F(NeutralParallelDiffusionTest, NeutralZeroCollisionFrequency) {
  Options options{{"test", {{"dneut", 1.0}}}};
  NeutralParallelDiffusion component("test", options, nullptr);

  // State without any collision frequencies
  Options state{{"species",
                 {{"d",
                   {{"AA", 2.0},
                    {"density", 1e19},
                    {"temperature", 1.0},
                    {"velocity", 0.0},
                    {"momentum", 0.0},
                    {"collision_frequencies", {{"d_d+_iz", 0.0}, {"d_d+_cx", 0.0}}}}}}}};

  component.declareAllSpecies(SpeciesInformation({}, {"d"}, {}, {}));
  component.transform(state);

  ASSERT_TRUE(state["species"]["d"].isSet("density_source"));
  ASSERT_TRUE(state["species"]["d"].isSet("energy_source"));
  ASSERT_TRUE(state["species"]["d"].isSet("momentum_source"));

  // Check that density source is finite
  const Field3D density_source = get<Field3D>(state["species"]["d"]["density_source"]);
  const auto& region = density_source.getRegion("RGN_NOBNDRY");
  ASSERT_TRUE(std::all_of(region.begin(), region.end(), [&](const auto& i) {
    return std::isfinite(density_source[i]);
  }));
}

TEST_F(NeutralParallelDiffusionTest, InvalidCollisionsMode) {
  Options options{
      {"test",
       {{"dneut", 1.0}, {"diffusion_collisions_mode", "multiecies"}}}}; // Mis-spelled

  NeutralParallelDiffusion component("test", options, nullptr);
  // State without any collision frequencies
  Options state{{"species",
                 {{"d",
                   {{"AA", 2.0},
                    {"density", 1e19},
                    {"temperature", 1.0},
                    {"velocity", 0.0},
                    {"momentum", 0.0},
                    {"collision_frequencies", {{"d_d+_iz", 0.0}, {"d_d+_cx", 0.0}}}}}}}};

  component.declareAllSpecies(SpeciesInformation({}, {"d"}, {}, {}));

  EXPECT_THROW(component.transform(state);, BoutException);
}

TEST_F(NeutralParallelDiffusionTest, Multispecies) {
  Options options{
      {"test", {{"dneut", 1.0}, {"diffusion_collisions_mode", "multispecies"}}}};
  NeutralParallelDiffusion component("test", options, nullptr);

  // State without any collision frequencies
  Options state{{"species",
                 {{"d",
                   {{"AA", 2.0},
                    {"density", 1e19},
                    {"temperature", 1.0},
                    {"velocity", 0.0},
                    {"momentum", 0.0},
                    {"collision_frequencies", {{"d_d+_iz", 0.0}, {"d_d+_cx", 0.0}}}}}}}};

  component.declareAllSpecies(SpeciesInformation({}, {"d"}, {}, {}));
  component.transform(state);

  ASSERT_TRUE(state["species"]["d"].isSet("density_source"));
  ASSERT_TRUE(state["species"]["d"].isSet("energy_source"));
  ASSERT_TRUE(state["species"]["d"].isSet("momentum_source"));

  // Check that density source is finite
  const Field3D density_source = get<Field3D>(state["species"]["d"]["density_source"]);
  const auto& region = density_source.getRegion("RGN_NOBNDRY");
  ASSERT_TRUE(std::all_of(region.begin(), region.end(), [&](const auto& i) {
    return std::isfinite(density_source[i]);
  }));
}

TEST_F(NeutralParallelDiffusionTest, Diagnose) {
  Options options{{"test", {{"dneut", 1.0}, {"diagnose", true}}}};

  NeutralParallelDiffusion component("test", options, nullptr);

  // State without any collision frequencies
  Options state{{"species",
                 {{"d",
                   {{"AA", 2.0},
                    {"density", 1e19},
                    {"temperature", 1.0},
                    {"velocity", 0.0},
                    {"momentum", 0.0},
                    {"collision_frequencies", {{"d_d+_iz", 0.0}, {"d_d+_cx", 0.0}}}}}}}};

  component.declareAllSpecies(SpeciesInformation({}, {"d"}, {}, {}));
  component.transform(state);

  Options outputs = {
      {"Tnorm", 1.0}, {"Nnorm", 1.0}, {"Omega_ci", 1.0}, {"Cs0", 1.0}, {"rho_s0", 1.0}};

  component.outputVars(outputs);

  ASSERT_TRUE(outputs.isSet("Dd_Dpar"));
  ASSERT_TRUE(outputs.isSet("Sd_Dpar"));
  ASSERT_TRUE(outputs.isSet("Ed_Dpar"));
}
