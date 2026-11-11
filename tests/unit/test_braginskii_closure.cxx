
#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/braginskii_closure.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using BraginskiiClosureTest = FakeMeshFixture;

std::set<ComponentInformation> makeExpected(std::initializer_list<std::string> names) {
  std::set<ComponentInformation> result;
  for (const auto& name : names) {
    result.emplace(name, name);
  }
  return result;
}

std::set<ComponentInformation> toSet(std::vector<ComponentInformation> input) {
  return std::set(input.begin(), input.end());
}

TEST_F(BraginskiiClosureTest, CreateDefault) {
  Options options;
  BraginskiiClosure component("test", options, nullptr);
  EXPECT_EQ(toSet(component.additionalComponents()),
            makeExpected({"braginskii_collisions", "braginskii_friction",
                          "braginskii_heat_exchange", "braginskii_conduction",
                          "braginskii_electron_viscosity", "braginskii_ion_viscosity",
                          "braginskii_thermal_force"}));
}

TEST_F(BraginskiiClosureTest, CreateMinimal) {
  Options options = {{"test",
                      {{"electron_viscosity", false},
                       {"ion_viscosity", false},
                       {"thermal_force", false}}}};
  BraginskiiClosure component("test", options, nullptr);
  EXPECT_EQ(toSet(component.additionalComponents()),
            makeExpected({"braginskii_collisions", "braginskii_friction",
                          "braginskii_heat_exchange", "braginskii_conduction"}));
}

TEST_F(BraginskiiClosureTest, CreateThermalForce) {
  Options options = {{"test",
                      {{"electron_viscosity", false},
                       {"ion_viscosity", false},
                       {"thermal_force", true}}}};
  BraginskiiClosure component("test", options, nullptr);
  EXPECT_EQ(toSet(component.additionalComponents()),
            makeExpected({"braginskii_collisions", "braginskii_friction",
                          "braginskii_heat_exchange", "braginskii_conduction",
                          "braginskii_thermal_force"}));
}

TEST_F(BraginskiiClosureTest, CreateIonViscosity) {
  Options options = {{"test",
                      {{"electron_viscosity", false},
                       {"ion_viscosity", true},
                       {"thermal_force", false}}}};
  BraginskiiClosure component("test", options, nullptr);
  EXPECT_EQ(toSet(component.additionalComponents()),
            makeExpected({"braginskii_collisions", "braginskii_friction",
                          "braginskii_heat_exchange", "braginskii_conduction",
                          "braginskii_ion_viscosity"}));
}

TEST_F(BraginskiiClosureTest, CreateElectronViscosity) {
  Options options = {{"test",
                      {{"electron_viscosity", true},
                       {"ion_viscosity", false},
                       {"thermal_force", false}}}};
  BraginskiiClosure component("test", options, nullptr);
  EXPECT_EQ(toSet(component.additionalComponents()),
            makeExpected({"braginskii_collisions", "braginskii_friction",
                          "braginskii_heat_exchange", "braginskii_conduction",
                          "braginskii_electron_viscosity"}));
}
