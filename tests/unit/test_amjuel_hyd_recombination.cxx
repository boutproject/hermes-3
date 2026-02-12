
#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/amjuel_hydrogen.hxx"

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
using HydrogenRCTest = FakeMeshFixture;

TEST_F(HydrogenRCTest, CreateComponent) {
  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"test", {{"type", "h+ + e -> h"}}}};

  AmjuelHydRecombinationIsotope<'h'> component("test", options, nullptr);
}

// Check that recombination is a sink of ions, source of neutrals
TEST_F(HydrogenRCTest, DensitySourceSigns) {
  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"test", {{"type", "h+ + e -> h"}}}};

  AmjuelHydRecombinationIsotope<'h'> component("test", options, nullptr);

  Field3D electron_dens(1.0), electron_temp(1.0), electron_vel(1.0);
  Field3D atom_dens(1.0), atom_temp(1.0), atom_vel(1.0);
  Field3D ion_dens(1.0), ion_temp(1.0), ion_vel(1.0);
  Options state{{"species",
                 {{"e",
                   {{"AA", 1.0},
                    {"density", electron_dens},
                    {"temperature", electron_temp},
                    {"velocity", electron_vel}}},
                  {"h",
                   {{"AA", 1.0},
                    {"density", atom_dens},
                    {"temperature", atom_temp},
                    {"velocity", atom_vel}}},
                  {"h+",
                   {{"AA", 1.0},
                    {"charge", 1.0},
                    {"density", ion_dens},
                    {"temperature", ion_temp},
                    {"velocity", ion_vel}}}}},
                {"test", {{"type", "h+ + e -> h"}}}};

  component.transform(state);

  ASSERT_TRUE(state["species"]["e"].isSet("energy_source"));

  auto atom_density_source = get<Field3D>(state["species"]["h"]["density_source"]);
  auto ion_density_source = get<Field3D>(state["species"]["h+"]["density_source"]);
  auto electron_density_source = get<Field3D>(state["species"]["e"]["density_source"]);

  auto atom_momentum_source = get<Field3D>(state["species"]["h"]["momentum_source"]);
  auto ion_momentum_source = get<Field3D>(state["species"]["h+"]["momentum_source"]);

  auto atom_energy_source = get<Field3D>(state["species"]["h"]["energy_source"]);
  auto ion_energy_source = get<Field3D>(state["species"]["h+"]["energy_source"]);

  BOUT_FOR_SERIAL(i, atom_density_source.getRegion("RGN_NOBNDRY")) {
    output.write("{}: {}\n", i.ind, atom_density_source[i]);
    ASSERT_TRUE(atom_density_source[i] > 0.0);
    ASSERT_TRUE(ion_density_source[i] < 0.0);
    ASSERT_DOUBLE_EQ(electron_density_source[i], ion_density_source[i]);

    ASSERT_TRUE(atom_momentum_source[i] > 0.0);
    ASSERT_TRUE(ion_momentum_source[i] < 0.0);

    ASSERT_TRUE(atom_energy_source[i] > 0.0);
    ASSERT_TRUE(ion_energy_source[i] < 0.0);
  }
}
