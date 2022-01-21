
#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh

#include "../../include/amjuel_hyd_recombination.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

#include <bout/constants.hxx>
#include <field_factory.hxx> // For generating functions

// Reuse the "standard" fixture for FakeMesh
using HydrogenRCTest = FakeMeshFixture;

TEST_F(HydrogenRCTest, CreateComponent) {
  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}}};

  AmjuelHydRecombinationIsotope<'h'> component("test", options, nullptr);
}

// Check that recombination is a sink of ions, source of neutrals
TEST_F(HydrogenRCTest, DensitySourceSigns) {
  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}}};

  AmjuelHydRecombinationIsotope<'h'> component("test", options, nullptr);

  Options state{{"species",
                 {{"e",
                   {{"density", 1.0},
                    {"temperature", 1.0}}},
                  {"h",
                   {{"AA", 1.0},
                    {"density", 1.0},
                    {"temperature", 1.0},
                    {"velocity", 1.0}}},
                  {"h+",
                   {{"AA", 1.0},
                    {"charge", 1.0},
                    {"density", 1.0},
                    {"temperature", 1.0},
                    {"velocity", 1.0}}}}}};

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
    ASSERT_FLOAT_EQ(electron_density_source[i], ion_density_source[i]);

    ASSERT_TRUE(atom_momentum_source[i] > 0.0);
    ASSERT_TRUE(ion_momentum_source[i] < 0.0);

    ASSERT_TRUE(atom_energy_source[i] > 0.0);
    ASSERT_TRUE(ion_energy_source[i] < 0.0);
  }
}
