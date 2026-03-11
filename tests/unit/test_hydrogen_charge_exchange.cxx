
#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/hydrogen_charge_exchange.hxx"

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
using HydrogenCXTest = FakeMeshFixture;

TEST_F(HydrogenCXTest, CreateComponent) {
  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"test", {{"type", "h + h+ -> h+ + h"}}}};

  HydrogenChargeExchange<'h', 'h'> component("test", options, nullptr);
}

TEST_F(HydrogenCXTest, RateAt1eV) {
  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"test", {{"type", "h + h+ -> h+ + h"}}}};

  HydrogenChargeExchange<'h', 'h'> component("test", options, nullptr);

  Field3D atom_dens(1.0), atom_temp(0.0), atom_vel(0.0);
  Field3D ion_dens(1.0), ion_temp(1.0), ion_vel(0.0);
  Options state{{"species",
                 {{"h",
                   {{"AA", 1.0},
                    {"density", atom_dens},
                    {"temperature", atom_temp},
                    {"velocity", atom_vel}}},
                  {"h+",
                   {{"AA", 1.0},
                    {"density", ion_dens},
                    {"temperature", ion_temp}, // lnT = 0.0
                    {"velocity", ion_vel}}}}},
                {"units", {{"eV", 1.0}}}};

  component.transform(state);

  // Should have set momentum and energy sources
  ASSERT_TRUE(state["species"]["h"].isSet("momentum_source"));
  ASSERT_TRUE(state["species"]["h"].isSet("energy_source"));
  ASSERT_TRUE(state["species"]["h+"].isSet("momentum_source"));
  ASSERT_TRUE(state["species"]["h+"].isSet("energy_source"));

  // Since here lnT = 0, we're only testing the b0 coefficient
  // <σv> should be exp(-18.5028) cm^3/s

  // Should appear as a source of energy in the atoms
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["h"]["energy_source"]),
                           (3. / 2) * exp(-18.5028) * 1e-6, "RGN_NOBNDRY"));
}
