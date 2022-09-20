#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh

#include "../../include/snb_conduction.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

#include <field_factory.hxx>  // For generating functions

// Reuse the "standard" fixture for FakeMesh
using SNBConductionTest = FakeMeshFixture;

TEST_F(SNBConductionTest, CreateComponent) {
  Options options;
  SNBConduction component("test", options, nullptr);
}

TEST_F(SNBConductionTest, Transform) {
  Options options;
  SNBConduction component("test", options, nullptr);

  Options state {{"units", {{"meters", 1.0},
                            {"eV", 1.0},
                            {"inv_meters_cubed", 1e19},
                            {"seconds", 1e-6}}},
                 {"species", {{"e", {{"temperature", 1.0},
                                     {"density", 1.0}}}}}};
  component.transform(state);

  ASSERT_TRUE(state["species"]["e"].isSet("energy_source"));
  // Zero temperature gradient everywhere -> No divergence of heat flux
  auto source = get<Field3D>(state["species"]["e"]["energy_source"]);
  BOUT_FOR_SERIAL(i, source.getRegion("RGN_NOBNDRY")) {
    ASSERT_LT(abs(source[i]), 1e-20);
  }
}
