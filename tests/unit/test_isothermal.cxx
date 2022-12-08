
#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh

#include "../../include/isothermal.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using IsothermalTest = FakeMeshFixture;

TEST_F(IsothermalTest, CreateComponent) {
  Options options;
  options["units"]["eV"] = 1.0;
  options["test"]["temperature"] = 1.0;

  Isothermal component("test", options, nullptr);
}

TEST_F(IsothermalTest, NoDensity) {
  Options options;
  options["units"]["eV"] = 5.0; // Tnorm
  options["test"]["temperature"] = 1.0;  // 1eV

  Isothermal component("test", options, nullptr);

  Options state;  // No density
  component.transform(state);

  auto T = get<Field3D>(state["species"]["test"]["temperature"]);
  BOUT_FOR_SERIAL(i, T.getRegion("RGN_ALL")) {
    ASSERT_DOUBLE_EQ(T[i], 0.2);   // Normalised 1eV/Tnorm
  }

  // No pressure
  ASSERT_FALSE(state["species"]["test"]["pressure"].isSet());
}

TEST_F(IsothermalTest, GivenTandN) {
  Options options;
  options["units"]["eV"] = 5.0;
  options["e"]["temperature"] = 15.0; // In eV -> Normalised Te = 3

  Isothermal component("e", options, nullptr);

  Field3D Ne = 2.0;
  
  Options state;
  state["species"]["e"]["density"] = Ne;

  component.transform(state);

  // Temperature should be a scalar, and normalised to 1
  auto Te = get<Field3D>(state["species"]["e"]["temperature"]);
  BOUT_FOR_SERIAL(i, Te.getRegion("RGN_ALL")) {
    ASSERT_DOUBLE_EQ(Te[i], 3.0);
  }

  auto Pe = get<Field3D>(state["species"]["e"]["pressure"]);
  
  BOUT_FOR_SERIAL(i, Pe.getRegion("RGN_ALL")) {
    ASSERT_DOUBLE_EQ(Pe[i], 6.0);
  }
}
