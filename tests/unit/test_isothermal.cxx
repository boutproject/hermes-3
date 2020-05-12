
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

  Isothermal component("test", options, nullptr);
}

TEST_F(IsothermalTest, NeedsDensity) {
  Options options;
  options["units"]["eV"] = 5.0;

  // No temperature specified -> Uses normalised temperature of 1
  Isothermal component("test", options, nullptr);

  Options state;  // No density
  ASSERT_THROW(component.transform(state), BoutException);
}
       
TEST_F(IsothermalTest, DefaultTemperature) {
  Options options;
  options["units"]["eV"] = 5.0;

  // No temperature specified -> Uses normalised temperature of 1
  Isothermal component("e", options, nullptr);

  Field3D Ne = 2.0;
  
  Options state;
  state["species"]["e"]["density"] = Ne;

  component.transform(state);

  // Temperature should be a scalar, and normalised to 1
  auto Te = get<BoutReal>(state["species"]["e"]["temperature"]);
  ASSERT_DOUBLE_EQ(Te, 1.0);

  auto Pe = get<Field3D>(state["species"]["e"]["pressure"]);
  
  BOUT_FOR(i, Pe.getRegion("RGN_ALL")) {
    ASSERT_DOUBLE_EQ(Pe[i], 2.0);
  }
}

TEST_F(IsothermalTest, GivenTemperature) {
  Options options;
  options["units"]["eV"] = 5.0;
  options["e"]["temperature"] = 15.0; // In eV -> Normalised Te = 3
  
  // No temperature specified -> Uses normalised temperature of 1
  Isothermal component("e", options, nullptr);

  Field3D Ne = 2.0;
  
  Options state;
  state["species"]["e"]["density"] = Ne;

  component.transform(state);

  // Temperature should be a scalar, and normalised to 1
  auto Te = get<BoutReal>(state["species"]["e"]["temperature"]);
  ASSERT_DOUBLE_EQ(Te, 3.0);

  auto Pe = get<Field3D>(state["species"]["e"]["pressure"]);
  
  BOUT_FOR(i, Pe.getRegion("RGN_ALL")) {
    ASSERT_DOUBLE_EQ(Pe[i], 6.0);
  }
}
