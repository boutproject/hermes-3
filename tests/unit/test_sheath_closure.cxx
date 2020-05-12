
#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh

#include "../../include/sheath_closure.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using SheathClosureTest = FakeMeshFixture;

TEST_F(SheathClosureTest, CreateComponent) {
  Options options;
  options["units"]["meters"] = 1.0;
  options["test"]["connection_length"] = 10;
  
  SheathClosure component("test", options, nullptr);
}

TEST_F(SheathClosureTest, NeedsDensity) {
  Options options;
  options["units"]["meters"] = 1.0;
  options["test"]["connection_length"] = 10;
  
  SheathClosure component("test", options, nullptr);

  Options state;
  state["fields"]["phi"] = Field3D(2.0);

  // Needs electron density
  ASSERT_THROW(component.transform(state), BoutException);
}

TEST_F(SheathClosureTest, PhiAndDensity) {
  Options options;
  options["units"]["meters"] = 1.0;
  options["test"]["connection_length"] = 10;
  
  SheathClosure component("test", options, nullptr);

  Options state;
  state["fields"]["phi"] = Field3D(2.0);
  state["species"]["e"]["density"] = Field3D(1.5);
  component.transform(state);

  ASSERT_TRUE(state["fields"].isSet("DivJextra"));
  ASSERT_TRUE(state["species"]["e"].isSet("density_source"));
}

TEST_F(SheathClosureTest, Temperature) {
  Options options;
  options["units"]["meters"] = 1.0;
  options["test"]["connection_length"] = 10;
  
  SheathClosure component("test", options, nullptr);

  Options state;
  state["fields"]["phi"] = Field3D(2.0);
  state["species"]["e"]["density"] = Field3D(1.5);
  state["species"]["e"]["temperature"] = Field3D(1.2);
  
  component.transform(state);

  ASSERT_TRUE(state["fields"].isSet("DivJextra"));
  ASSERT_TRUE(state["species"]["e"].isSet("density_source"));
  ASSERT_TRUE(state["species"]["e"].isSet("energy_source"));

  Field3D energy_source = get<Field3D>(state["species"]["e"]["energy_source"]);
  BOUT_FOR(i, energy_source.getRegion("RGN_NOBNDRY")) {
    ASSERT_TRUE(energy_source[i] <= 0.0); // Always a sink
  }
}
