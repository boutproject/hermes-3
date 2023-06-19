
#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh

#include "../../include/zero_current.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

#include <bout/field_factory.hxx>  // For generating functions

// Reuse the "standard" fixture for FakeMesh
using ZeroCurrentTest = FakeMeshFixture;

TEST_F(ZeroCurrentTest, CreateComponent) {
  Options options;

  options["test"]["charge"] = 1.0; // Must be a charged species
  ZeroCurrent component("test", options, nullptr);
}

TEST_F(ZeroCurrentTest, ElectronFlowVelocity) {
  Options options;

  options["e"]["charge"] = -1.0;

  ZeroCurrent component("e", options, nullptr);

  options["species"]["e"]["charge"] = -1.0;
  options["species"]["e"]["density"] = 2.0;

  options["species"]["ion"]["density"] = 1.0;
  options["species"]["ion"]["charge"] = 2.0;

  // Set ion velocity
  Field3D Vi =  FieldFactory::get()->create3D("y - x", &options, mesh);
  options["species"]["ion"]["velocity"] = Vi;

  component.transform(options);

  // Electron velocity should be equal to ion velocity
  Field3D Ve = get<Field3D>(options["species"]["e"]["velocity"]);
  BOUT_FOR_SERIAL(i, Ve.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(Ve[i], Vi[i]) << "Electron velocity not equal to ion velocity at " << i;
  }
}
