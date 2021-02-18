
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

#include <field_factory.hxx>  // For generating functions

// Reuse the "standard" fixture for FakeMesh
using ZeroCurrentTest = FakeMeshFixture;

TEST_F(ZeroCurrentTest, CreateComponent) {
  Options options;
  
  ZeroCurrent component("test", options, nullptr);
}

TEST_F(ZeroCurrentTest, MissingElectronPressure) {
  Options options;
  
  ZeroCurrent component("test", options, nullptr);

  ASSERT_THROW(component.transform(options), BoutException);
}

TEST_F(ZeroCurrentTest, ZeroPressureGradient) {
  Options options;
  
  ZeroCurrent component("test", options, nullptr);

  options["species"]["e"]["pressure"] = 1.0;
  options["species"]["e"]["density"] = 1.0;
  options["species"]["e"]["charge"] = -1.0;

  options["species"]["h+"]["density"] = 1.0;
  options["species"]["h+"]["charge"] = 1.0;

  component.transform(options);
  
  // Should have a momentum source, but zero because no pressure gradient
  ASSERT_TRUE(options["species"]["h+"].isSet("momentum_source"));
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(options["species"]["h+"]["momentum_source"]), 0.0,
                           "RGN_NOBNDRY"));
}

TEST_F(ZeroCurrentTest, WithPressureGradient) {
  Options options;
  
  ZeroCurrent component("test", options, nullptr);

  // Create a function which increases in y
  options["species"]["e"]["pressure"] =
      FieldFactory::get()->create3D("y", &options, mesh);
  options["species"]["e"]["density"] = 1.0;
  options["species"]["e"]["charge"] = -1.0;

  options["species"]["h+"]["density"] = 1.0;
  options["species"]["h+"]["charge"] = 1.0;

  component.transform(options);
  
  // Should have a momentum source
  ASSERT_TRUE(options["species"]["h+"].isSet("momentum_source"));

  // Force on ions should be in negative y direction
  Field3D F = get<Field3D>(options["species"]["h+"]["momentum_source"]);

  BOUT_FOR_SERIAL(i, F.getRegion("RGN_NOBNDRY")) {
    ASSERT_LT(F[i], 0.0) << "Force on ions (h+) not negative at " << i;
  }
}

TEST_F(ZeroCurrentTest, ElectronFlowVelocity) {
  Options options;
  
  ZeroCurrent component("test", options, nullptr);

  options["species"]["e"]["pressure"] =
      FieldFactory::get()->create3D("y", &options, mesh);
  options["species"]["e"]["density"] = 2.0;
  options["species"]["e"]["charge"] = -1.0;

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

TEST_F(ZeroCurrentTest, ElectronForceBalance) {
  Options options;
  ZeroCurrent component("test", options, nullptr);

  options["species"]["e"]["pressure"] =1.0;
  options["species"]["e"]["density"] = 2.0;
  options["species"]["e"]["charge"] = -1.0;

  options["species"]["e"]["momentum_source"] = 0.5;

  // Should have E = momentum_source / density = 0.5 / 2.0

  options["species"]["ion"]["density"] = 1.0;
  options["species"]["ion"]["charge"] = 3.0;

  component.transform(options);

  // Should give ion momentum source charge * E = 3 * 0.5 / 2.0
  Field3D F = get<Field3D>(options["species"]["ion"]["momentum_source"]);

  BOUT_FOR_SERIAL(i, F.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(F[i], 3. * 0.5 / 2.0) << "Momentum source not correct at " << i;
  }
}
