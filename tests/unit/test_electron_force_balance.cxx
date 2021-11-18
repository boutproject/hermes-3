
#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh

#include "../../include/electron_force_balance.hxx"

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
using ElectronForceBalanceTest = FakeMeshFixture;

TEST_F(ElectronForceBalanceTest, CreateComponent) {
  Options options;

  ElectronForceBalance component("test", options, nullptr);
}

TEST_F(ElectronForceBalanceTest, MissingElectronPressure) {
  Options options;

  ElectronForceBalance component("test", options, nullptr);

  ASSERT_THROW(component.transform(options), BoutException);
}

TEST_F(ElectronForceBalanceTest, ZeroPressureGradient) {
  Options options;

  ElectronForceBalance component("test", options, nullptr);

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

TEST_F(ElectronForceBalanceTest, WithPressureGradient) {
  Options options;

  ElectronForceBalance component("test", options, nullptr);

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


TEST_F(ElectronForceBalanceTest, ForceBalance) {
  Options options;
  ElectronForceBalance component("test", options, nullptr);

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
