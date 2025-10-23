
#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/braginskii_electron_viscosity.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
class BraginskiiElectronViscosityTest : public FakeMeshFixture {
public:
  BraginskiiElectronViscosityTest()
      : FakeMeshFixture(), options({{"units",
                                     {{"eV", 1.0},
                                      {"meters", 1.0},
                                      {"seconds", 1.0},
                                      {"inv_meters_cubed", 1e19}}}}),
        component("test", options, nullptr) {}
  Options options;
  BraginskiiElectronViscosity component;

  static Field3D constantGradient(BoutReal a, BoutReal b_x, BoutReal b_y, BoutReal b_z) {
    Field3D result(a);
    BOUT_FOR_SERIAL(i, result.getRegion("RGN_ALL")) {
      result[i] += b_x * i.x() + b_y * i.y() + b_z * i.z();
    }
    return result;
  }

  static Field3D linearGradient(BoutReal a, BoutReal b1_x, BoutReal b2_x, BoutReal b1_y,
                                BoutReal b2_y, BoutReal b1_z, BoutReal b2_z) {
    Field3D result(a);
    BOUT_FOR_SERIAL(i, result.getRegion("RGN_ALL")) {
      const int x = i.x(), y = i.y(), z = i.z();
      result[i] +=
          b1_x * x + b2_x * SQ(x) + b1_y * y + b2_y * SQ(y) + b1_z * z + b2_z * SQ(z);
    }
    return result;
  }
};

TEST_F(BraginskiiElectronViscosityTest, ViscosityPressureScaling) {
  Options state1, state2;
  state1["species"]["e"]["density"] = 1;
  state1["species"]["e"]["charge"] = -1;
  state1["species"]["e"]["AA"] = 1. / 1836;
  state1["species"]["e"]["velocity"] = linearGradient(1., 0., 1., 1., 1., 0., 0.);
  state1["species"]["e"]["collision_frequency"] = 1.;

  state2 = state1.copy();

  state1["species"]["e"]["pressure"] = constantGradient(0, 0., 1., 0.);
  state2["species"]["e"]["pressure"] =
      2 * state1["species"]["e"]["pressure"].as<Field3D>();

  component.transform(state1);
  component.transform(state2);

  Field3D visc1 = state1["species"]["e"]["momentum_source"],
          visc2 = state2["species"]["e"]["momentum_source"];
  BOUT_FOR_SERIAL(i, visc1.getRegion("RGN_NOBNDRY")) {
    // Viscosity is proportional to the pressure.
    ASSERT_NE(visc1[i], 0.);
    ASSERT_FLOAT_EQ(2 * visc1[i], visc2[i]);
  }
}

TEST_F(BraginskiiElectronViscosityTest, ViscosityCollisionScaling) {
  Options state1, state2;
  state1["species"]["e"]["density"] = 1;
  state1["species"]["e"]["charge"] = -1;
  state1["species"]["e"]["AA"] = 1. / 1836;
  state1["species"]["e"]["velocity"] = linearGradient(1., 0., 1., 1., 1., 0., 0.);
  state1["species"]["e"]["pressure"] = constantGradient(0, 0., 1., 0.);

  state2 = state1.copy();

  state1["species"]["e"]["collision_frequency"] =
      linearGradient(1., 0.1, 0.2, -0.1, -0.4, 0., 1.);
  state2["species"]["e"]["collision_frequency"] =
      2 * state1["species"]["e"]["collision_frequency"].as<Field3D>();

  component.transform(state1);
  component.transform(state2);

  Field3D visc1 = state1["species"]["e"]["momentum_source"],
          visc2 = state2["species"]["e"]["momentum_source"];
  BOUT_FOR_SERIAL(i, visc1.getRegion("RGN_NOBNDRY")) {
    // Viscosity is inversely proportional to the collision frequency.
    ASSERT_NE(visc1[i], 0.);
    ASSERT_FLOAT_EQ(visc1[i], 2 * visc2[i]);
  }
}

TEST_F(BraginskiiElectronViscosityTest, ViscosityVelocityScaling) {
  Options state0, state1, state2;
  state1["species"]["e"]["density"] = 1;
  state1["species"]["e"]["charge"] = -1;
  state1["species"]["e"]["AA"] = 1. / 1836;
  state1["species"]["e"]["collision_frequency"] = 1.;
  state1["species"]["e"]["pressure"] = constantGradient(0, 0., 1., 0.);

  state0 = state1.copy();
  state2 = state1.copy();

  state0["species"]["e"]["velocity"] = constantGradient(1., 1., 0., 1.);
  state1["species"]["e"]["velocity"] = linearGradient(1., 0., 1., 1., 1., 0., 0.);
  state2["species"]["e"]["velocity"] =
      2 * state1["species"]["e"]["velocity"].as<Field3D>();

  component.transform(state0);
  component.transform(state1);
  component.transform(state2);

  Field3D visc0 = state0["species"]["e"]["momentum_source"],
          visc1 = state1["species"]["e"]["momentum_source"],
          visc2 = state2["species"]["e"]["momentum_source"];
  BOUT_FOR_SERIAL(i, visc1.getRegion("RGN_NOBNDRY")) {
    // There will be no viscosity if there is no parallel velocity gradient.
    ASSERT_FLOAT_EQ(visc0[i], 0.);
    // Viscosity is proportional to the parallel gradient of the velocity.
    ASSERT_NE(visc1[i], 0.);
    ASSERT_FLOAT_EQ(2 * visc1[i], visc2[i]);
  }
}
