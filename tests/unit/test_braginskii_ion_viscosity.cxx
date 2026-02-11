#include <bout/bout_types.hxx>
#include <bout/region.hxx>
#include <bout/utils.hxx>
#include "gtest/gtest.h"

#include "../../include/braginskii_ion_viscosity.hxx"
#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
class BraginskiiIonViscosityTest : public FakeMeshFixture {
public:
  BraginskiiIonViscosityTest()
      : FakeMeshFixture(), options({{"units",
                                     {{"eV", 1.0},
                                      {"meters", 1.0},
                                      {"seconds", 1.0},
                                      {"inv_meters_cubed", 1e19}}}}),
        component("test", options, nullptr) {}
  Options options;
  BraginskiiIonViscosity component;

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
      const int x = i.x();
      const int y = i.y();
      const int z = i.z();
      result[i] +=
          b1_x * x + b2_x * SQ(x) + b1_y * y + b2_y * SQ(y) + b1_z * z + b2_z * SQ(z);
    }
    return result;
  }
};

TEST_F(BraginskiiIonViscosityTest, ViscosityPressureScaling) {
  Options state1;
  Options state2;
  state1["species"]["d+"]["density"] = 1;
  state1["species"]["d+"]["temperature"] = 1;
  state1["species"]["d+"]["charge"] = 1;
  state1["species"]["d+"]["AA"] = 2.;
  state1["species"]["d+"]["velocity"] = linearGradient(1., 0., 1., 1., 1., 0., 0.);
  state1["species"]["d+"]["collision_frequencies"]["d+_d+_coll"] = 0.5;
  state1["species"]["d+"]["collision_frequencies"]["d+_c+_coll"] = 0.5;

  state2 = state1.copy();

  state1["species"]["d+"]["pressure"] = constantGradient(0, 0., 1., 0.);
  state2["species"]["d+"]["pressure"] =
      2 * state1["species"]["d+"]["pressure"].as<Field3D>();

  component.declareAllSpecies({"d+", "c+"});
  component.transform(state1);
  component.transform(state2);

  Field3D visc1 = state1["species"]["d+"]["momentum_source"];
  Field3D visc2 = state2["species"]["d+"]["momentum_source"];
  BOUT_FOR_SERIAL(i, visc1.getRegion("RGN_NOBNDRY")) {
    // Viscosity is proportional to the pressure.
    ASSERT_NE(visc1[i], 0.);
    ASSERT_DOUBLE_EQ(2 * visc1[i], visc2[i]);
  }
}

TEST_F(BraginskiiIonViscosityTest, ViscosityCollisionScaling) {
  Options state1;
  Options state2;
  state1["species"]["d+"]["density"] = 1;
  state1["species"]["d+"]["temperature"] = 1;
  state1["species"]["d+"]["charge"] = 1;
  state1["species"]["d+"]["AA"] = 2.;
  state1["species"]["d+"]["velocity"] = linearGradient(1., 0., 1., 1., 1., 0., 0.);
  state1["species"]["d+"]["pressure"] = constantGradient(0, 0., 1., 0.);

  state2 = state1.copy();
  state1["species"]["d+"]["collision_frequencies"]["d+_d+_coll"] =
      linearGradient(1., 0.1, 0.2, 0.1, 0.4, 1., 0.2);
  state1["species"]["d+"]["collision_frequencies"]["d+_he+_coll"] = 0.5;
  state2["species"]["d+"]["collision_frequencies"]["d+_d+_coll"] =
      2 * state1["species"]["d+"]["collision_frequencies"]["d+_d+_coll"].as<Field3D>();
  state2["species"]["d+"]["collision_frequencies"]["d+_he+_coll"] =
      2 * state1["species"]["d+"]["collision_frequencies"]["d+_he+_coll"].as<Field3D>();

  component.declareAllSpecies({"d+", "he+"});
  component.transform(state1);
  component.transform(state2);

  Field3D visc1 = state1["species"]["d+"]["momentum_source"];
  Field3D visc2 = state2["species"]["d+"]["momentum_source"];
  BOUT_FOR_SERIAL(i, visc1.getRegion("RGN_NOBNDRY")) {
    // Viscosity is inversely proportional to the collision frequency.
    ASSERT_NE(visc1[i], 0.);
    ASSERT_DOUBLE_EQ(visc1[i], 2 * visc2[i]);
  }
}

TEST_F(BraginskiiIonViscosityTest, ViscosityVelocityScaling) {
  Options state0;
  Options state1;
  Options state2;
  state1["species"]["d+"]["density"] = 1;
  state1["species"]["d+"]["temperature"] = 1;
  state1["species"]["d+"]["charge"] = 1;
  state1["species"]["d+"]["AA"] = 2.;
  state1["species"]["d+"]["collision_frequencies"]["d+_d+_coll"] = 0.5;
  state1["species"]["d+"]["collision_frequencies"]["d+_c+_coll"] = 0.5;
  state1["species"]["d+"]["pressure"] = constantGradient(0, 0., 1., 0.);

  state0 = state1.copy();
  state2 = state1.copy();

  state0["species"]["d+"]["velocity"] = constantGradient(1., 1., 0., 1.);
  state1["species"]["d+"]["velocity"] = linearGradient(1., 0., 1., 1., 1., 0., 0.);
  state2["species"]["d+"]["velocity"] =
      2 * state1["species"]["d+"]["velocity"].as<Field3D>();

  component.declareAllSpecies({"d+", "c+"});
  component.transform(state0);
  component.transform(state1);
  component.transform(state2);

  Field3D visc0 = state0["species"]["d+"]["momentum_source"];
  Field3D visc1 = state1["species"]["d+"]["momentum_source"];
  Field3D visc2 = state2["species"]["d+"]["momentum_source"];
  BOUT_FOR_SERIAL(i, visc1.getRegion("RGN_NOBNDRY")) {
    // There will be no viscosity if there is no parallel velocity gradient.
    ASSERT_DOUBLE_EQ(visc0[i], 0.);
    // Viscosity is proportional to the parallel gradient of the velocity.
    ASSERT_NE(visc1[i], 0.);
    ASSERT_DOUBLE_EQ(2 * visc1[i], visc2[i]);
  }
}

TEST_F(BraginskiiIonViscosityTest, ViscosityCollisionMode) {
  Options state1;
  Options state2;
  state1["species"]["d+"]["density"] = 1;
  state1["species"]["d+"]["temperature"] = 1;
  state1["species"]["d+"]["charge"] = 1;
  state1["species"]["d+"]["AA"] = 2.;
  state1["species"]["d+"]["velocity"] = linearGradient(1., 0., 1., 1., 1., 0., 0.);
  state1["species"]["d+"]["collision_frequencies"]["d+_d+_coll"] = 0.5;
  state1["species"]["d+"]["collision_frequencies"]["d+_c+_coll"] = 0.5;
  state1["species"]["d+"]["pressure"] = constantGradient(0, 0., 1., 0.);
  state2 = state1.copy();

  Options options2 = options.copy();
  options2["test2"]["viscosity_collisions_mode"] = "braginskii";
  BraginskiiIonViscosity component2("test2", options2, nullptr);

  component.declareAllSpecies({"d+", "c+"});
  component.transform(state1);
  component2.declareAllSpecies({"d+", "c+"});
  component2.transform(state2);

  Field3D visc1 = state1["species"]["d+"]["momentum_source"];
  Field3D visc2 = state2["species"]["d+"]["momentum_source"];
  BOUT_FOR_SERIAL(i, visc1.getRegion("RGN_NOBNDRY")) {
    // Collision frequency calculated using multispecies mode is twice that calculated
    // with Braginskii mode, so viscosity should be half.
    ASSERT_NE(visc1[i], 0.);
    ASSERT_DOUBLE_EQ(2 * visc1[i], visc2[i]);
  }
}

// TODO: Add tests for more types of collisions (e.g., electron-ion, charge-exchante, ion
// neutral...)

// TODO: Add tests for perpendicular flows
