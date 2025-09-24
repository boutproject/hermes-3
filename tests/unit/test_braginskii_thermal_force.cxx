
#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh
#include "fake_mesh_fixture.hxx"

#include "../../include/braginskii_thermal_force.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
class BraginskiiThermalForceTest
  : public FakeMeshFixture {
public:
  BraginskiiThermalForceTest()
      : FakeMeshFixture(), options({{"units",
                                     {{"eV", 1.0},
                                      {"meters", 1.0},
                                      {"seconds", 1.0},
                                      {"inv_meters_cubed", 1e19}}}}),
        component("test", options, nullptr) {
    grad1 = constantGradient(1., 0., 1., 0.);
    grad2 = constantGradient(0., 0., 2., 0.);
    grad_perp = constantGradient(0., 1., 0., 0.);
  }
  Options options;
  BraginskiiThermalForce component;
  Field3D grad1, grad2, grad_perp;

  static Field3D constantGradient(BoutReal a, BoutReal b_x, BoutReal b_y, BoutReal b_z) {
    Field3D result(a);
    BOUT_FOR_SERIAL(i, result.getRegion("RGN_ALL")) {
      result[i] += b_x * i.x() + b_y * i.y() + b_z * i.z();
    }
    return result;
  }
};

TEST_F(BraginskiiThermalForceTest, OnlyElectrons) {
  Options state;
  state["species"]["e"]["density"] = 1;
  state["species"]["e"]["temperature"] = grad1;
  state["species"]["e"]["charge"] = -1;
  state["species"]["e"]["AA"] = 1./1836;

  component.transform(state);
  EXPECT_FALSE(state["species"]["e"]["momentum_source"].isSet());
}

TEST_F(BraginskiiThermalForceTest, OnlyOneIon) {
  Options state;
  state["species"]["d+"]["density"] = 1;
  state["species"]["d+"]["temperature"] = grad1;
  state["species"]["d+"]["charge"] = 1;
  state["species"]["d+"]["AA"] = 2.;

  component.transform(state);
  EXPECT_FALSE(state["species"]["d+"]["momentum_source"].isSet());
}

TEST_F(BraginskiiThermalForceTest, ElectronIonBalance) {
  Options state;
  state["species"]["e"]["density"] = 1;
  state["species"]["e"]["temperature"] = grad1;
  state["species"]["e"]["charge"] = -1;
  state["species"]["e"]["AA"] = 1./1836;
  state["species"]["d+"]["density"] = 1;
  state["species"]["d+"]["temperature"] = grad1;
  state["species"]["d+"]["charge"] = 1;
  state["species"]["d+"]["AA"] = 2.;

  component.transform(state);
  
  Field3D mom_e = state["species"]["e"]["momentum_source"];
  Field3D mom_d = state["species"]["d+"]["momentum_source"];
  BOUT_FOR_SERIAL(i, mom_e.getRegion("RGN_NOBNDRY")) {
    // Forces on the species should be equal and opposite
    ASSERT_FLOAT_EQ(mom_e[i], -mom_d[i]);
  }
}

TEST_F(BraginskiiThermalForceTest, IonIonBalance) {
  Options state;
  state["species"]["c"]["density"] = 0.01;
  state["species"]["c"]["temperature"] = grad1;
  state["species"]["c"]["charge"] = 1;
  state["species"]["c"]["AA"] = 12.;
  state["species"]["d+"]["density"] = 1;
  state["species"]["d+"]["temperature"] = grad1;
  state["species"]["d+"]["charge"] = 1;
  state["species"]["d+"]["AA"] = 2.;

  component.transform(state);
  
  Field3D mom_c = state["species"]["c"]["momentum_source"];
  Field3D mom_d = state["species"]["d+"]["momentum_source"];
  BOUT_FOR_SERIAL(i, mom_c.getRegion("RGN_NOBNDRY")) {
    // Forces on the species should be equal and opposite
    ASSERT_FLOAT_EQ(mom_c[i], -mom_d[i]);
  }
}

TEST_F(BraginskiiThermalForceTest, NoNetForce) {
  Options state;
  state["species"]["c"]["density"] = 0.01;
  state["species"]["c"]["temperature"] = grad1;
  state["species"]["c"]["charge"] = 1;
  state["species"]["c"]["AA"] = 12.;
  state["species"]["ar"]["density"] = 0.01;
  state["species"]["ar"]["temperature"] = grad2;
  state["species"]["ar"]["charge"] = 1;
  state["species"]["ar"]["AA"] = 40.;
  state["species"]["d"]["density"] = 1;
  state["species"]["d"]["temperature"] = grad1;
  state["species"]["d"]["charge"] = 0;
  state["species"]["d"]["AA"] = 2.;
  state["species"]["d+"]["density"] = 1;
  state["species"]["d+"]["temperature"] = grad1;
  state["species"]["d+"]["charge"] = 1;
  state["species"]["d+"]["AA"] = 2.;
  state["species"]["e"]["density"] = 1;
  state["species"]["e"]["temperature"] = grad2;
  state["species"]["e"]["charge"] = -1;
  state["species"]["e"]["AA"] = 1./1836;

  component.transform(state);
  Field3D force(0.);
  for (const auto& [name, species] : state["species"].subsections()) {
    if ((*species)["charge"].as<int>() != 0) {
      force += (*species)["momentum_source"].as<Field3D>();
    }
  }
  // There is no external force on plasma, so all of the thermal
  // forces should balance out.
  BOUT_FOR_SERIAL(i, force.getRegion("RGN_NOBNDRY")) {
    ASSERT_FLOAT_EQ(force[i], 0.0);
  }
}

TEST_F(BraginskiiThermalForceTest, ElectronForceDensityScaling) {
  Options state;
  state["species"]["d1+"]["density"] = 1;
  state["species"]["d1+"]["temperature"] = grad1;
  state["species"]["d1+"]["charge"] = 1;
  state["species"]["d1+"]["AA"] = 2.;
  state["species"]["d2+"]["density"] = 2;
  state["species"]["d2+"]["temperature"] = grad1;
  state["species"]["d2+"]["charge"] = 1;
  state["species"]["d2+"]["AA"] = 2.;
  state["species"]["e"]["density"] = 1;
  state["species"]["e"]["temperature"] = grad1;
  state["species"]["e"]["charge"] = -1;
  state["species"]["e"]["AA"] = 1./1836;

  component.transform(state);
  Field3D mom1 = state["species"]["d1+"]["momentum_source"], mom2 = state["species"]["d2+"]["momentum_source"];
  BOUT_FOR_SERIAL(i, mom1.getRegion("RGN_NOBNDRY")) {
    // Force is directly proportional to ion density
    ASSERT_NE(mom1[i], 0.);
    ASSERT_FLOAT_EQ(2*mom1[i], mom2[i]);
  }
}

TEST_F(BraginskiiThermalForceTest, ElectronForceChargeScaling) {
  Options state;
  state["species"]["d1+"]["density"] = 1;
  state["species"]["d1+"]["temperature"] = grad1;
  state["species"]["d1+"]["charge"] = 1;
  state["species"]["d1+"]["AA"] = 2.;
  state["species"]["d2+"]["density"] = 1;
  state["species"]["d2+"]["temperature"] = grad1;
  state["species"]["d2+"]["charge"] = 2;
  state["species"]["d2+"]["AA"] = 2.;
  state["species"]["e"]["density"] = 1;
  state["species"]["e"]["temperature"] = grad1;
  state["species"]["e"]["charge"] = -1;
  state["species"]["e"]["AA"] = 1./1836;

  component.transform(state);
  Field3D mom1 = state["species"]["d1+"]["momentum_source"], mom2 = state["species"]["d2+"]["momentum_source"];
  BOUT_FOR_SERIAL(i, mom1.getRegion("RGN_NOBNDRY")) {
    // Force is proportional to square of ion density
    ASSERT_NE(mom1[i], 0.);
    ASSERT_FLOAT_EQ(4*mom1[i], mom2[i]);
  }
}

TEST_F(BraginskiiThermalForceTest, ElectronForceTemperatureGradScaling) {
  Options state0, state1, state2;
  state1["species"]["d+"]["density"] = 1;
  state1["species"]["d+"]["temperature"] = grad1;
  state1["species"]["d+"]["charge"] = 1;
  state1["species"]["d+"]["AA"] = 2.;
  state1["species"]["e"]["density"] = 1;
  state1["species"]["e"]["charge"] = -1;
  state1["species"]["e"]["AA"] = 1./1836;

  state0 = state1.copy();
  state2 = state1.copy();

  state0["species"]["e"]["temperature"] = grad_perp;
  state1["species"]["e"]["temperature"] = grad1;
  state2["species"]["e"]["temperature"] = grad2;

  state1["species"]["d2+"]["density"] = 1;
  state1["species"]["d2+"]["temperature"] = grad1;
  state1["species"]["d2+"]["charge"] = 1;
  state1["species"]["d2+"]["AA"] = 2.;  

  component.transform(state0);
  component.transform(state1);
  component.transform(state2);
  
  Field3D mom0 = state0["species"]["d+"]["momentum_source"], mom1 = state1["species"]["d+"]["momentum_source"], mom2 = state2["species"]["d+"]["momentum_source"], mom1_prime = state1["species"]["d2+"]["momentum_source"];
BOUT_FOR_SERIAL(i, mom1.getRegion("RGN_NOBNDRY")) {
    // Temperature gradient of the ion doesn't influence force
    ASSERT_FLOAT_EQ(mom1[i], mom1_prime[i]);
    // Force is proportional to the parallel electron temperature gradient
    ASSERT_NE(mom1[i], 0.);
    ASSERT_FLOAT_EQ(mom0[i], 0.);
    ASSERT_FLOAT_EQ(mom2[i], 2*mom1[i]);
  }
}

TEST_F(BraginskiiThermalForceTest, IonIonForceTemperatureGradScaling) {
  Options state1, state2;
  state1["species"]["c1"]["density"] = 0.01;
  state1["species"]["c1"]["temperature"] = grad1;
  state1["species"]["c1"]["charge"] = 1;
  state1["species"]["c1"]["AA"] = 12.;
  state1["species"]["c2"]["density"] = 0.01;
  state1["species"]["c2"]["temperature"] = grad2;
  state1["species"]["c2"]["charge"] = 1;
  state1["species"]["c2"]["AA"] = 12.;
  state1["species"]["d+"]["density"] = 1;
  state1["species"]["d+"]["charge"] = 1;
  state1["species"]["d+"]["AA"] = 2.;

  state2 = state1.copy();
  state1["species"]["d+"]["temperature"] = grad1;
  state2["species"]["d+"]["temperature"] = grad2;

  component.transform(state1);
  component.transform(state2);

  Field3D mom1_1 = state1["species"]["c1"]["momentum_source"], mom1_2 = state1["species"]["c2"]["momentum_source"], mom2_1 = state2["species"]["c1"]["momentum_source"], mom2_2 = state2["species"]["c2"]["momentum_source"];
  BOUT_FOR_SERIAL(i, mom1_1.getRegion("RGN_NOBNDRY")) {
    // Changing temperature gradient of the heavy ion should not change the force
    EXPECT_FLOAT_EQ(mom1_1[i], mom1_2[i]);
    EXPECT_FLOAT_EQ(mom2_1[i], mom2_2[i]);
    // The force is proportional to the temperature gradient of the light ion
    EXPECT_FLOAT_EQ(2*mom1_1[i], mom2_1[i]);
  }
}

class BraginskiiThermalForceTest_MassRatio
    : public BraginskiiThermalForceTest,
      public testing::WithParamInterface<std::tuple<int, int, bool>> {};

TEST_P(BraginskiiThermalForceTest_MassRatio, CheckForIonMasses) {
  auto [aa1, aa2, thermal_force_present] = GetParam();
  Options state{
      {"species",
       {{"M", {{"density", 1}, {"temperature", grad1}, {"charge", 1}, {"AA", aa1}}},
        {"N", {{"density", 1}, {"temperature", grad2}, {"charge", 2}, {"AA", aa2}}}}}};
  component.transform(state);
  if (thermal_force_present) {
    Field3D momentum_source1 = state["species"]["M"]["momentum_source"];
    Field3D momentum_source2 = state["species"]["N"]["momentum_source"];

    BOUT_FOR_SERIAL(i, momentum_source1.getRegion("RGN_NOBNDRY")) {
      // The masses of the ions are different enough for thermal force to be present.
      EXPECT_NE(momentum_source1[i], 0.);
      EXPECT_NE(momentum_source2[i], 0.);
    }
  } else {
    EXPECT_FALSE(state["species"]["M"]["momentum_source"].isSet());
    EXPECT_FALSE(state["species"]["N"]["momentum_source"].isSet());
  }
}

INSTANTIATE_TEST_SUITE_P(BraginskiiThermalIonMasses, BraginskiiThermalForceTest_MassRatio,
                         testing::Values(std::make_tuple(1, 50, true), std::make_tuple(50, 1, true), std::make_tuple(3, 11, true), std::make_tuple(11, 3, true), std::make_tuple(4, 10, false), std::make_tuple(10, 4, false), std::make_tuple(5, 50, false), std::make_tuple(50, 5, false), std::make_tuple(1, 7, false), std::make_tuple(7, 1, false)));
