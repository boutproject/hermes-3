#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh
#include "fake_mesh_fixture.hxx"

#include "../../include/evolve_pressure.hxx"
#include "../../include/evolve_energy.hxx"
#include <bout/solver.hxx>

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
template <typename T>
class BraginskiiConductionTest
  : public FakeMeshFixture {
public:
  BraginskiiConductionTest()
      : FakeMeshFixture(),
        options({{"units",
                  {{"eV", 1.0},
                   {"meters", 1.0},
                   {"seconds", 1.0},
                   {"inv_meters_cubed", 1e19}}},
                 {"test+", {{"thermal_conduction", true}, {"diagnose", true}}},
                 {"solver", {{"type", "euler"}}}}),
        output_base({{"Nnorm", 1e19}, {"Tnorm", 1.}, {"Omega_ci", 1.}, {"rho_s0", 1.}}),
        state_base({{"fastest_wave", 0.}, {"species", {{"test+", {{"velocity", 0.}, {"density", 1.}, {"pressure", 1.}, { "AA", 2. }}}}}}) {
    Options::root()["Ptest+"]["function"] = "1.0";
    temp1 = linearGradient(1., 0., 0., 1., 1., 0., 0.);
    temp2 = linearGradient(0., 0., 0., 1., 2., 0., 0.);
    temp_perp_variation = linearGradient(0., 1., 2., 0., 0., 0., 0.);
  }
  Options options, output_base, state_base;
  Field3D temp1, temp2, temp_perp_variation;

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

  T makeComponent(BoutReal kappa, std::string mode = "braginskii") {
    Options opts = options.copy();
    opts["test+"]["kappa_coefficient"] = kappa;
    opts["test+"]["conduction_collisions_mode"] = mode;
    // FIXME: The first time this is called it will result in a call
    // to MPI_Init, which is slow. Is there some way to fake/avoid
    // that?
    auto solver = Solver::create(&(options["solver"]));
    return T("test+", opts, solver.get());
  }

  Field3D getDeriv(const Options& state);
};

template <>
Field3D BraginskiiConductionTest<EvolveEnergy>::getDeriv(const Options& state) {
  return state["ddt(Etest+)"];
}

template <>
Field3D BraginskiiConductionTest<EvolvePressure>::getDeriv(const Options& state) {
  return state["ddt(Ptest+)"];
}

using EnergyTypes = ::testing::Types<EvolvePressure, EvolveEnergy>;
TYPED_TEST_SUITE(BraginskiiConductionTest, EnergyTypes);

TYPED_TEST(BraginskiiConductionTest, ConductionGradientScaling) {
  TypeParam component = this->makeComponent(1.);
  Options state0, state1 = this->state_base.copy(), state2;
  state1["species"]["test+"]["charge"] = 1;
  state1["species"]["test+"]["collision_frequencies"]["test+_test+_coll"] = 1.;
  state1["species"]["test+"]["collision_frequency"] = 1.;

  state0 = state1.copy();
  state2 = state1.copy();

  state0["species"]["test+"]["temperature"] = this->temp_perp_variation;
  state1["species"]["test+"]["temperature"] = this->temp1;
  state2["species"]["test+"]["temperature"] = this->temp2;

  Options output0 = this->output_base.copy(), output1 = this->output_base.copy(), output2 = this->output_base.copy();

  component.finally(state0);
  component.outputVars(output0);
  component.finally(state1);
  component.outputVars(output1);
  component.finally(state2);
  component.outputVars(output2);
  
  Field3D conduction0 = this->getDeriv(output0), conduction1 = this->getDeriv(output1), conduction2 = this->getDeriv(output2);
BOUT_FOR_SERIAL(i, conduction1.getRegion("RGN_NOBNDRY")) {
    // Conduction is proportional to the parallel second derivative of temperature
    EXPECT_NE(conduction1[i], 0.);
    EXPECT_FLOAT_EQ(conduction0[i], 0.);
    EXPECT_FLOAT_EQ(conduction2[i], 2*conduction1[i]);
  }
}

TYPED_TEST(BraginskiiConductionTest, ConductionKappaScaling) {
  TypeParam component1 = this->makeComponent(1.), component2 = this->makeComponent(2.);
  Options state = this->state_base.copy();
  state["species"]["test+"]["charge"] = 1;
  state["species"]["test+"]["collision_frequencies"]["test+_test+_coll"] = 1.;
  state["species"]["test+"]["collision_frequency"] = 1.;
  state["species"]["test+"]["temperature"] = this->temp1;

  Options  output1 = this->output_base.copy(), output2 = this->output_base.copy();

  component1.finally(state);
  component1.outputVars(output1);
  component2.finally(state);
  component2.outputVars(output2);
  
  Field3D conduction1 = this->getDeriv(output1), conduction2 = this->getDeriv(output2);
BOUT_FOR_SERIAL(i, conduction1.getRegion("RGN_NOBNDRY")) {
    // Conduction is proportional to the parallel diffusivity
    EXPECT_NE(conduction1[i], 0.);
    EXPECT_FLOAT_EQ(conduction2[i], 2*conduction1[i]);
  }
}


TYPED_TEST(BraginskiiConductionTest, ConductionCollisionScaling) {
  TypeParam component = this->makeComponent(1.);
  Options state0, state1 = this->state_base.copy(), state2;
  state1["species"]["test+"]["charge"] = 1;
  state1["species"]["test+"]["temperature"] = this->temp1;

  state0 = state1.copy();
  state2 = state1.copy();

  state0["species"]["test+"]["collision_frequencies"]["test+_test+_coll"] = 0.5;
  state0["species"]["test+"]["collision_frequency"] = 0.5;
  state1["species"]["test+"]["collision_frequencies"]["test+_test+_coll"] = 1.;
  state1["species"]["test+"]["collision_frequency"] = 1.;
  state2["species"]["test+"]["collision_frequencies"]["test+_test+_coll"] = 2.;
  state2["species"]["test+"]["collision_frequency"] = 2.;

  Options output0 = this->output_base.copy(), output1 = this->output_base.copy(), output2 = this->output_base.copy();

  component.finally(state0);
  component.outputVars(output0);
  component.finally(state1);
  component.outputVars(output1);
  component.finally(state2);
  component.outputVars(output2);
  
  Field3D conduction0 = this->getDeriv(output0), conduction1 = this->getDeriv(output1), conduction2 = this->getDeriv(output2);
BOUT_FOR_SERIAL(i, conduction1.getRegion("RGN_NOBNDRY")) {
    // Conduction is inversely proportional to the collision frequency
    EXPECT_NE(conduction1[i], 0.);
    EXPECT_FLOAT_EQ(conduction0[i], 2*conduction1[i]);
    EXPECT_FLOAT_EQ(2*conduction2[i], conduction1[i]);
  }
}


TYPED_TEST(BraginskiiConductionTest, ConductionCollisionsMode) {
  TypeParam component_braginskii = this->makeComponent(1., "braginskii"), component_multispecies = this->makeComponent(1., "multispecies");
  Options state = this->state_base.copy();
  state["species"]["test+"]["charge"] = 1;
  state["species"]["test+"]["temperature"] = this->temp1;
  state["species"]["test+"]["collision_frequencies"]["test+_test+_coll"] = 0.5;
  state["species"]["test+"]["collision_frequencies"]["test+_e_coll"] = 0.5;
  state["species"]["test+"]["collision_frequency"] = 1.0;

  Options output_brag = this->output_base.copy(), output_multi = this->output_base.copy();

  component_braginskii.finally(state);
  component_braginskii.outputVars(output_brag);
  component_multispecies.finally(state);
  component_multispecies.outputVars(output_multi);
  
  Field3D conduction_brag = this->getDeriv(output_brag), conduction_multi = this->getDeriv(output_multi);
  BOUT_FOR_SERIAL(i, conduction_brag.getRegion("RGN_NOBNDRY")) {
    // Conduction is inversely proportional to the collision frequency
    EXPECT_NE(conduction_brag[i], 0.);
    EXPECT_FLOAT_EQ(conduction_brag[i], 2*conduction_multi[i]);
  }
}
