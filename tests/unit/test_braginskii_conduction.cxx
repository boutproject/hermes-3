#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include <bout/bout_types.hxx>
#include <bout/region.hxx>
#include <bout/solver.hxx>
#include <bout/utils.hxx>
#include "gtest/gtest.h"

#include "../../include/braginskii_conduction.hxx"
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
class BraginskiiConductionTest : public FakeMeshFixture {
public:
  BraginskiiConductionTest()
      : FakeMeshFixture(), options({
                               {"units",
                                {{"eV", 1.0},
                                 {"meters", 1.0},
                                 {"seconds", 1.0},
                                 {"inv_meters_cubed", 1e19}}},
                               {"test+",
                                {{"type", "evolve_energy"},
                                 {"thermal_conduction", true},
                                 {"diagnose", true}}},
                           }),
        output_base({{"Nnorm", 1e19}, {"Tnorm", 1.}, {"Omega_ci", 1.}, {"rho_s0", 1.}}),
        state_base(
            {{"fastest_wave", 0.},
             {"species",
              {{"test+",
                {{"velocity", 0.}, {"density", 1.}, {"pressure", 1.}, {"AA", 2.}}}}}}) {
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
      const int x = i.x();
      const int y = i.y();
      const int z = i.z();
      result[i] +=
          b1_x * x + b2_x * SQ(x) + b1_y * y + b2_y * SQ(y) + b1_z * z + b2_z * SQ(z);
    }
    return result;
  }

  BraginskiiConduction makeComponent(BoutReal kappa, std::string mode = "braginskii",
                                     BoutReal kappa_timescale_limit = -1.) {
    Options opts = options.copy();
    opts["test+"]["kappa_coefficient"] = kappa;
    opts["test+"]["conduction_collisions_mode"] = std::move(mode);
    opts["test+"]["kappa_timescale_limit"] = kappa_timescale_limit;
    return BraginskiiConduction("test+", opts, nullptr);
  }

  // Helper to make state for the conduction timescale limiter
  Options makeStateTimescaleLimiter(BoutReal density = 1.) {
    Options state = state_base.copy();
    state["species"]["test+"]["charge"] = 1;
    state["species"]["test+"]["density"].force(density);
    state["species"]["test+"]["collision_frequencies"]["test+_test+_coll"] = 1.;
    state["species"]["test+"]["collision_frequency"] = 1.;
    state["species"]["test+"]["temperature"] = temp1;
    return state;
  }

  Field3D getDeriv(const Options& state) {
    return state["species"]["test+"]["energy_source"];
  }
};

TEST_F(BraginskiiConductionTest, ConductionGradientScaling) {
  BraginskiiConduction component = this->makeComponent(1.);
  Options state0;
  Options state1 = this->state_base.copy();
  Options state2;
  state1["species"]["test+"]["charge"] = 1;
  state1["species"]["test+"]["collision_frequencies"]["test+_test+_coll"] = 1.;
  state1["species"]["test+"]["collision_frequency"] = 1.;

  state0 = state1.copy();
  state2 = state1.copy();

  state0["species"]["test+"]["temperature"] = this->temp_perp_variation;
  state1["species"]["test+"]["temperature"] = this->temp1;
  state2["species"]["test+"]["temperature"] = this->temp2;

  component.declareAllSpecies({"test+"});
  component.transform(state0);
  component.transform(state1);
  component.transform(state2);

  Field3D conduction0 = this->getDeriv(state0);
  Field3D conduction1 = this->getDeriv(state1);
  Field3D conduction2 = this->getDeriv(state2);
  BOUT_FOR_SERIAL(i, conduction1.getRegion("RGN_NOBNDRY")) {
    // Conduction is proportional to the parallel second derivative of temperature
    EXPECT_NE(conduction1[i], 0.);
    EXPECT_DOUBLE_EQ(conduction0[i], 0.);
    EXPECT_DOUBLE_EQ(conduction2[i], 2 * conduction1[i]);
  }
}

TEST_F(BraginskiiConductionTest, ConductionKappaScaling) {
  BraginskiiConduction component1 = this->makeComponent(1.);
  BraginskiiConduction component2 = this->makeComponent(2.);
  Options state1 = this->state_base.copy();
  state1["species"]["test+"]["charge"] = 1;
  state1["species"]["test+"]["collision_frequencies"]["test+_test+_coll"] = 1.;
  state1["species"]["test+"]["collision_frequency"] = 1.;
  state1["species"]["test+"]["temperature"] = this->temp1;
  Options state2 = state1.copy();

  component1.declareAllSpecies({"test+"});
  component1.transform(state1);
  component2.declareAllSpecies({"test+"});
  component2.transform(state2);

  Field3D conduction1 = this->getDeriv(state1);
  Field3D conduction2 = this->getDeriv(state2);
  BOUT_FOR_SERIAL(i, conduction1.getRegion("RGN_NOBNDRY")) {
    // Conduction is proportional to the parallel diffusivity
    EXPECT_NE(conduction1[i], 0.);
    EXPECT_DOUBLE_EQ(conduction2[i], 2 * conduction1[i]);
  }
}

TEST_F(BraginskiiConductionTest, ConductionCollisionScaling) {
  BraginskiiConduction component = this->makeComponent(1.);
  Options state0;
  Options state1 = this->state_base.copy();
  Options state2;
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

  component.declareAllSpecies({"test+"});
  component.transform(state0);
  component.transform(state1);
  component.transform(state2);

  Field3D conduction0 = this->getDeriv(state0);
  Field3D conduction1 = this->getDeriv(state1);
  Field3D conduction2 = this->getDeriv(state2);
  BOUT_FOR_SERIAL(i, conduction1.getRegion("RGN_NOBNDRY")) {
    // Conduction is inversely proportional to the collision frequency
    EXPECT_NE(conduction1[i], 0.);
    EXPECT_DOUBLE_EQ(conduction0[i], 2 * conduction1[i]);
    EXPECT_DOUBLE_EQ(2 * conduction2[i], conduction1[i]);
  }
}

TEST_F(BraginskiiConductionTest, ConductionCollisionsMode) {
  BraginskiiConduction component_braginskii = this->makeComponent(1., "braginskii");
  BraginskiiConduction component_multispecies = this->makeComponent(1., "multispecies");
  Options state_brag = this->state_base.copy();
  state_brag["species"]["test+"]["charge"] = 1;
  state_brag["species"]["test+"]["temperature"] = this->temp1;
  state_brag["species"]["test+"]["collision_frequencies"]["test+_test+_coll"] = 0.5;
  state_brag["species"]["test+"]["collision_frequencies"]["test+_e_coll"] = 0.5;
  state_brag["species"]["test+"]["collision_frequency"] = 1.0;
  Options state_multi = state_brag.copy();

  component_braginskii.declareAllSpecies({"test+"});
  component_braginskii.transform(state_brag);
  component_multispecies.declareAllSpecies({"test+", "e"});
  component_multispecies.transform(state_multi);

  Field3D conduction_brag = this->getDeriv(state_brag);
  Field3D conduction_multi = this->getDeriv(state_multi);
  BOUT_FOR_SERIAL(i, conduction_brag.getRegion("RGN_NOBNDRY")) {
    EXPECT_NE(conduction_brag[i], 0.);
    EXPECT_DOUBLE_EQ(conduction_brag[i], 2 * conduction_multi[i]);
  }
}

// Test timescale limiter. Fake mesh is closed field lines only.
TEST_F(BraginskiiConductionTest, KappaTimescaleLimitFormula) {
  // The density enters only the limiter, so varying it tests the limiter alone
  for (const auto& [tau_limit, density] :
       std::vector<std::pair<BoutReal, BoutReal>>{{1., 1.}, {2., 1.}, {1., 2.}}) {
    BraginskiiConduction component = this->makeComponent(1., "braginskii", tau_limit);
    Options state = this->makeStateTimescaleLimiter(density);

    component.declareAllSpecies({"test+"});
    component.transform(state);
    Options outputs = this->output_base.copy();
    component.outputVars(outputs);

    ASSERT_TRUE(outputs.isSet("kappa_par_test+_timelim_max"));
    ASSERT_TRUE(outputs.isSet("kappa_par_test+_unlimited"));

    // kappa_par / kappa_par_unlimited = 1 / (1 + kappa_par_unlimited / kappa_par_max)
    //                                 = 1 / (1 + 0.5 / (n / tau_limit))
    //                                 = 1 / (1 + 0.5 * tau_limit / n)
    const BoutReal limit_factor = 1. / (1. + 0.5 * tau_limit / density);
    Field3D kappa_par_max = outputs["kappa_par_test+_timelim_max"];
    Field3D kappa_par_unlimited = outputs["kappa_par_test+_unlimited"];
    Field3D kappa_par = state["species"]["test+"]["kappa_par"];
    BOUT_FOR_SERIAL(i, kappa_par.getRegion("RGN_NOBNDRY")) {
      // kappa_par_max = n * dl^2 / tau_limit
      //               = n * (dy * sqrt(g_22))^2 / tau_limit
      //               = n * (1 * sqrt(1))^2 / tau_limit
      //               = n / tau_limit
      EXPECT_DOUBLE_EQ(kappa_par_max[i], density / tau_limit);
      // kappa_par_unlimited = kappa_coefficient * P * (1 / nu) / AA
      //                     = 1 * 1 * (1 / 1) / 2
      //                     = 0.5
      EXPECT_DOUBLE_EQ(kappa_par_unlimited[i], 0.5);
      // kappa_par = kappa_par_unlimited / (1 + kappa_par_unlimited / kappa_par_max)
      //           = 0.5 / (1 + 0.5 / (n / tau_limit))
      //           = 0.5 * limit_factor
      EXPECT_DOUBLE_EQ(kappa_par[i], 0.5 * limit_factor);
      // Harmonic combination: 1 / kappa_par = 1 / kappa_par_unlimited + 1 / kappa_par_max
      EXPECT_DOUBLE_EQ(kappa_par[i], kappa_par_unlimited[i] * kappa_par_max[i]
                                         / (kappa_par_unlimited[i] + kappa_par_max[i]));
    }
  }
}

// Test that timescale limiter is inactive on open field lines
TEST_F(BraginskiiConductionTest, KappaTimescaleLimitOpenFieldLines) {
  static_cast<FakeMesh*>(mesh)->ix_separatrix = 0; // Open field lines
  BraginskiiConduction component_off = this->makeComponent(1.);
  BraginskiiConduction component_limit = this->makeComponent(1., "braginskii", 1.);
  Options state_off = this->makeStateTimescaleLimiter();
  Options state_limit = this->makeStateTimescaleLimiter();

  component_off.declareAllSpecies({"test+"});
  component_off.transform(state_off);
  component_limit.declareAllSpecies({"test+"});
  component_limit.transform(state_limit);
  Options outputs = this->output_base.copy();
  component_limit.outputVars(outputs);

  Field3D kappa_par_max = outputs["kappa_par_test+_timelim_max"];
  Field3D kappa_par_unlimited = outputs["kappa_par_test+_unlimited"];
  Field3D kappa_par = state_limit["species"]["test+"]["kappa_par"];
  Field3D conduction_off = this->getDeriv(state_off);
  Field3D conduction_limit = this->getDeriv(state_limit);
  BOUT_FOR_SERIAL(i, kappa_par.getRegion("RGN_NOBNDRY")) {
    EXPECT_TRUE(std::isnan(kappa_par_max[i]));
    EXPECT_DOUBLE_EQ(kappa_par[i], kappa_par_unlimited[i]);
    EXPECT_NE(conduction_off[i], 0.);
    EXPECT_DOUBLE_EQ(conduction_limit[i], conduction_off[i]);
  }
}

// TODO: Add tests for more types of collisions (e.g., electron-ion, charge-exchante, ion
// neutral...)
