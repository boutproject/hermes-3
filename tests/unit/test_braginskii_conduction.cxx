#include <cmath>
#include <limits>
#include <string>
#include <utility>

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

  BraginskiiConduction makeComponent(BoutReal kappa, std::string mode = "braginskii") {
    Options opts = options.copy();
    opts["test+"]["kappa_coefficient"] = kappa;
    opts["test+"]["conduction_collisions_mode"] = std::move(mode);
    return BraginskiiConduction("test+", opts, nullptr);
  }

  Field3D getDeriv(const Options& state) {
    return state["species"]["test+"]["energy_source"];
  }

  Options uniformState(BoutReal nu = 1.0, BoutReal temperature = 2.0,
                       BoutReal density = 1.0) {
    auto state = state_base.copy();
    auto& species = state["species"]["test+"];
    species["charge"] = 1;
    species["temperature"] = Field3D(temperature);
    species["density"] = Field3D(density);
    species["pressure"] = Field3D(density * temperature);
    species["collision_frequencies"]["test+_test+_coll"] = Field3D(nu);
    return state;
  }

  void run(Options& state, Options& diagnostics) {
    BraginskiiConduction component("braginskii_conduction", options, nullptr);
    component.declareAllSpecies({"test+"});
    component.transform(state);
    component.outputVars(diagnostics);
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

TEST_F(BraginskiiConductionTest, TransitDisabledByDefault) {
  auto state = uniformState();
  run(state, output_base);
  EXPECT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["test+"]["kappa_par"]), 3.9,
                           "RGN_NOBNDRY"));
  EXPECT_FALSE(output_base.isSet("kappa_transit_factor_test+"));
}

TEST_F(BraginskiiConductionTest, TransitNoGridDefaultsAndZeroGradient) {
  options["test+"]["kappa_transit_limit"] = true;
  auto state = uniformState();
  run(state, output_base);
  EXPECT_TRUE(IsFieldEqual(get<Field3D>(output_base["kappa_transit_length_test+"]), 8.0));
  EXPECT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["test+"]["kappa_par"]),
                           3.9 * 8.0 / (3.9 + 8.0), "RGN_NOBNDRY"));
  EXPECT_TRUE(IsFieldEqual(getDeriv(state), 0.0, "RGN_NOBNDRY"));
  EXPECT_TRUE(IsFieldEqual(get<Field3D>(output_base["kappa_transit_factor_test+"]),
                           8.0 / (3.9 + 8.0), "RGN_NOBNDRY"));
}

TEST_F(BraginskiiConductionTest, TransitMeshGeometryAndNormalisation) {
  options["test+"]["kappa_transit_limit"] = true;
  options["units"]["meters"].force(0.5);
  Field3D q = linearGradient(-2., 0., 0., -1., 0., 0., 0.);
  static_cast<FakeMesh*>(mesh)->setGridDataSource(
      new FakeGridDataSource({{"q", q}, {"Rxy", 1.5}}));
  auto state = uniformState();
  run(state, output_base);
  const Field3D length = get<Field3D>(output_base["kappa_transit_length_test+"]);
  const Field3D kappa = get<Field3D>(state["species"]["test+"]["kappa_par"]);
  BOUT_FOR_SERIAL(i, length.getRegion("RGN_NOBNDRY")) {
    EXPECT_DOUBLE_EQ(length[i], std::abs(q[i]) * 3.0);
    EXPECT_NEAR(kappa[i], 3.9 * length[i] / (3.9 + length[i]), 1e-14);
  }
}

TEST_F(BraginskiiConductionTest, TransitCustomMeshNamesAndPartialFallback) {
  options["test+"]["kappa_transit_limit"] = true;
  options["test+"]["kappa_transit_q_field"] = "q95";
  static_cast<FakeMesh*>(mesh)->setGridDataSource(new FakeGridDataSource({{"q95", 3.0}}));
  auto state = uniformState();
  run(state, output_base);
  EXPECT_TRUE(IsFieldEqual(get<Field3D>(output_base["kappa_transit_length_test+"]), 6.0));
}

TEST_F(BraginskiiConductionTest, TransitExplicitGeometryOverridesMesh) {
  options["test+"]["kappa_transit_limit"] = true;
  options["test+"]["kappa_transit_q"] = -3.0;
  options["test+"]["kappa_transit_R"] = 2.0;
  static_cast<FakeMesh*>(mesh)->setGridDataSource(
      new FakeGridDataSource({{"q", 9.0}, {"Rxy", 4.0}}));
  auto state = uniformState();
  run(state, output_base);
  EXPECT_TRUE(IsFieldEqual(get<Field3D>(output_base["kappa_transit_length_test+"]), 6.0));
}

TEST_F(BraginskiiConductionTest, TransitExplicitLengthOverridesGeometry) {
  options["test+"]["kappa_transit_limit"] = true;
  options["test+"]["kappa_transit_length"] = 3.0;
  options["test+"]["kappa_transit_q"] = 0.0; // Unused when length is supplied
  options["test+"]["kappa_transit_alpha"] = 0.5;
  auto state = uniformState();
  run(state, output_base);
  EXPECT_TRUE(IsFieldEqual(get<Field3D>(output_base["kappa_transit_length_test+"]), 3.0));
  EXPECT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["test+"]["kappa_par"]),
                           3.9 * 1.5 / (3.9 + 1.5), "RGN_NOBNDRY"));
}

TEST_F(BraginskiiConductionTest, TransitCollisionalAndCollisionlessLimits) {
  options["test+"]["kappa_transit_limit"] = true;
  for (const BoutReal nu : {1e8, 1e-8}) {
    auto state = uniformState(nu);
    auto diagnostics = output_base.copy();
    run(state, diagnostics);
    const Field3D kappa = get<Field3D>(state["species"]["test+"]["kappa_par"]);
    const BoutReal limit = nu > 1.0 ? 3.9 / nu : 8.0;
    BOUT_FOR_SERIAL(i, kappa.getRegion("RGN_NOBNDRY")) {
      EXPECT_NEAR(kappa[i] / limit, 1.0, 1e-7);
      EXPECT_LE(kappa[i], 8.0);
    }
  }
}

TEST_F(BraginskiiConductionTest, TransitVacuumAndColdCellsAreFinite) {
  options["test+"]["kappa_transit_limit"] = true;
  for (const auto& values : {std::pair{0.0, 1.0}, std::pair{2.0, 0.0}}) {
    auto state = uniformState(1.0, values.first, values.second);
    auto diagnostics = output_base.copy();
    run(state, diagnostics);
    const Field3D kappa = get<Field3D>(state["species"]["test+"]["kappa_par"]);
    BOUT_FOR_SERIAL(i, kappa.getRegion("RGN_NOBNDRY")) {
      EXPECT_TRUE(std::isfinite(kappa[i]));
      EXPECT_DOUBLE_EQ(kappa[i], 0.0);
    }
  }
}

TEST_F(BraginskiiConductionTest, TransitRejectsInvalidParameters) {
  options["test+"]["kappa_transit_limit"] = true;
  for (const auto& key : {"kappa_transit_alpha", "kappa_transit_length"}) {
    for (const BoutReal value : {0.0, -1.0, std::numeric_limits<BoutReal>::infinity()}) {
      auto invalid = options.copy();
      invalid["test+"][key] = value;
      EXPECT_THROW(BraginskiiConduction("braginskii_conduction", invalid, nullptr),
                   BoutException);
    }
  }
}

TEST_F(BraginskiiConductionTest, TransitLimitsSmallPerturbationResponse) {
  // The existing flux limiter alone leaves the Braginskii linear response.
  options["test+"]["kappa_limit_alpha"] = 0.2;
  const BoutReal factor = 8.0 / (3.9 + 8.0);
  for (const BoutReal amplitude : {1e-5, 1e-7}) {
    auto state = uniformState();
    const Field3D temperature = 2.0 + amplitude * temp1;
    state["species"]["test+"]["temperature"] = temperature;
    state["species"]["test+"]["pressure"] = temperature;
    auto limited = state.copy();
    auto diagnostics = output_base.copy();
    options["test+"]["kappa_transit_limit"].force(false);
    run(state, diagnostics);
    diagnostics = output_base.copy();
    options["test+"]["kappa_transit_limit"].force(true);
    run(limited, diagnostics);
    const Field3D baseline = getDeriv(state);
    const Field3D response = getDeriv(limited);
    BOUT_FOR_SERIAL(i, response.getRegion("RGN_NOBNDRY")) {
      ASSERT_NE(baseline[i], 0.0);
      EXPECT_NEAR(response[i] / baseline[i], factor, 100.0 * amplitude);
    }
  }
}

TEST_F(BraginskiiConductionTest, TransitMassAndDensityScaling) {
  options["test+"]["kappa_transit_limit"] = true;
  for (const BoutReal mass : {1.0 / 1836.0, 2.0, 12.0}) {
    auto state = uniformState(1e-8, 2.0, 3.0);
    state["species"]["test+"]["AA"].force(mass);
    auto diagnostics = output_base.copy();
    run(state, diagnostics);
    const Field3D kappa = get<Field3D>(state["species"]["test+"]["kappa_par"]);
    const BoutReal cap = 3.0 * std::sqrt(2.0 / mass) * 8.0;
    BOUT_FOR_SERIAL(i, kappa.getRegion("RGN_NOBNDRY")) {
      EXPECT_NEAR(kappa[i] / cap, 1.0, 1e-7);
      EXPECT_LE(kappa[i], cap);
    }
  }
}

// TODO: Add tests for more types of collisions (e.g., electron-ion, charge-exchante, ion
// neutral...)
