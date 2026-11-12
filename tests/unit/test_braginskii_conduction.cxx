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

TEST_F(BraginskiiConductionTest, TransitGeometryAndDefaults) {
  auto state = uniformState();
  run(state, output_base);
  EXPECT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["test+"]["kappa_par"]), 3.9,
                           "RGN_NOBNDRY"));
  EXPECT_FALSE(output_base.isSet("kappa_transit_factor_test+"));

  options["test+"]["kappa_transit_limit"].force(true);
  auto check = [&](const char* name, const Field3D& length, BoutReal alpha = 1.) {
    SCOPED_TRACE(name);
    options["test+"]["kappa_transit_alpha"].force(alpha);
    auto state = uniformState();
    auto diagnostics = output_base.copy();
    run(state, diagnostics);
    const Field3D factor = alpha * length / (3.9 + alpha * length);
    const Field3D kappa = 3.9 * factor;
    EXPECT_TRUE(
        IsFieldEqual(get<Field3D>(diagnostics["kappa_transit_length_test+"]), length));
    EXPECT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["test+"]["kappa_par"]), kappa,
                             "RGN_NOBNDRY", 1e-14));
    EXPECT_TRUE(IsFieldEqual(get<Field3D>(diagnostics["kappa_transit_factor_test+"]),
                             factor, "RGN_NOBNDRY"));
    EXPECT_TRUE(IsFieldEqual(getDeriv(state), 0., "RGN_NOBNDRY"));
  };
  check("no grid", Field3D(8.));

  const Field3D q = linearGradient(-2., 0., 0., -1., 0., 0., 0.);
  auto* grid = static_cast<FakeMesh*>(mesh);
  grid->setGridDataSource(new FakeGridDataSource({{"q", q}, {"Rxy", 1.5}}));
  options["units"]["meters"].force(0.5);
  check("mesh and normalisation", abs(q) * 3.);

  options["units"]["meters"].force(1.);
  options["test+"]["kappa_transit_q_field"].force("q95");
  grid->setGridDataSource(new FakeGridDataSource({{"q95", 3.}}));
  check("custom name and partial fallback", Field3D(6.));

  options["test+"]["kappa_transit_q_field"].force("q");
  options["test+"]["kappa_transit_q"] = -3.;
  options["test+"]["kappa_transit_R"] = 2.;
  grid->setGridDataSource(new FakeGridDataSource({{"q", 9.}, {"Rxy", 4.}}));
  check("geometry overrides", Field3D(6.));

  options["test+"]["kappa_transit_length"] = 3.;
  options["test+"]["kappa_transit_q"].force(0.);
  check("length overrides invalid q", Field3D(3.), 0.5);
}

TEST_F(BraginskiiConductionTest, TransitLimitsAndScaling) {
  options["test+"]["kappa_transit_limit"] = true;
  const struct {
    const char* name;
    BoutReal nu, temperature, density, mass, expected;
  } cases[] = {
      {"collisional", 1e8, 2., 1., 2., 3.9e-8},
      {"collisionless", 1e-8, 2., 1., 2., 8.},
      {"cold", 1., 0., 1., 2., 0.},
      {"vacuum", 1., 2., 0., 2., 0.},
      {"electron mass", 1e-8, 2., 3., 1. / 1836., 24. * std::sqrt(2. * 1836.)},
      {"deuterium mass", 1e-8, 2., 3., 2., 24.},
      {"carbon mass", 1e-8, 2., 3., 12., 24. * std::sqrt(2. / 12.)},
  };
  for (const auto& test : cases) {
    SCOPED_TRACE(test.name);
    auto state = uniformState(test.nu, test.temperature, test.density);
    state["species"]["test+"]["AA"].force(test.mass);
    auto diagnostics = output_base.copy();
    run(state, diagnostics);
    const Field3D kappa = get<Field3D>(state["species"]["test+"]["kappa_par"]);
    const BoutReal cap = test.density * std::sqrt(test.temperature / test.mass) * 8.;
    BOUT_FOR_SERIAL(i, kappa.getRegion("RGN_NOBNDRY")) {
      EXPECT_TRUE(std::isfinite(kappa[i]));
      EXPECT_NEAR(kappa[i], test.expected, 1e-7 * test.expected);
      EXPECT_LE(kappa[i], cap);
    }
  }
}

TEST_F(BraginskiiConductionTest, TransitRejectsInvalidParameters) {
  // Expected exceptions do not need expensive backtraces. Restore on scope exit.
  struct BacktraceGuard {
    BacktraceGuard() { BoutException::disableBacktrace(); }
    ~BacktraceGuard() { BoutException::enableBacktrace(); }
  } guard;
  options["test+"]["kappa_transit_limit"] = true;
  for (const auto& key : {"kappa_transit_alpha", "kappa_transit_length"}) {
    SCOPED_TRACE(key);
    for (const BoutReal value : {0.0, -1.0, std::numeric_limits<BoutReal>::infinity()}) {
      SCOPED_TRACE(value);
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

// TODO: Add tests for more types of collisions (e.g., electron-ion, charge-exchante, ion
// neutral...)
