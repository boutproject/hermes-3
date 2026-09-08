
#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/neutral_mixed.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using NeutralMixedTest = FakeMeshFixture;
// Component test checking only that the state
// includes Nd, Pd, NVd as evolved variables.
TEST_F(NeutralMixedTest, CreateComponent) {
  FakeSolver solver;

  Options options{
      {"units",
       {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}, {"meters", 1.0}}},
      {"d", {{"type", "neutral_mixed"}, {"AA", 2.0}, {"evolve_momentum", true}}}};
  NeutralMixed component("d", options, &solver);

  Options state = solver.getState();

  EXPECT_TRUE(state.isSet("Nd"));
  EXPECT_TRUE(state.isSet("Pd"));
  EXPECT_TRUE(state.isSet("NVd"));
}
// Component test checking only that the state
// includes Nd, Pd as evolved variables when evolve_momentum = false.
TEST_F(NeutralMixedTest, CreateComponentEvolveMomentumFalse) {
  FakeSolver solver;

  Options options{
      {"units",
       {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}, {"meters", 1.0}}},
      {"d", {{"type", "neutral_mixed"}, {"AA", 2.0}, {"evolve_momentum", false}}}};
  NeutralMixed component("d", options, &solver);

  Options state = solver.getState();

  EXPECT_TRUE(state.isSet("Nd"));
  EXPECT_TRUE(state.isSet("Pd"));
  EXPECT_FALSE(state.isSet("NVd"));
}
// Transform test checking only that the state has data set in
// the the required auxiliary variables. Note that boundary conditions
// are also set in the evolved variables.
TEST_F(NeutralMixedTest, Transform) {
  FakeSolver solver;

  Options options{
      {"units",
       {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}, {"meters", 1.0}}},
      {"d", {{"type", "neutral_mixed"}, {"AA", 2.0}, {"evolve_momentum", true}}}};
  NeutralMixed component("d", options, &solver);

  Options state;
  component.transform(state);

  Options& species = state["species"]["d"];
  EXPECT_TRUE(species.isSet("density"));
  EXPECT_TRUE(species.isSet("AA"));
  EXPECT_TRUE(species.isSet("pressure"));
  EXPECT_TRUE(species.isSet("momentum"));
  EXPECT_TRUE(species.isSet("velocity"));
  EXPECT_TRUE(species.isSet("temperature"));
}
// Test of finally() for this component following
// tests of evolve_density and evolve_pressure.
// We provide a state vector where all fields are constants.
// We use finally to compute the ddt() of the evolved
// state variables. For these inputs, the only non-zero
// contribution comes from the sources, which may be tested
// by checking that the ddt() variables match expected constants.
// Note that this test does not check the definition of the
// differential operators in the right-hand-side of the equations.
TEST_F(NeutralMixedTest, Finally) {
  FakeSolver solver;

  Options options{
      {"units",
       {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}, {"meters", 1.0}}},
      {"d", {{"type", "neutral_mixed"}, {"AA", 2.0}, {"evolve_momentum", true}}}};
  NeutralMixed component("d", options, &solver);

  // Call the finally() method with a density, energy, and momentum source
  const Options state = {{"species",
                          {{"d",
                            {{"density", 1.0},
                             {"density_source", 0.5},
                             {"pressure", 1.0},
                             {"energy_source", 1.5},
                             {"momentum", 1.0},
                             {"momentum_source", 0.75},
                             {"temperature", 1.0},
                             {"velocity", 1.0}}}}}};
  component.finally(state);

  Options ddt = solver.getTimeDerivs();

  EXPECT_TRUE(ddt.isSet("Nd"));
  Field3D ddt_Nd = ddt["Nd"].as<Field3D>();

  BOUT_FOR_SERIAL(i, ddt_Nd.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(0.5, ddt_Nd[i]);
  }

  EXPECT_TRUE(ddt.isSet("Pd"));
  Field3D ddt_Pd = ddt["Pd"].as<Field3D>();

  BOUT_FOR_SERIAL(i, ddt_Pd.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(1.0, ddt_Pd[i]);
  }

  EXPECT_TRUE(ddt.isSet("NVd"));
  Field3D ddt_NVd = ddt["NVd"].as<Field3D>();

  BOUT_FOR_SERIAL(i, ddt_NVd.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(0.75, ddt_NVd[i]);
  }
}
// Identical to the test above, but using the collisionality_override variable.
TEST_F(NeutralMixedTest, FinallyCollisionalityOverride) {
  FakeSolver solver;

  Options options{
      {"units",
       {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}, {"meters", 1.0}}},
      {"d",
       {{"type", "neutral_mixed"},
        {"AA", 2.0},
        {"evolve_momentum", true},
        {"collisionality_override", 1.0}}}};
  NeutralMixed component("d", options, &solver);

  // Call the finally() method with a density, energy, and momentum source
  const Options state = {{"species",
                          {{"d",
                            {{"density", 1.0},
                             {"density_source", 0.5},
                             {"pressure", 1.0},
                             {"energy_source", 1.5},
                             {"momentum", 1.0},
                             {"momentum_source", 0.75},
                             {"temperature", 1.0},
                             {"velocity", 1.0}}}}}};
  component.finally(state);

  Options ddt = solver.getTimeDerivs();

  EXPECT_TRUE(ddt.isSet("Nd"));
  Field3D ddt_Nd = ddt["Nd"].as<Field3D>();

  BOUT_FOR_SERIAL(i, ddt_Nd.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(0.5, ddt_Nd[i]);
  }

  EXPECT_TRUE(ddt.isSet("Pd"));
  Field3D ddt_Pd = ddt["Pd"].as<Field3D>();

  BOUT_FOR_SERIAL(i, ddt_Pd.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(1.0, ddt_Pd[i]);
  }

  EXPECT_TRUE(ddt.isSet("NVd"));
  Field3D ddt_NVd = ddt["NVd"].as<Field3D>();

  BOUT_FOR_SERIAL(i, ddt_NVd.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(0.75, ddt_NVd[i]);
  }
}
// Identical to the test above, but using evolve_momentum = false.
TEST_F(NeutralMixedTest, FinallyEvolveMomentumFalse) {
  FakeSolver solver;

  Options options{
      {"units",
       {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}, {"meters", 1.0}}},
      {"d", {{"type", "neutral_mixed"}, {"AA", 2.0}, {"evolve_momentum", false}}}};
  NeutralMixed component("d", options, &solver);

  // Call the finally() method with a density, energy, and momentum source
  const Options state = {{"species",
                          {{"d",
                            {{"density", 1.0},
                             {"density_source", 0.5},
                             {"pressure", 1.0},
                             {"energy_source", 1.5},
                             {"momentum", 1.0},
                             {"momentum_source", 0.75},
                             {"temperature", 1.0},
                             {"velocity", 1.0}}}}}};
  component.finally(state);

  Options ddt = solver.getTimeDerivs();

  EXPECT_TRUE(ddt.isSet("Nd"));
  Field3D ddt_Nd = ddt["Nd"].as<Field3D>();

  BOUT_FOR_SERIAL(i, ddt_Nd.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(0.5, ddt_Nd[i]);
  }

  EXPECT_TRUE(ddt.isSet("Pd"));
  Field3D ddt_Pd = ddt["Pd"].as<Field3D>();

  BOUT_FOR_SERIAL(i, ddt_Pd.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(1.0, ddt_Pd[i]);
  }

  EXPECT_FALSE(ddt.isSet("NVd"));
}
// Identical to the test above, but using nonorthogonal_operators = true.
TEST_F(NeutralMixedTest, FinallyNonorthogonalOperators) {
  FakeSolver solver;

  Options options{
      {"units",
       {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}, {"meters", 1.0}}},
      {"d",
       {{"type", "neutral_mixed"},
        {"AA", 2.0},
        {"evolve_momentum", true},
        {"nonorthogonal_operators", true}}}};
  NeutralMixed component("d", options, &solver);

  // Call the finally() method with a density, energy, and momentum source
  const Options state = {{"species",
                          {{"d",
                            {{"density", 1.0},
                             {"density_source", 0.5},
                             {"pressure", 1.0},
                             {"energy_source", 1.5},
                             {"momentum", 1.0},
                             {"momentum_source", 0.75},
                             {"temperature", 1.0},
                             {"velocity", 1.0}}}}}};
  component.finally(state);

  Options ddt = solver.getTimeDerivs();

  EXPECT_TRUE(ddt.isSet("Nd"));
  Field3D ddt_Nd = ddt["Nd"].as<Field3D>();

  BOUT_FOR_SERIAL(i, ddt_Nd.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(0.5, ddt_Nd[i]);
  }

  EXPECT_TRUE(ddt.isSet("Pd"));
  Field3D ddt_Pd = ddt["Pd"].as<Field3D>();

  BOUT_FOR_SERIAL(i, ddt_Pd.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(1.0, ddt_Pd[i]);
  }

  EXPECT_TRUE(ddt.isSet("NVd"));
  Field3D ddt_NVd = ddt["NVd"].as<Field3D>();

  BOUT_FOR_SERIAL(i, ddt_NVd.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(0.75, ddt_NVd[i]);
  }
}

// General purpose test helper. Has options to enable or disable collision frequencies
// as needed by the tests. Features radial and poloidal variations in P, T and V.
namespace {

Options runNeutralMixedTest(Options options, bool with_collisions = false) {
  FakeSolver solver;

  options["units"] = {
      {"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}, {"meters", 1.0}};

  // Effectively disable the gradient floors and ceilings - this is necessary
  // because those are tuned to realistic conditions and the tests only check
  // limiter behaviour at the moment. A test that needs a particular value sets
  // it itself, and what it sets is left alone here.
  const auto loosen = [&options](const std::string& key, BoutReal value) {
    if (!options["d"].isSet(key)) {
      options["d"][key] = value;
    }
  };
  loosen("limiter_gradient_floor", 1e-10);
  loosen("limiter_gradient_ceiling", 1e10);
  loosen("limiter_gradient_floor_eta", 1e-10);
  loosen("limiter_gradient_ceiling_eta", 1e10);

  NeutralMixed component("d", options, &solver);

  Field3D Pn =
      makeField<Field3D>([](Ind3D& i) { return 1.0 + i.x() + 0.5 * i.y(); }, mesh);
  Field3D Tn =
      makeField<Field3D>([](Ind3D& i) { return 1.0 + 0.25 * i.x() + 0.5 * i.y(); }, mesh);
  Field3D Vn =
      makeField<Field3D>([](Ind3D& i) { return 1.0 + 0.5 * i.x() + 0.25 * i.y(); }, mesh);

  // Call the finally() method with a density, energy, and momentum source
  Options state = {{"species",
                    {{"d",
                      {{"density", 1.0},
                       {"density_source", 0.5},
                       {"pressure", Pn},
                       {"energy_source", 1.5},
                       {"momentum", 1.0},
                       {"momentum_source", 0.75},
                       {"temperature", Tn},
                       {"velocity", Vn}}}}}};

  // Simulate collision frequencies if needed
  if (with_collisions) {
    state["species"]["d"]["collision_frequency"] = 1.0;
    state["species"]["d"]["collision_frequencies"]["d_d+_cx"] = 100.0;
  }

  component.finally(state);

  // Construct state with norms for outputVars to add diagnostics to
  Options out{
      {"Nnorm", 1.0}, {"Tnorm", 1.0}, {"Omega_ci", 1.0}, {"Cs0", 1.0}, {"rho_s0", 1.0}};

  component.outputVars(out);

  return out;
}
} // namespace

// Check that flux limiter reduces to a simple harmonic mean when sharpness = 1.0.
TEST_F(NeutralMixedTest, DnnHarmonicLimiter) {

  Options out = runNeutralMixedTest({{"d",
                                      {
                                          {"type", "neutral_mixed"},
                                          {"diagnose", true},
                                          {"AA", 2.0},
                                          {"flux_limiter_sharpness", 1.0},
                                      }}});

  Field3D Dnn = out["Dnnd"].as<Field3D>();
  Field3D Dunl = out["Dnnd_unlimited"].as<Field3D>();
  Field3D Dmax = out["Dnnd_max"].as<Field3D>();

  BOUT_FOR_SERIAL(i, Dnn.getRegion("RGN_NOBNDRY")) {
    EXPECT_DOUBLE_EQ(Dnn[i], Dunl[i] * Dmax[i] / (Dunl[i] + Dmax[i]));
  }
}

// Check that an aggressive flux limit limits the flux.
TEST_F(NeutralMixedTest, DnnTightLimit) {

  Options out = runNeutralMixedTest({{"d",
                                      {{"type", "neutral_mixed"},
                                       {"diagnose", true},
                                       {"AA", 2.0},
                                       {"flux_limit", 1e-5}}}});

  Field3D Dnn = out["Dnnd"].as<Field3D>();
  Field3D Dunl = out["Dnnd_unlimited"].as<Field3D>();
  Field3D Dmax = out["Dnnd_max"].as<Field3D>();

  BOUT_FOR_SERIAL(i, Dnn.getRegion("RGN_NOBNDRY")) {
    EXPECT_LT(Dmax[i], Dunl[i]);
    EXPECT_LT(Dnn[i], Dunl[i]);
    EXPECT_NEAR(Dnn[i], Dmax[i], Dnn[i] * 1e-3);
  }
}

// Check that a loose flux limit doesn't limit the flux.
TEST_F(NeutralMixedTest, DnnLooseLimit) {

  Options out = runNeutralMixedTest({{"d",
                                      {{"type", "neutral_mixed"},
                                       {"diagnose", true},
                                       {"AA", 2.0},
                                       {"flux_limit", 1e6}}}});

  Field3D Dnn = out["Dnnd"].as<Field3D>();
  Field3D Dunl = out["Dnnd_unlimited"].as<Field3D>();
  Field3D Dmax = out["Dnnd_max"].as<Field3D>();

  BOUT_FOR_SERIAL(i, Dnn.getRegion("RGN_NOBNDRY")) {
    EXPECT_GT(Dmax[i], Dunl[i]);
    EXPECT_NEAR(Dnn[i], Dunl[i], Dnn[i] * 1e-3);
  }
}

// Check that a tight conduction flux limit limits the conductive flux.
TEST_F(NeutralMixedTest, ConductionTightLimit) {

  Options out = runNeutralMixedTest({{"d",
                                      {{"type", "neutral_mixed"},
                                       {"diagnose", true},
                                       {"AA", 2.0},
                                       {"combined_limiters", false},
                                       {"flux_limit_cond_perp", 1e-5},
                                       {"flux_limit_cond_par", 1e-5}}}});

  Field3D kappa_perp = out["kappa_d_perp"].as<Field3D>();
  Field3D kappa_par = out["kappa_d_par"].as<Field3D>();
  Field3D kappa_unlimited = out["kappa_d_unlimited"].as<Field3D>();
  Field3D kappa_max_perp = out["kappa_d_max_perp"].as<Field3D>();
  Field3D kappa_max_par = out["kappa_d_max_par"].as<Field3D>();

  BOUT_FOR_SERIAL(i, kappa_unlimited.getRegion("RGN_NOBNDRY")) {
    EXPECT_LT(kappa_max_perp[i], kappa_unlimited[i]);
    EXPECT_LT(kappa_perp[i], kappa_unlimited[i]);
    EXPECT_NEAR(kappa_perp[i], kappa_max_perp[i], kappa_perp[i] * 1e-3);
    EXPECT_NEAR(kappa_par[i], kappa_max_par[i], kappa_par[i] * 1e-3);
  }
}

// Check that a loose conduction flux limit does not limit the conductive flux.
TEST_F(NeutralMixedTest, ConductionLooseLimit) {

  Options out = runNeutralMixedTest({{"d",
                                      {{"type", "neutral_mixed"},
                                       {"diagnose", true},
                                       {"AA", 2.0},
                                       {"combined_limiters", false},
                                       {"flux_limit_cond_perp", 1e6},
                                       {"flux_limit_cond_par", 1e6}}}});

  Field3D kappa_perp = out["kappa_d_perp"].as<Field3D>();
  Field3D kappa_par = out["kappa_d_par"].as<Field3D>();
  Field3D kappa_unlimited = out["kappa_d_unlimited"].as<Field3D>();
  Field3D kappa_max_perp = out["kappa_d_max_perp"].as<Field3D>();

  BOUT_FOR_SERIAL(i, kappa_unlimited.getRegion("RGN_NOBNDRY")) {
    EXPECT_GT(kappa_max_perp[i], kappa_unlimited[i]);
    EXPECT_NEAR(kappa_perp[i], kappa_unlimited[i], kappa_perp[i] * 1e-3);
    EXPECT_NEAR(kappa_par[i], kappa_unlimited[i], kappa_par[i] * 1e-3);
  }
}

// Check that a tight viscosity flux limit limits the viscous flux.
TEST_F(NeutralMixedTest, ViscosityTightLimit) {

  Options out = runNeutralMixedTest({{"d",
                                      {{"type", "neutral_mixed"},
                                       {"diagnose", true},
                                       {"AA", 2.0},
                                       {"combined_limiters", false},
                                       {"flux_limit_visc_perp", 1e-5},
                                       {"flux_limit_visc_par", 1e-5}}}});

  Field3D eta_perp = out["eta_d_perp"].as<Field3D>();
  Field3D eta_par = out["eta_d_par"].as<Field3D>();
  Field3D eta_unlimited = out["eta_d_unlimited"].as<Field3D>();
  Field3D eta_max_perp = out["eta_d_max_perp"].as<Field3D>();
  Field3D eta_max_par = out["eta_d_max_par"].as<Field3D>();

  BOUT_FOR_SERIAL(i, eta_unlimited.getRegion("RGN_NOBNDRY")) {
    EXPECT_LT(eta_max_perp[i], eta_unlimited[i]);
    EXPECT_LT(eta_perp[i], eta_unlimited[i]);
    EXPECT_NEAR(eta_perp[i], eta_max_perp[i], eta_perp[i] * 1e-3);
    EXPECT_NEAR(eta_par[i], eta_max_par[i], eta_par[i] * 1e-3);
  }
}

// Check that a loose viscosity flux limit does not limit the viscous flux.
TEST_F(NeutralMixedTest, ViscosityLooseLimit) {

  Options out = runNeutralMixedTest({{"d",
                                      {{"type", "neutral_mixed"},
                                       {"diagnose", true},
                                       {"AA", 2.0},
                                       {"combined_limiters", false},
                                       {"flux_limit_visc_perp", 1e6},
                                       {"flux_limit_visc_par", 1e6}}}});

  Field3D eta_perp = out["eta_d_perp"].as<Field3D>();
  Field3D eta_par = out["eta_d_par"].as<Field3D>();
  Field3D eta_unlimited = out["eta_d_unlimited"].as<Field3D>();
  Field3D eta_max_perp = out["eta_d_max_perp"].as<Field3D>();

  BOUT_FOR_SERIAL(i, eta_unlimited.getRegion("RGN_NOBNDRY")) {
    EXPECT_GT(eta_max_perp[i], eta_unlimited[i]);
    EXPECT_NEAR(eta_perp[i], eta_unlimited[i], eta_perp[i] * 1e-3);
    EXPECT_NEAR(eta_par[i], eta_unlimited[i], eta_par[i] * 1e-3);
  }
}

// Check that the explicit diffusion limit can override the flux limitation.
// Dmax is a harmonic mean of the flux limit and the explicit limit, so it will never
// be exactly equal to the explicit limit.
// This test checks if both Dmax and Dnn are within 1% of an explicit limit, and whether
// Dmax is always less than the explicit limit (as expected for a harmonic mean).
TEST_F(NeutralMixedTest, DnnExplicitLimit) {

  Options out = runNeutralMixedTest({{"d",
                                      {{"type", "neutral_mixed"},
                                       {"diagnose", true},
                                       {"AA", 2.0},
                                       {"flux_limit", 0.2},
                                       {"diffusion_limit", 1e-5}}}});

  Field3D Dnn = out["Dnnd"].as<Field3D>();
  Field3D Dmax = out["Dnnd_max"].as<Field3D>();

  BOUT_FOR_SERIAL(i, Dnn.getRegion("RGN_NOBNDRY")) {
    EXPECT_NEAR(Dmax[i], 1e-5, 1e-2 * Dmax[i]);
    EXPECT_NEAR(Dnn[i], 1e-5, 1e-2 * Dnn[i]);
    EXPECT_LT(Dmax[i], 1e-5);
  }
}

// Check that adding collisionality reduces Dnn.
TEST_F(NeutralMixedTest, DnnCollisionalityImpact) {

  Options options = {{"d", {{"type", "neutral_mixed"}, {"diagnose", true}, {"AA", 2.0}}}};

  Field3D Dnn = runNeutralMixedTest(options.copy(), false)["Dnnd"].as<Field3D>();
  Field3D Dnn_coll = runNeutralMixedTest(options.copy(), true)["Dnnd"].as<Field3D>();

  BOUT_FOR_SERIAL(i, Dnn.getRegion("RGN_NOBNDRY")) { EXPECT_LT(Dnn_coll[i], Dnn[i]); }
}

// Check that reducing neutral_lmax raises collisionality floor.
// Lower neutral_lmax leads to a higher nu_pseudo_mfp,
// which increases total nu. and reduces Dnn.
TEST_F(NeutralMixedTest, DnnCollisionalityFloor) {

  Field3D Dnn_lo_lmax = runNeutralMixedTest({{"d",
                                              {{"type", "neutral_mixed"},
                                               {"diagnose", true},
                                               {"AA", 2.0},
                                               {"neutral_lmax", 0.01}}}})["Dnnd"]
                            .as<Field3D>();

  Field3D Dnn_hi_lmax = runNeutralMixedTest({{"d",
                                              {{"type", "neutral_mixed"},
                                               {"diagnose", true},
                                               {"AA", 2.0},
                                               {"neutral_lmax", 100}}}})["Dnnd"]
                            .as<Field3D>();

  BOUT_FOR_SERIAL(i, Dnn_lo_lmax.getRegion("RGN_NOBNDRY")) {
    EXPECT_LT(Dnn_lo_lmax[i], Dnn_hi_lmax[i]);
  }
}

// With combined limiters, conduction and viscosity are derived from limited
// Dnn instead of being separately calculated. There is no separate limitation
// for perpendicular and parallel limiters.
TEST_F(NeutralMixedTest, CombinedLimitersDeriveFromDnn) {
  const BoutReal AA = 2.0;
  const BoutReal Nn = 1.0; // density set in the state below

  Options out = runNeutralMixedTest({{"d",
                                      {{"type", "neutral_mixed"},
                                       {"diagnose", true},
                                       {"AA", AA},
                                       {"combined_limiters", true}}}});

  Field3D Dnn = out["Dnnd"].as<Field3D>();
  Field3D kappa_perp = out["kappa_d_perp"].as<Field3D>();
  Field3D kappa_par = out["kappa_d_par"].as<Field3D>();
  Field3D eta_perp = out["eta_d_perp"].as<Field3D>();
  Field3D eta_par = out["eta_d_par"].as<Field3D>();

  BOUT_FOR_SERIAL(i, Dnn.getRegion("RGN_NOBNDRY")) {
    EXPECT_DOUBLE_EQ(kappa_perp[i], (5. / 2) * Nn * Dnn[i]);
    EXPECT_DOUBLE_EQ(kappa_par[i], kappa_perp[i]);
    EXPECT_DOUBLE_EQ(eta_perp[i], (2. / 5) * AA * kappa_perp[i]);
    EXPECT_DOUBLE_EQ(eta_par[i], eta_perp[i]);
  }
}

// Disabling the viscosity limiter leaves eta unlimited.
TEST_F(NeutralMixedTest, ViscosityLimiterOffLeavesEtaUnlimited) {

  Options out = runNeutralMixedTest({{"d",
                                      {{"type", "neutral_mixed"},
                                       {"diagnose", true},
                                       {"AA", 2.0},
                                       {"combined_limiters", false},
                                       {"flux_limit_visc_perp", -1.0},
                                       {"flux_limit_visc_par", -1.0}}}});

  Field3D eta_perp = out["eta_d_perp"].as<Field3D>();
  Field3D eta_par = out["eta_d_par"].as<Field3D>();
  Field3D eta_unlimited = out["eta_d_unlimited"].as<Field3D>();

  BOUT_FOR_SERIAL(i, eta_unlimited.getRegion("RGN_NOBNDRY")) {
    EXPECT_GT(eta_unlimited[i], 0.0);
    EXPECT_DOUBLE_EQ(eta_perp[i], eta_unlimited[i]);
    EXPECT_DOUBLE_EQ(eta_par[i], eta_unlimited[i]);
  }
}

// Disabling the conduction limiter leaves kappa unlimited.
TEST_F(NeutralMixedTest, ConductionLimiterOffLeavesKappaUnlimited) {

  Options out = runNeutralMixedTest({{"d",
                                      {{"type", "neutral_mixed"},
                                       {"diagnose", true},
                                       {"AA", 2.0},
                                       {"combined_limiters", false},
                                       {"flux_limit_cond_perp", -1.0},
                                       {"flux_limit_cond_par", -1.0}}}});

  Field3D kappa_perp = out["kappa_d_perp"].as<Field3D>();
  Field3D kappa_par = out["kappa_d_par"].as<Field3D>();
  Field3D kappa_unlimited = out["kappa_d_unlimited"].as<Field3D>();

  BOUT_FOR_SERIAL(i, kappa_unlimited.getRegion("RGN_NOBNDRY")) {
    EXPECT_GT(kappa_unlimited[i], 0.0);
    EXPECT_DOUBLE_EQ(kappa_perp[i], kappa_unlimited[i]);
    EXPECT_DOUBLE_EQ(kappa_par[i], kappa_unlimited[i]);
  }
}

// Setting a per-channel option while the limiters are combined throws exception.
TEST_F(NeutralMixedTest, CombinedLimitersRejectPerChannelOptions) {
  EXPECT_THROW(runNeutralMixedTest({{"d",
                                     {{"type", "neutral_mixed"},
                                      {"diagnose", true},
                                      {"AA", 2.0},
                                      {"combined_limiters", true},
                                      {"flux_limit_cond_perp", 0.5}}}}),
               BoutException);
}

/////////////////////////////////////////////////////////////////////////////////
// SHEATH TESTS
/////////////////////////////////////////////////////////////////////////////////

namespace {

// Test for zero at boundary
void expectZeroedAtTarget(const Field3D& f) {
  Mesh* localmesh = f.getMesh();
  for (int x = localmesh->xstart; x <= localmesh->xend; ++x) {
    for (int z = 0; z < localmesh->LocalNz; ++z) {
      const BoutReal interior = f(x, localmesh->ystart, z);
      ASSERT_GT(interior, 0.0);
      EXPECT_DOUBLE_EQ(f(x, localmesh->ystart - 1, z), -interior);
    }
  }
}

// Test for greater than zero at boundary
void expectNotZeroedAtTarget(const Field3D& f) {
  Mesh* localmesh = f.getMesh();
  for (int x = localmesh->xstart; x <= localmesh->xend; ++x) {
    for (int z = 0; z < localmesh->LocalNz; ++z) {
      ASSERT_GT(f(x, localmesh->ystart, z), 0.0);
      EXPECT_GT(f(x, localmesh->ystart - 1, z), 0.0);
    }
  }
}
} // namespace

// Default behaviour: conductivity and viscosity are zeroed at sheath.
// Combined limiters
TEST_F(NeutralMixedTest, CombinedCoefficientsZeroedAtSheath) {

  Options out = runNeutralMixedTest({{"d",
                                      {{"type", "neutral_mixed"},
                                       {"diagnose", true},
                                       {"AA", 2.0},
                                       {"combined_limiters", true},
                                       {"zero_sheath_conductivity", true},
                                       {"zero_sheath_viscosity", true}}}});

  expectZeroedAtTarget(out["Dnnd"].as<Field3D>());

  expectZeroedAtTarget(out["kappa_d_perp"].as<Field3D>());
  expectZeroedAtTarget(out["kappa_d_par"].as<Field3D>());
  expectZeroedAtTarget(out["kappa_d_unlimited"].as<Field3D>());

  expectZeroedAtTarget(out["eta_d_perp"].as<Field3D>());
  expectZeroedAtTarget(out["eta_d_par"].as<Field3D>());
  expectZeroedAtTarget(out["eta_d_unlimited"].as<Field3D>());
}

// Default behaviour: conductivity and viscosity are zeroed at sheath.
// Separate limiters
TEST_F(NeutralMixedTest, SeparateCoefficientsZeroedAtSheath) {

  Options out = runNeutralMixedTest({{"d",
                                      {{"type", "neutral_mixed"},
                                       {"diagnose", true},
                                       {"AA", 2.0},
                                       {"combined_limiters", false},
                                       {"zero_sheath_conductivity", true},
                                       {"zero_sheath_viscosity", true}}}});

  expectZeroedAtTarget(out["Dnnd"].as<Field3D>());

  expectZeroedAtTarget(out["kappa_d_perp"].as<Field3D>());
  expectZeroedAtTarget(out["kappa_d_par"].as<Field3D>());
  expectZeroedAtTarget(out["kappa_d_unlimited"].as<Field3D>());

  expectZeroedAtTarget(out["eta_d_perp"].as<Field3D>());
  expectZeroedAtTarget(out["eta_d_par"].as<Field3D>());
  expectZeroedAtTarget(out["eta_d_unlimited"].as<Field3D>());
}

// Conductivity and viscosity are not zeroed at sheath.
// Combined limiters
TEST_F(NeutralMixedTest, CombinedCoefficientsNotZeroedAtSheath) {

  Options out = runNeutralMixedTest({{"d",
                                      {{"type", "neutral_mixed"},
                                       {"diagnose", true},
                                       {"AA", 2.0},
                                       {"combined_limiters", true},
                                       {"zero_sheath_conductivity", false},
                                       {"zero_sheath_viscosity", false}}}});

  expectZeroedAtTarget(out["Dnnd"].as<Field3D>());

  expectNotZeroedAtTarget(out["kappa_d_perp"].as<Field3D>());
  expectNotZeroedAtTarget(out["kappa_d_par"].as<Field3D>());
  expectNotZeroedAtTarget(out["kappa_d_unlimited"].as<Field3D>());

  expectNotZeroedAtTarget(out["eta_d_perp"].as<Field3D>());
  expectNotZeroedAtTarget(out["eta_d_par"].as<Field3D>());
  expectNotZeroedAtTarget(out["eta_d_unlimited"].as<Field3D>());
}

// Conductivity and viscosity are not zeroed at sheath.
// Separate limiters
TEST_F(NeutralMixedTest, SeparateCoefficientsNotZeroedAtSheath) {

  Options out = runNeutralMixedTest({{"d",
                                      {{"type", "neutral_mixed"},
                                       {"diagnose", true},
                                       {"AA", 2.0},
                                       {"combined_limiters", false},
                                       {"zero_sheath_conductivity", false},
                                       {"zero_sheath_viscosity", false}}}});

  expectZeroedAtTarget(out["Dnnd"].as<Field3D>());

  expectNotZeroedAtTarget(out["kappa_d_perp"].as<Field3D>());
  expectNotZeroedAtTarget(out["kappa_d_par"].as<Field3D>());
  expectNotZeroedAtTarget(out["kappa_d_unlimited"].as<Field3D>());

  expectNotZeroedAtTarget(out["eta_d_perp"].as<Field3D>());
  expectNotZeroedAtTarget(out["eta_d_par"].as<Field3D>());
  expectNotZeroedAtTarget(out["eta_d_unlimited"].as<Field3D>());
}
