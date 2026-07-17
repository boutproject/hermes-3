
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

// Function to test cross-field diffusion in presence of a radial pressure gradient.
namespace {
auto runDnnTest(Options options, bool with_collisions = false) {
  FakeSolver solver;

  options["units"] = {
      {"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}, {"meters", 1.0}};
  NeutralMixed component("d", options, &solver);

  // Make pressure gradient in X direction
  Field3D Pn = makeField<Field3D>([](Ind3D& i) { return 1.0 + 1 * i.x(); }, mesh);

  // Call the finally() method with a density, energy, and momentum source
  Options state = {{"species",
                    {{"d",
                      {{"density", 1.0},
                       {"density_source", 0.5},
                       {"pressure", Pn},
                       {"energy_source", 1.5},
                       {"momentum", 1.0},
                       {"momentum_source", 0.75},
                       {"temperature", 1.0},
                       {"velocity", 1.0}}}}}};

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

  Field3D Dnn = out["Dnnd"].as<Field3D>();
  Field3D Dunl = out["Dnnd_unlimited"].as<Field3D>();
  Field3D Dmax = out["Dnnd_max"].as<Field3D>();

  return std::make_tuple(Dnn, Dunl, Dmax);
}
} // namespace

// Check that flux limiter reduces to a simple harmonic mean when sharpness = 1.0.
TEST_F(NeutralMixedTest, DnnHarmonicLimiter) {

  auto [Dnn, Dunl, Dmax] = runDnnTest({{"d",
                                        {
                                            {"type", "neutral_mixed"},
                                            {"diagnose", true},
                                            {"AA", 2.0},
                                            {"flux_limiter_sharpness", 1.0},
                                        }}});

  BOUT_FOR_SERIAL(i, Dnn.getRegion("RGN_NOBNDRY")) {
    EXPECT_DOUBLE_EQ(Dnn[i], Dunl[i] * Dmax[i] / (Dunl[i] + Dmax[i]));
  }
}

// Check that an aggressive flux limit limits the flux.
TEST_F(NeutralMixedTest, DnnTightLimit) {

  auto [Dnn, Dunl, Dmax] = runDnnTest({{"d",
                                        {{"type", "neutral_mixed"},
                                         {"diagnose", true},
                                         {"AA", 2.0},
                                         {"flux_limit", 1e-5}}}});

  BOUT_FOR_SERIAL(i, Dnn.getRegion("RGN_NOBNDRY")) {
    EXPECT_LT(Dmax[i], Dunl[i]);
    EXPECT_LT(Dnn[i], Dunl[i]);
    EXPECT_NEAR(Dnn[i], Dmax[i], Dnn[i] * 1e-3);
  }
}

// Check that a loose flux limit doesn't limit the flux.
TEST_F(NeutralMixedTest, DnnLooseLimit) {

  auto [Dnn, Dunl, Dmax] = runDnnTest({{"d",
                                        {{"type", "neutral_mixed"},
                                         {"diagnose", true},
                                         {"AA", 2.0},
                                         {"flux_limit", 1000}}}});

  BOUT_FOR_SERIAL(i, Dnn.getRegion("RGN_NOBNDRY")) {
    EXPECT_GT(Dmax[i], Dunl[i]);
    EXPECT_NEAR(Dnn[i], Dunl[i], Dnn[i] * 1e-3);
  }
}

// Check that the explicit diffusion limit can override the flux limitation.
TEST_F(NeutralMixedTest, DnnExplicitLimit) {

  auto [Dnn, Dunl, Dmax] = runDnnTest({{"d",
                                        {{"type", "neutral_mixed"},
                                         {"diagnose", true},
                                         {"AA", 2.0},
                                         {"flux_limit", 0.2},
                                         {"diffusion_limit", 1e-5}}}});

  BOUT_FOR_SERIAL(i, Dnn.getRegion("RGN_NOBNDRY")) {
    EXPECT_DOUBLE_EQ(Dmax[i], 1e-5);
    EXPECT_NEAR(Dnn[i], 1e-5, 1e-3 * Dnn[i]);
  }
}

// Check that adding collisionality reduces Dnn.
TEST_F(NeutralMixedTest, DnnCollisionalityImpact) {

  Options options = {{"d", {{"type", "neutral_mixed"}, {"diagnose", true}, {"AA", 2.0}}}};

  auto [Dnn, Dunl, Dmax] = runDnnTest(options.copy(), false);
  auto [Dnn_coll, Dunl_coll, Dmax_coll] = runDnnTest(options.copy(), true);

  BOUT_FOR_SERIAL(i, Dnn.getRegion("RGN_NOBNDRY")) { EXPECT_LT(Dnn_coll[i], Dnn[i]); }
}

// Check that reducing neutral_lmax raises collisionality floor.
// Lower neutral_lmax leads to a higher nu_pseudo_mfp,
// which increases total nu. and reduces Dnn.
TEST_F(NeutralMixedTest, DnnCollisionalityFloor) {

  auto [Dnn_lo_lmax, Dunl_lo_lmax, Dmax_lo_lmax] =
      runDnnTest({{"d",
                   {{"type", "neutral_mixed"},
                    {"diagnose", true},
                    {"AA", 2.0},
                    {"neutral_lmax", 0.01}}}});

  auto [Dnn_hi_lmax, Dunl_hi_lmax, Dmax_hi_lmax] =
      runDnnTest({{"d",
                   {{"type", "neutral_mixed"},
                    {"diagnose", true},
                    {"AA", 2.0},
                    {"neutral_lmax", 100}}}});

  BOUT_FOR_SERIAL(i, Dnn_lo_lmax.getRegion("RGN_NOBNDRY")) {
    EXPECT_LT(Dnn_lo_lmax[i], Dnn_hi_lmax[i]);
  }
}
