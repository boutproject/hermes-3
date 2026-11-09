#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include "test_extras.hxx" // FakeMesh

#include <bout/bout_types.hxx>
#include <bout/constants.hxx>
#include <bout/coordinates.hxx>
#include <bout/difops.hxx>
#include <bout/field_factory.hxx> // For generating functions
#include <bout/mesh.hxx>

#include "../../include/component.hxx"
#include "../../include/guarded_options.hxx"
#include "../../include/permissions.hxx"
#include "../../include/relax_potential.hxx"

#include <cmath>

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using RelaxPotentialTest = FakeMeshFixture;

namespace {

Options baseOptions() {
  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";
  return {{"units",
           {{"seconds", 1.0},
            {"Tesla", 1.0},
            {"meters", 1.0},
            {"eV", 100.0},
            {"inv_meters_cubed", 1e19}}}};
}

Field3D makeXProfile(BoutReal scale = 1.0) {
  Field3D result{0.0, mesh};
  for (int i = mesh->xstart; i <= mesh->xend; i++) {
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      for (int k = 0; k < mesh->LocalNz; k++) {
        result(i, j, k) = scale * i;
      }
    }
  }
  return result;
}

} // namespace

/// Check a RelaxPotential object can be successfully constructed and
/// registers vorticity and phi variables with the solver
TEST_F(RelaxPotentialTest, CreateComponent) {
  FakeSolver solver;
  Options options = baseOptions();

  RelaxPotential component("test", options, &solver);

  const Options state = solver.getState();

  ASSERT_TRUE(state.isSet("Vort"));
  ASSERT_TRUE(state.isSet("phi1"));
}

/// Check that with `evolve_vorticity = false` only `phi1` is solver-managed,
/// while `Vort` is left as a derived field and is therefore not registered.
TEST_F(RelaxPotentialTest, CreateComponentRelaxMode) {
  FakeSolver solver;
  Options options = baseOptions();
  options["test"]["evolve_vorticity"] = false;

  const RelaxPotential component("test", options, &solver);

  const Options state = solver.getState();

  ASSERT_FALSE(state.isSet("Vort"));
  ASSERT_TRUE(state.isSet("phi1"));
}

/// Check that the RelaxPotential::transform method sets voricity and
/// phi when the component is instantiated with default parameters
TEST_F(RelaxPotentialTest, Transform) {
  FakeSolver solver;
  Options options = baseOptions();

  RelaxPotential component("test", options, &solver);

  Options state;
  component.transform(state);

  ASSERT_TRUE(state["fields"].isSet("vorticity"));
  ASSERT_TRUE(state["fields"].isSet("phi"));
  // TODO: Check that the results are actually correct
}

/// Check that the RelaxPotential::transform method sets voricity and
/// phi when the component is instantiated with phi_boundary_relax = true
TEST_F(RelaxPotentialTest, PhiBoundaryRelax) {
  FakeSolver solver;

  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";
  Options options{{"units",
                   {{"seconds", 1.0},
                    {"Tesla", 1.0},
                    {"meters", 1.0},
                    {"eV", 100.},
                    {"inv_meters_cubed", 1e19}}},
                  {"test", {{"phi_boundary_relax", true}}}};

  RelaxPotential component("test", options, &solver);

  Options state{{"time", 0.1}};
  component.transform(state);

  ASSERT_TRUE(state["fields"].isSet("vorticity"));
  ASSERT_TRUE(state["fields"].isSet("phi"));
  // TODO: Check that the results are actually correct
}

/// Check that the RelaxPotential::transform method sets voricity,
/// phi, and DivJcol when the component is instantiated with
/// collisional_friction = true
TEST_F(RelaxPotentialTest, CollisionalFriction) {
  FakeSolver solver;

  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";
  Options options{{"units",
                   {{"seconds", 1.0},
                    {"Tesla", 1.0},
                    {"meters", 1.0},
                    {"eV", 100.},
                    {"inv_meters_cubed", 1e19}}},
                  {"test", {{"collisional_friction", true}}}};

  RelaxPotential component("test", options, &solver);
  component.declareAllSpecies(SpeciesInformation({"i"}));
  Coordinates* coords = mesh->getCoordinates();
  coords->setBxy(1.0);

  Options state{
      {"time", 0.1},
      {"species",
       {{"i",
         {{"AA", 1.}, {"charge", 1.}, {"density", 1.}, {"collision_frequency", 1.}}}}}};
  component.transform(state);

  ASSERT_TRUE(state["fields"].isSet("vorticity"));
  ASSERT_TRUE(state["fields"].isSet("phi"));
  ASSERT_TRUE(state["fields"].isSet("DivJcol"));
  // TODO: Check that the results are actually correct
}

TEST_F(RelaxPotentialTest, CalculatePihat) {
  FakeSolver solver;
  Options options = baseOptions();
  options["test"]["average_atomic_mass"] = 2.5;

  RelaxPotential component("test", options, &solver);

  Options state{{"species", {{"d+", {{"pressure", 1.2}, {"AA", 2.0}, {"charge", 1.4}}}}}};
  Permissions permissions{{readOnly("species")}};
  const GuardedOptions guarded_state{&state, &permissions};

  const Field3D Pi_hat = component.calculatePihat(guarded_state["species"]);

  ASSERT_TRUE(IsFieldEqual(Pi_hat, 1.2 * 2.0 / 2.5 / 1.4, "RGN_NOBNDRY"));
}

/// Check that `applyPhiBoundary()` uses the sheath-potential formula when
/// `sheath_boundary = true`, so the X-guard cells are set from electron temperature.
TEST_F(RelaxPotentialTest, ApplyPhiBoundarySheath) {
  FakeSolver solver;
  Options options = baseOptions();
  options["test"]["sheath_boundary"] = true;

  RelaxPotential component("test", options, &solver);

  Field3D phi{1.0};
  Options state{
      {"species", {{"e", {{"AA", 1.0}, {"charge", -1.0}, {"temperature", 2.0}}}}}};
  Permissions permissions{{readOnly("species")}};
  const GuardedOptions guarded_state{&state, &permissions};

  component.applyPhiBoundary(phi, guarded_state);

  const BoutReal pi = std::acos(-1.0);
  const BoutReal sheathmult = std::log(0.5 * std::sqrt(1.0 / pi));
  const BoutReal expected = (2.0 * (2.0 * sheathmult)) - 1.0;

  ASSERT_NEAR(phi(mesh->xstart - 1, mesh->ystart, 0), expected, 1e-12);
  ASSERT_NEAR(phi(mesh->xend + 1, mesh->ystart, 0), expected, 1e-12);
}

/// Check that `applyPhiBoundary()` relaxes the X-boundary exponentially toward
/// zero-gradient when `phi_boundary_relax = true`, using the configured timescale.
TEST_F(RelaxPotentialTest, ApplyPhiBoundaryRelax) {
  FakeSolver solver;
  Options options = baseOptions();
  options["test"]["phi_boundary_relax"] = true;
  options["test"]["phi_boundary_timescale"] = 1.0;

  RelaxPotential component("test", options, &solver);

  Field3D phi{2.0};
  for (int j = mesh->ystart; j <= mesh->yend; j++) {
    for (int k = 0; k < mesh->LocalNz; k++) {
      phi(mesh->xstart - 1, j, k) = 0.0;
      phi(mesh->xend + 1, j, k) = 0.0;
    }
  }

  Options state;
  state["time"].assign(0.0, "initial");
  Permissions permissions{{readOnly("time")}};
  const GuardedOptions guarded_state{&state, &permissions};

  component.applyPhiBoundary(phi, guarded_state);
  Options state2;
  state2["time"].assign(1.0, "next");
  const GuardedOptions guarded_state2{&state2, &permissions};
  component.applyPhiBoundary(phi, guarded_state2);

  const BoutReal expected = 2.0 - (2.0 * std::exp(-1.0));
  ASSERT_NEAR(phi(mesh->xstart - 1, mesh->ystart, 0), expected, 1e-12);
  ASSERT_NEAR(phi(mesh->xend + 1, mesh->ystart, 0), expected, 1e-12);
}

/// Check that `calculateDivJdia()` returns zero for uniform-zero `phi`, and that
/// the corresponding `energy_source` update is also zero in that case.
TEST_F(RelaxPotentialTest, CalculateDivJdia) {
  FakeSolver solver;
  Options options = baseOptions();

  RelaxPotential component("test", options, &solver);

  Options state{{"species", {{"d+", {{"pressure", 1.2}, {"AA", 2.0}, {"charge", 1.4}}}}}};
  Permissions permissions{{readOnly("species"), readWrite("species:d+:energy_source")}};
  const GuardedOptions guarded_state{&state, &permissions};

  Field3D phi{0.0}; // Y boundary modified in calculateDivJdia
  const Field3D DivJdia = component.calculateDivJdia(phi, guarded_state["species"]);

  ASSERT_TRUE(IsFieldEqual(DivJdia, 0.0, "RGN_NOBNDRY"));
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["d+"]["energy_source"]), 0.0,
                           "RGN_NOBNDRY"));
}

/// Check that `calculateDivJcol()` returns zero when `phi + pi_hat` is spatially
/// uniform, because there is then no perpendicular gradient to drive the term.
TEST_F(RelaxPotentialTest, CalculateDivJcol) {
  FakeSolver solver;
  Options options = baseOptions();

  RelaxPotential component("test", options, &solver);

  Options state{{"species",
                 {{"d+",
                   {{"density", 3.0},
                    {"collision_frequency", 0.4},
                    {"AA", 2.0},
                    {"charge", 1.0}}}}}};
  Permissions permissions{{readOnly("species")}};
  const GuardedOptions guarded_state{&state, &permissions};

  const Field3D DivJcol = component.calculateDivJcol(1.0, 0.5, guarded_state["species"]);

  ASSERT_TRUE(IsFieldEqual(DivJcol, 0.0, "RGN_NOBNDRY"));
}

/// Check that `calculateParallelCurrentSource()` forwards an externally supplied
/// `fields:DivJextra` contribution unchanged into the returned vorticity source.
TEST_F(RelaxPotentialTest, CalculateParallelCurrentSourceIncludesDivJextra) {
  FakeSolver solver;
  Options options = baseOptions();

  RelaxPotential component("test", options, &solver);

  const Options state{{"fields", {{"DivJextra", 1.25}}}};
  const Field3D source = component.calculateParallelCurrentSource(state);

  ASSERT_TRUE(IsFieldEqual(source, 1.25, "RGN_NOBNDRY"));
}

/// Check the species parallel current reaches the vorticity source, not only DivJextra
TEST_F(RelaxPotentialTest, CalculateParallelCurrentSourceIncludesSpeciesCurrent) {
  FakeSolver solver;
  Options options = baseOptions();

  RelaxPotential component("test", options, &solver);

  Field3D NV{0.0, mesh};
  for (int i = 0; i < mesh->LocalNx; i++) {
    for (int j = 0; j < mesh->LocalNy; j++) {
      for (int k = 0; k < mesh->LocalNz; k++) {
        NV(i, j, k) = 0.1 * j * j; // varies along y, so Div_par isn't zero
      }
    }
  }
  // fields has to be there (it always is in a run, phi lives in it)
  const Options state{
      {"fields", {{"DivJextra", 1.25}}},
      {"species", {{"d+", {{"charge", 1.0}, {"AA", 2.0}, {"momentum", NV}}}}}};
  const Field3D source = component.calculateParallelCurrentSource(state);

  const Field3D div_jpar = Div_par(Field3D{0.5 * NV}); // jpar = (Z / A) NV
  ASSERT_GT(max(abs(div_jpar), false, "RGN_NOBNDRY"), 0.0);
  ASSERT_TRUE(IsFieldEqual(source, Field3D{1.25 + div_jpar}, "RGN_NOBNDRY"));
}

/// Check that in evolving-vorticity mode `calculatePhi1Source()` depends only on
/// `vort_from_phi - vort`, and is therefore independent of the supplied `vort_rhs`.
TEST_F(RelaxPotentialTest, CalculatePhi1SourceEvolveVorticityUsesVortDifference) {
  FakeSolver solver;
  Options options = baseOptions();

  RelaxPotential component("test", options, &solver);

  const Field3D source1 = component.calculatePhi1Source(2.0, 9.0, 5.0);
  const Field3D source2 = component.calculatePhi1Source(2.0, 0.0, 5.0);

  ASSERT_TRUE(IsFieldEqual(source1, source2, "RGN_NOBNDRY"));
}

/// Check that in relax-only mode `calculatePhi1Source()` depends on `vort_rhs`
/// alone, and not on the separate `vort` or `vort_from_phi` arguments.
TEST_F(RelaxPotentialTest, CalculatePhi1SourceRelaxModeUsesVortRhs) {
  FakeSolver solver;
  Options options = baseOptions();
  options["test"]["evolve_vorticity"] = false;

  const RelaxPotential component("test", options, &solver);

  const Field3D source1 = component.calculatePhi1Source(2.0, 9.0, 5.0);
  const Field3D source2 = component.calculatePhi1Source(8.0, 9.0, 1.0);
  const Field3D source3 = component.calculatePhi1Source(2.0, 3.0, 5.0);

  ASSERT_TRUE(IsFieldEqual(source1, source2, "RGN_NOBNDRY"));
  ASSERT_FALSE(IsFieldEqual(source1, source3, "RGN_NOBNDRY"));
}

/// Check the Boussinesq diamagnetic-polarisation regression: species with the same
/// `A/Z * pressure` contribution should produce the same vorticity contribution.
TEST_F(RelaxPotentialTest, BoussinesqVorticityUsesChargeWeighting) {
  FakeSolver solver;
  Options options = baseOptions();
  options["test"]["diamagnetic"] = false;

  RelaxPotential component("test", options, &solver);

  const Field3D phi{0.0, mesh};
  const Field3D pressure = makeXProfile();

  Options state_charge_two{
      {"species", {{"d+", {{"pressure", pressure}, {"AA", 2.0}, {"charge", 2.0}}}}}};
  Options state_charge_one{
      {"species",
       {{"d+", {{"pressure", 0.5 * pressure}, {"AA", 2.0}, {"charge", 1.0}}}}}};
  Permissions permissions{{readOnly("species")}};
  const GuardedOptions guarded_charge_two{&state_charge_two, &permissions};
  const GuardedOptions guarded_charge_one{&state_charge_one, &permissions};
  auto species_charge_two = guarded_charge_two["species"];
  auto species_charge_one = guarded_charge_one["species"];

  const Field3D vort_charge_two = component.vorticity(phi, species_charge_two);
  const Field3D vort_charge_one = component.vorticity(phi, species_charge_one);

  ASSERT_TRUE(IsFieldEqual(vort_charge_two, vort_charge_one, "RGN_NOBNDRY"));
}

/// Check the boundary-ordering regression fixed in `transform()`: `finally()` should
/// use `Vort_from_phi` computed from boundary-conditioned `phi`, not from pre-boundary `phi`.
/// The expected value is reconstructed by mirroring the `transform()` sequence:
/// apply the `phi1` boundary condition, convert to `phi`, apply the configured
/// radial `phi` boundary, then calculate the vorticity from that field.
TEST_F(RelaxPotentialTest, FinallyUsesBoundaryConditionedPhiVorticity) {
  FakeSolver solver;
  Options options = baseOptions();
  options["test"]["diamagnetic"] = false;
  options["test"]["exb_advection"] = false;
  options["test"]["phi_dissipation"] = false;
  options["test"]["lambda_1"] = 1.0;
  options["test"]["lambda_2"] = 1.0;
  options["test"]["sheath_boundary"] = true;

  RelaxPotential component("test", options, &solver);

  solver.setState({{"Vort", Field3D(0.0, mesh)}, {"phi1", Field3D(1.0, mesh)}});

  Options state{
      {"species", {{"e", {{"AA", 1.0}, {"charge", -1.0}, {"temperature", 2.0}}}}}};
  component.declareAllSpecies(SpeciesInformation({"e"}));
  component.transform(state);
  component.finally(state);

  Permissions permissions{{readOnly("species")}};
  const GuardedOptions guarded_state{&state, &permissions};
  auto guarded_species = guarded_state["species"];

  Field3D expected_phi1{1.0, mesh};
  expected_phi1.applyBoundary("neumann");
  Field3D expected_phi = expected_phi1;
  component.applyPhiBoundary(expected_phi, guarded_state);
  const Field3D expected = component.calculatePhi1Source(
      0.0, 0.0, component.vorticity(expected_phi, guarded_species));
  const Options ddt = solver.getTimeDerivs();

  ASSERT_TRUE(IsFieldEqual(ddt["phi1"].as<Field3D>(), expected, "RGN_NOBNDRY"));
}

/// Check that diagnostic output reports the physical potential derivative
/// `ddt(phi) = ddt(phi1) / lambda_2`, while still outputting the derived `Vort` field.
TEST_F(RelaxPotentialTest, OutputVarsScalesPhysicalPhiDerivative) {
  FakeSolver solver;
  Options options = baseOptions();
  options["test"]["evolve_vorticity"] = false;
  options["test"]["diamagnetic"] = false;
  options["test"]["exb_advection"] = false;
  options["test"]["phi_dissipation"] = false;
  options["test"]["phi_boundary_relax"] = true;
  options["test"]["lambda_1"] = 2.0;
  options["test"]["lambda_2"] = 4.0;
  options["test"]["diagnose"] = true;

  RelaxPotential component("test", options, &solver);

  Options state{{"time", 0.0}, {"fields", {{"DivJextra", 1.25}}}};
  component.transform(state);
  component.finally(state);

  const Options ddt = solver.getTimeDerivs();
  Options outputs{{"Nnorm", 1.0}, {"Tnorm", 1.0}, {"Omega_ci", 1.0}};
  component.outputVars(outputs);

  ASSERT_TRUE(outputs.isSet("Vort"));
  ASSERT_TRUE(IsFieldEqual(outputs["ddt(phi)"].as<Field3D>(),
                           Field3D(ddt["phi1"].as<Field3D>() / 4.0), "RGN_NOBNDRY"));
}

/// Check that `finally()` adds the extracted source terms into `ddt(Vort)`:
/// with nonlinear and dissipation terms disabled, `DivJextra` should appear unchanged.
TEST_F(RelaxPotentialTest, FinallyUsesExtractedSources) {
  FakeSolver solver;
  Options options = baseOptions();
  options["test"]["exb_advection"] = false;
  options["test"]["phi_dissipation"] = false;

  RelaxPotential component("test", options, &solver);

  Options state{{"fields", {{"DivJextra", 1.25}}}};
  component.transform(state);
  component.finally(state);

  Options ddt = solver.getTimeDerivs();
  ASSERT_TRUE(IsFieldEqual(ddt["Vort"].as<Field3D>(), 1.25, "RGN_NOBNDRY"));
}

namespace {
// the ExB operator reaches two cells out in x, and it needs more than one interior x
// cell (phi is zero on both x faces, so with one cell every ExB velocity vanishes)
using RelaxPotentialChannelsTest = FakeMeshFixture_tmpl<9, 7, 8, false, 2>;

constexpr BoutReal me_mp = 1. / 1836;

// every term in Vort_rhs switched on, so each channel has something in it
Options channelOptions(bool evolve_vorticity, bool diagnose) {
  Options options = baseOptions();
  auto& opt = options["test"];
  opt["evolve_vorticity"] = evolve_vorticity;
  opt["diagnose"] = diagnose;
  opt["collisional_friction"] = true;
  // the poloidal ExB part reads DDY(phi) in corner guard cells that are never set
  opt["poloidal_flows"] = false;
  opt["viscosity"] = 0.1;
  opt["viscosity_par"] = 0.05;
  opt["vort_dissipation"] = true;
  opt["hyper_z"] = 0.01;
  opt["lambda_1"] = 1.0;
  opt["lambda_2"] = 2.0;
  return options;
}

// curvature that varies in x and y, so Div(Jdia) isn't zero
void setCurvature() {
  dynamic_cast<FakeMesh*>(mesh)->setGridDataSource(
      new FakeGridDataSource{{{"bxcvx", "0.1 + 0.05 * y"},
                              {"bxcvy", "0.02 * x"},
                              {"bxcvz", "0.3 + 0.1 * x * y"}}});
}

Field3D make(const std::string& expr) {
  Options field_options;
  return FieldFactory::get()->create3D(expr, &field_options, mesh);
}

Options channelState() {
  return Options{
      {"sound_speed", make("1 + 0.1 * x + 0.05 * y")},
      {"fields", {{"DivJextra", make("0.01 * (1 + x * y)")}}},
      {"species",
       {{"d+",
         {{"AA", 2.0},
          {"charge", 1.0},
          {"density", make("1 + 0.1 * x + 0.2 * y + 0.05 * sin(z)")},
          {"pressure", make("2 + 0.3 * x * y + 0.1 * cos(z)")},
          {"momentum", make("0.6 + 0.2 * y + 0.1 * x * cos(z)")},
          {"collision_frequency", make("0.5 + 0.1 * y")}}},
        {"e",
         {{"AA", me_mp},
          {"charge", -1.0},
          {"density", make("1 + 0.1 * x + 0.2 * y + 0.05 * sin(z)")},
          {"pressure", make("1.5 + 0.2 * x + 0.1 * y * y")},
          {"temperature", make("1.2 + 0.1 * y")},
          {"momentum", make("-0.002 + 0.001 * x * y")}}}}}};
}

struct ChannelRun {
  Options state;
  Options outputs;
  Options ddt;
};

ChannelRun runChannels(bool evolve_vorticity, bool diagnose) {
  setCurvature();
  FakeSolver solver;
  Options options = channelOptions(evolve_vorticity, diagnose);
  RelaxPotential component("test", options, &solver);
  component.declareAllSpecies(SpeciesInformation({"d+", "e"}));

  Options init{{"phi1", make("0.5 * x * x + 0.3 * sin(z) + 0.4 * y * (1 + 0.1 * x)")}};
  if (evolve_vorticity) {
    init["Vort"] = make("0.2 + 0.1 * x * y + 0.05 * cos(z)");
  }
  solver.setState(init);

  ChannelRun run{channelState(), {{"Nnorm", 1.0}, {"Tnorm", 1.0}, {"Omega_ci", 1.0}}, {}};
  component.transform(run.state);
  component.finally(run.state);
  component.outputVars(run.outputs);
  run.ddt = solver.getTimeDerivs();
  return run;
}

// Pi_hat as calculatePihat builds it (average_atomic_mass is 2)
Field3D piHat(const Options& state) {
  return get<Field3D>(state["species"]["d+"]["pressure"]) * (2.0 / 2.0 / 1.0)
         + get<Field3D>(state["species"]["e"]["pressure"]) * (me_mp / 2.0 / -1.0);
}

BoutReal scaleOf(const Field3D& f) { return max(abs(f), false, "RGN_NOBNDRY"); }

const std::vector<std::string> species_names{"d+", "e"};
} // namespace

/// Closure D: each channel is its weight times the exact Vort_rhs term, the unweighted
/// terms sum to Vort_rhs and the weighted ones to Evort_relax, in both modes
TEST_F(RelaxPotentialChannelsTest, ChannelsCloseVorticity) {
  for (const bool evolve_vorticity : {true, false}) {
    SCOPED_TRACE(evolve_vorticity ? "evolve_vorticity" : "relax only");
    setCurvature();
    FakeSolver solver;
    Options options = channelOptions(evolve_vorticity, true);
    RelaxPotential component("test", options, &solver);
    component.declareAllSpecies(SpeciesInformation({"d+", "e"}));
    Options init{{"phi1", make("0.5 * x * x + 0.3 * sin(z) + 0.4 * y * (1 + 0.1 * x)")}};
    if (evolve_vorticity) {
      init["Vort"] = make("0.2 + 0.1 * x * y + 0.05 * cos(z)");
    }
    solver.setState(init);

    Options state = channelState();
    component.transform(state);
    component.finally(state);
    Options outputs{{"Nnorm", 1.0}, {"Tnorm", 1.0}, {"Omega_ci", 1.0}};
    component.outputVars(outputs);
    const Options ddt = solver.getTimeDerivs();

    const Field3D phi = get<Field3D>(outputs["phi"]);
    const Field3D vort = get<Field3D>(outputs["Vort"]);
    const Field3D vort_rhs = get<Field3D>(outputs["ddt(Vort)"]);
    const Field3D pi_hat = piHat(state);
    const Field3D weight = -(phi + pi_hat);

    // the unweighted terms, taken from the outputs or the component's own functions
    const Field3D div_jdia = get<Field3D>(outputs["DivJdia"]);
    const Field3D div_jcol = get<Field3D>(outputs["DivJcol"]);
    const Field3D exb = component.calculateExBAdvectionSource(vort, phi, pi_hat);
    const Field3D par = component.calculateParallelCurrentSource(state);
    const Field3D diss = component.calculateDissipationSource(state, vort, phi);
    const Field3D div_jextra = get<Field3D>(state["fields"]["DivJextra"]);

    const BoutReal rhs_scale = scaleOf(vort_rhs);
    ASSERT_GT(rhs_scale, 0.0);
    for (const Field3D& term : {div_jdia, div_jcol, exb, par, diss}) {
      EXPECT_GT(scaleOf(term), 0.0);
    }
    EXPECT_TRUE(IsFieldEqual(Field3D{div_jdia + div_jcol + exb + par + diss}, vort_rhs,
                             "RGN_NOBNDRY", 1e-12 * rhs_scale));

    // field level channels are the weight times exactly that term
    auto expect_weighted = [&](const std::string& name, const Field3D& expected) {
      const Field3D channel = get<Field3D>(outputs[name]);
      EXPECT_GT(scaleOf(channel), 0.0) << name;
      EXPECT_TRUE(IsFieldEqual(channel, expected, "RGN_NOBNDRY", 1e-12 * scaleOf(expected)))
          << name;
    };
    expect_weighted("Evort_DivJcol", weight * div_jcol);
    expect_weighted("Evort_ExB", weight * exb);
    expect_weighted("Evort_numerical", weight * diss);
    expect_weighted("Evort_relax", weight * vort_rhs);

    // per species: the phi part of each current, and Jdia.Grad(phi) is all of energy_source
    Field3D sum_phi_jdia{zeroFrom(phi)};
    Field3D sum_phi_jpar{zeroFrom(phi)};
    Field3D sum_div_jpar{zeroFrom(phi)};
    for (const auto& sp : species_names) {
      SCOPED_TRACE(sp);
      const Options& species = state["species"][sp];
      const Field3D jpar = (get<BoutReal>(species["charge"]) / get<BoutReal>(species["AA"]))
                           * get<Field3D>(species["momentum"]);
      const Field3D div_jpar = Div_par(jpar);
      sum_div_jpar += div_jpar;
      expect_weighted("E" + sp + "_phi_Div_jpar", -phi * div_jpar);
      sum_phi_jpar += get<Field3D>(outputs["E" + sp + "_phi_Div_jpar"]);
      sum_phi_jdia += get<Field3D>(outputs["E" + sp + "_phi_Div_Jdia"]);

      const Field3D jdia_grad_phi = get<Field3D>(outputs["E" + sp + "_Jdia_Grad_phi"]);
      EXPECT_GT(scaleOf(jdia_grad_phi), 0.0);
      EXPECT_TRUE(IsFieldEqual(jdia_grad_phi, get<Field3D>(species["energy_source"]),
                               "RGN_NOBNDRY", 0.0));
      EXPECT_EQ(sp, outputs["E" + sp + "_phi_Div_Jdia"].attributes["species"].as<std::string>());
    }
    // the per species split of Div(Jdia) only matches the sum to rounding
    expect_weighted("Evort_Pi_DivJ", -pi_hat * (div_jdia + sum_div_jpar));
    EXPECT_TRUE(IsFieldEqual(sum_phi_jdia, Field3D{-phi * div_jdia}, "RGN_NOBNDRY",
                             1e-12 * scaleOf(sum_phi_jdia)));

    // weighted terms add up to the residual (DivJextra has no channel of its own)
    const Field3D weighted_sum =
        sum_phi_jdia + sum_phi_jpar + get<Field3D>(outputs["Evort_DivJcol"])
        + get<Field3D>(outputs["Evort_ExB"]) + get<Field3D>(outputs["Evort_numerical"])
        + get<Field3D>(outputs["Evort_Pi_DivJ"]) + weight * div_jextra;
    const Field3D relax = get<Field3D>(outputs["Evort_relax"]);
    EXPECT_TRUE(IsFieldEqual(weighted_sum, relax, "RGN_NOBNDRY", 1e-12 * scaleOf(relax)));

    // Vort_mismatch is what drives phi1 when Vort evolves, and exactly zero when it doesn't
    const BoutReal lambda_1 = 1.0 / (100.0 * 1.0 / SI::qe / 1e19); // as normalised
    const Field3D ddt_phi1 = ddt["phi1"].as<Field3D>();
    const Field3D mismatch = get<Field3D>(outputs["Vort_mismatch"]);
    if (evolve_vorticity) {
      EXPECT_GT(scaleOf(mismatch), 0.0);
      EXPECT_TRUE(IsFieldEqual(Field3D{lambda_1 * mismatch}, ddt_phi1, "RGN_NOBNDRY",
                               1e-12 * scaleOf(ddt_phi1)));
    } else {
      EXPECT_TRUE(IsFieldEqual(mismatch, 0.0, "RGN_NOBNDRY", 0.0));
      // the residual is the weight on what drives phi1
      EXPECT_TRUE(IsFieldEqual(relax, Field3D{(phi + pi_hat) * ddt_phi1 / lambda_1},
                               "RGN_NOBNDRY", 1e-12 * scaleOf(relax)));
    }
  }
}

/// Channels don't touch ddt or energy_source, and nothing new is saved when off
TEST_F(RelaxPotentialChannelsTest, DiagnoseDoesNotChangeDdt) {
  for (const bool evolve_vorticity : {true, false}) {
    SCOPED_TRACE(evolve_vorticity ? "evolve_vorticity" : "relax only");
    const ChannelRun off = runChannels(evolve_vorticity, false);
    const ChannelRun on = runChannels(evolve_vorticity, true);

    EXPECT_TRUE(IsFieldEqual(off.ddt["phi1"].as<Field3D>(), on.ddt["phi1"].as<Field3D>(),
                             "RGN_NOBNDRY", 0.0));
    if (evolve_vorticity) {
      EXPECT_TRUE(IsFieldEqual(off.ddt["Vort"].as<Field3D>(), on.ddt["Vort"].as<Field3D>(),
                               "RGN_NOBNDRY", 0.0));
    }
    for (const auto& sp : species_names) {
      EXPECT_TRUE(IsFieldEqual(get<Field3D>(off.state["species"][sp]["energy_source"]),
                               get<Field3D>(on.state["species"][sp]["energy_source"]),
                               "RGN_NOBNDRY", 0.0))
          << sp;
    }
    for (const auto& name : {"Evort_relax", "Evort_Pi_DivJ", "Vort_mismatch",
                             "Ed+_Jdia_Grad_phi", "Ee_phi_Div_Jdia"}) {
      EXPECT_FALSE(off.outputs.isSet(name)) << name;
      EXPECT_TRUE(on.outputs.isSet(name)) << name;
    }
  }
}
