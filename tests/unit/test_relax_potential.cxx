#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include "test_extras.hxx" // FakeMesh

#include <bout/bout_types.hxx>
#include <bout/coordinates.hxx>
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
