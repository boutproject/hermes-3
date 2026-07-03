#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include "test_extras.hxx" // FakeMesh

#include <bout/bout_types.hxx>
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

} // namespace

TEST_F(RelaxPotentialTest, CreateComponent) {
  FakeSolver solver;
  Options options = baseOptions();

  const RelaxPotential component("test", options, &solver);

  const Options state = solver.getState();

  ASSERT_TRUE(state.isSet("Vort"));
  ASSERT_TRUE(state.isSet("phi1"));
}

TEST_F(RelaxPotentialTest, Transform) {
  FakeSolver solver;
  Options options = baseOptions();

  RelaxPotential component("test", options, &solver);

  Options state;
  component.transform(state);

  ASSERT_TRUE(state["fields"].isSet("vorticity"));
  ASSERT_TRUE(state["fields"].isSet("phi"));
}

TEST_F(RelaxPotentialTest, CalculatePihat) {
  FakeSolver solver;
  Options options = baseOptions();
  options["test"]["average_atomic_mass"] = 2.5;

  RelaxPotential component("test", options, &solver);

  Options state{{"species", {{"d+", {{"pressure", 1.2}, {"AA", 2.0}, {"charge", 1.4}}}}}};
  Permissions permissions{{readOnly("species")}};
  GuardedOptions guarded_state{&state, &permissions};

  const Field3D Pi_hat = component.calculatePihat(guarded_state["species"]);

  ASSERT_TRUE(IsFieldEqual(Pi_hat, 1.2 * 2.0 / 2.5 / 1.4, "RGN_NOBNDRY"));
}

TEST_F(RelaxPotentialTest, ApplyPhiBoundarySheath) {
  FakeSolver solver;
  Options options = baseOptions();
  options["test"]["sheath_boundary"] = true;

  RelaxPotential component("test", options, &solver);

  Field3D phi{1.0};
  Options state{{"species", {{"e", {{"AA", 1.0}, {"temperature", 2.0}}}}}};
  Permissions permissions{{readOnly("species")}};
  GuardedOptions guarded_state{&state, &permissions};

  component.applyPhiBoundary(phi, guarded_state);

  const BoutReal pi = std::acos(-1.0);
  const BoutReal sheathmult = std::log(0.5 * std::sqrt(1.0 / pi));
  const BoutReal expected = 2.0 * (2.0 * sheathmult) - 1.0;

  ASSERT_NEAR(phi(mesh->xstart - 1, mesh->ystart, 0), expected, 1e-12);
  ASSERT_NEAR(phi(mesh->xend + 1, mesh->ystart, 0), expected, 1e-12);
}

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

TEST_F(RelaxPotentialTest, CalculateDivJdia) {
  FakeSolver solver;
  Options options = baseOptions();

  RelaxPotential component("test", options, &solver);

  Options state{{"species", {{"d+", {{"pressure", 1.2}, {"AA", 2.0}, {"charge", 1.4}}}}}};
  Permissions permissions{{readOnly("species"), readWrite("species:d+:energy_source")}};
  GuardedOptions guarded_state{&state, &permissions};

  Field3D phi{0.0}; // Y boundary modified in calculateDivJdia
  const Field3D DivJdia = component.calculateDivJdia(phi, guarded_state["species"]);

  ASSERT_TRUE(IsFieldEqual(DivJdia, 0.0, "RGN_NOBNDRY"));
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["d+"]["energy_source"]), 0.0,
                           "RGN_NOBNDRY"));
}

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

TEST_F(RelaxPotentialTest, CalculateParallelCurrentSourceIncludesDivJextra) {
  FakeSolver solver;
  Options options = baseOptions();

  RelaxPotential component("test", options, &solver);

  const Options state{{"fields", {{"DivJextra", 1.25}}}};
  const Field3D source = component.calculateParallelCurrentSource(state);

  ASSERT_TRUE(IsFieldEqual(source, 1.25, "RGN_NOBNDRY"));
}

TEST_F(RelaxPotentialTest, CalculatePhi1SourceEvolveVorticityUsesVortDifference) {
  FakeSolver solver;
  Options options = baseOptions();

  const RelaxPotential component("test", options, &solver);

  const Field3D source1 = component.calculatePhi1Source(2.0, 9.0, 5.0);
  const Field3D source2 = component.calculatePhi1Source(2.0, 0.0, 5.0);

  ASSERT_TRUE(IsFieldEqual(source1, source2, "RGN_NOBNDRY"));
}

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
