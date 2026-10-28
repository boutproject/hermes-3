#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/vorticity.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using VorticityTest = FakeMeshFixture;

TEST_F(VorticityTest, CreateComponent) {
  FakeSolver solver;

  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";
  Options options{{"units", {{"seconds", 1.0}, {"Tesla", 1.0}, {"meters", 1.0}}}};

  Vorticity component("test", options, &solver);

  Options state = solver.getState();

  ASSERT_TRUE(state.isSet("Vort"));
}

TEST_F(VorticityTest, Transform) {
  FakeSolver solver;

  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";
  Options options{{"units", {{"seconds", 1.0}, {"Tesla", 1.0}, {"meters", 1.0}}}};

  Vorticity component("test", options, &solver);

  Options state;
  component.transform(state);

  ASSERT_TRUE(state["fields"].isSet("vorticity"));
  ASSERT_TRUE(state["fields"].isSet("phi"));
}

TEST_F(VorticityTest, TransformNoDiamagnetic) {
  FakeSolver solver;
  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";

  Options options{{"units", {{"seconds", 1.0}, {"Tesla", 1.0}, {"meters", 1.0}}},
                  {"test", {{"diamagnetic", false}}}};

  Vorticity component("test", options, &solver);
  component.declareAllSpecies({{"d+"}, {}});

  Options state{{"species", {{"d+", {{"AA", 2.0}, {"charge", 1.0}, {"pressure", 1.0}}}}}};

  component.transform(state);

  // No diamagnetic terms so values are not set
  ASSERT_FALSE(state["species"]["d+"].isSet("energy_source"));
  ASSERT_FALSE(state["fields"].isSet("DivJdia"));
}

TEST_F(VorticityTest, TransformWithCurlb_B) {
  FakeSolver solver;
  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";

  // Add curvature data
  static_cast<FakeMesh*>(bout::globals::mesh)
      ->setGridDataSource(
          new FakeGridDataSource{{{"bxcvx", 1.0}, {"bxcvy", 0.0}, {"bxcvz", 0.0}}});

  Options options{{"units", {{"seconds", 1.0}, {"Tesla", 1.0}, {"meters", 1.0}}},
                  {"test", {{"diamagnetic", true}}}};

  Vorticity component("test", options, &solver);
  component.declareAllSpecies({{"d+"}, {}});

  Options state{{"species", {{"d+", {{"AA", 2.0}, {"charge", 1.0}, {"pressure", 1.0}}}}}};

  component.transform(state);

  // Diamagnetic terms set both the divergence of the current,
  // and energy exchange with charged species
  ASSERT_TRUE(state["species"]["d+"].isSet("energy_source"));
  ASSERT_TRUE(state["fields"].isSet("DivJdia"));
}

TEST_F(VorticityTest, KinematicViscosity) {
  FakeSolver solver;

  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";
  Options options{{"units", {{"seconds", 1.0}, {"Tesla", 1.0}, {"meters", 1.0}}},
                  {"test", {{"viscosity", 1.0}, {"diamagnetic", false}}}};

  Vorticity component("test", options, &solver);
  component.declareAllSpecies({{"e"}, {"d"}, {"d+"}, {}});

  Options state{{"species",
                 {{"d+", {{"AA", 2.0}, {"charge", 1.0}, {"density", 1.0}}},
                  {"d", {{"AA", 2.0}, {"density", 1.0}}}}}};
  component.transform(state);

  ASSERT_TRUE(state["species"]["d+"].isSet("energy_source"));
  ASSERT_FALSE(state["species"]["d"].isSet("energy_source"));
}

TEST_F(VorticityTest, calculatePihat) {
  FakeSolver solver;
  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";

  Options options{{"units", {{"seconds", 1.0}, {"Tesla", 1.0}, {"meters", 1.0}}},
                  {"test", {{"average_atomic_mass", 2.5}}}};

  Vorticity component("test", options, &solver);

  Options state{{"species", {{"d+", {{"pressure", 1.2}, {"AA", 2.0}, {"charge", 1.4}}}}}};
  Permissions permissions{{readOnly("species"), readWrite("species:d+:energy_source")}};
  GuardedOptions guarded_state{&state, &permissions};

  Field3D Pi_hat = component.calculatePihat(guarded_state["species"]);

  ASSERT_TRUE(IsFieldEqual(Pi_hat, 1.2 * 2.0 / 2.5 / 1.4, "RGN_NOBNDRY"));
}

TEST_F(VorticityTest, calculateDivJdia) {
  FakeSolver solver;
  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";

  Options options{{"units", {{"seconds", 1.0}, {"Tesla", 1.0}, {"meters", 1.0}}},
                  {"test", {{"average_atomic_mass", 2.5}}}};

  Vorticity component("test", options, &solver);

  Options state{{"species", {{"d+", {{"pressure", 1.2}, {"AA", 2.0}, {"charge", 1.4}}}}}};
  Permissions permissions{{readOnly("species"), readWrite("species:d+:energy_source")}};
  GuardedOptions guarded_state{&state, &permissions};

  Field3D DivJdia = component.calculateDivJdia(0.0, guarded_state["species"]);

  ASSERT_TRUE(IsFieldEqual(DivJdia, 0.0, "RGN_NOBNDRY"));
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["d+"]["energy_source"]), 0.0,
                           "RGN_NOBNDRY"));
}

class VorticityBoundaryTest : public FakeMeshFixture_tmpl<10, 1, 8, false, 2> {
  void SetUp() override {
    Laplacian::cleanup();
    mesh->addBoundary(bout::boundary::NewBoundaryRegionXIn("core", 0, 0, mesh));
    mesh->addBoundary(bout::boundary::NewBoundaryRegionXOut("sol", 0, 0, mesh));
  }
  void TearDown() override { Laplacian::cleanup(); }
};

TEST_F(VorticityBoundaryTest, RelaxationIsPackedStateAndBacktrackingSafe) {
  Options::root()["mesh"]["paralleltransform"]["type"] = "identity";
  Options::root()["Vort"]["bndry_all"] = "neumann";
  Options options{{"units", {{"seconds", 0.2}, {"Tesla", 1.0}, {"meters", 1.0}}},
                  {"test",
                   {{"phi_boundary_relax", true},
                    {"phi_boundary_timescale", 0.4},
                    {"diamagnetic", false},
                    {"phi_dissipation", false},
                    {"viscosity", 0.1}}}};
  FakeSolver solver;
  Vorticity component("test", options, &solver);
  ASSERT_TRUE(solver.evolvesBoundary("Vort"));
  Field3D encoded{0.0};
  for (const auto& i : encoded.getRegion("RGN_ALL")) {
    if (i.x() < mesh->xstart) {
      encoded[i] = 0.08;
    }
    if (i.x() > mesh->xend) {
      encoded[i] = -0.03;
    }
  }
  solver.setState({{"Vort", encoded}});
  auto packed = solver.packState();
  for (BoutReal time : {5.0, 5.0, 1.0}) {
    solver.applyValueBoundaries(time);
    Options state{{"time", time}};
    component.transform(state);
    component.finally(state);
    solver.applyDerivativeBoundaries();
    EXPECT_EQ(solver.packState(), packed);
    const auto phi = state["fields"]["phi"].as<Field3D>();
    const auto physical = state["fields"]["vorticity"].as<Field3D>();
    const auto rate = solver.getTimeDerivs()["Vort"].as<Field3D>();
    for (const auto side : {std::pair{mesh->xstart, mesh->xstart - 1},
                            std::pair{mesh->xend, mesh->xend + 1}}) {
      const auto [cell, guard] = side;
      BoutReal target = 0.0;
      for (int z = 0; z < mesh->LocalNz; ++z) {
        target += phi(cell, 0, z) / mesh->LocalNz;
      }
      for (int z = 0; z < mesh->LocalNz; ++z) {
        EXPECT_NEAR((phi(cell, 0, z) + phi(guard, 0, z)) / 2, encoded(guard, 0, z),
                    1e-12);
        EXPECT_NEAR(rate(guard, 0, z), (target - encoded(guard, 0, z)) / 2, 1e-12);
        EXPECT_NEAR(physical(guard, 0, z), physical(cell, 0, z), 1e-12);
      }
    }
  }
  solver.unpackState(packed);
  EXPECT_EQ(solver.packState(), packed);
}

TEST_F(VorticityBoundaryTest, RestartMigrationIsOncePerInstance) {
  Options::root()["mesh"]["paralleltransform"]["type"] = "identity";
  Options options{{"units", {{"seconds", 1.0}, {"Tesla", 1.0}, {"meters", 1.0}}},
                  {"test", {{"phi_boundary_relax", true}}}};
  for (BoutReal initial : {0.2, 0.3}) {
    FakeSolver solver;
    Vorticity component("test", options, &solver);
    Options restart{{"phi", Field3D{initial}}};
    component.restartVars(restart);
    EXPECT_EQ(restart["vorticity_boundary_state_version"].as<int>(), 1);
    EXPECT_NEAR(solver.getState()["Vort"].as<Field3D>()(1, 0, 0), initial, 1e-15);
    solver.setState({{"Vort", Field3D{0.7}}});
    component.restartVars(restart);
    EXPECT_NEAR(solver.getState()["Vort"].as<Field3D>()(1, 0, 0), 0.7, 1e-15);
  }
}

TEST_F(VorticityBoundaryTest, RejectsInvalidBoundaryTimescale) {
  Options::root()["mesh"]["paralleltransform"]["type"] = "identity";
  Options options{
      {"units", {{"seconds", 1.0}, {"Tesla", 1.0}, {"meters", 1.0}}},
      {"test", {{"phi_boundary_relax", true}, {"phi_boundary_timescale", 0.0}}}};
  FakeSolver solver;
  EXPECT_THROW(Vorticity("test", options, &solver), BoutException);
}
