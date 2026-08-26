#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/component.hxx"
#include "../../include/relax_potential.hxx"

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

/// Check a RelaxPotential object can be successfully constructed and
/// registers vorticity and phi variables with the solver
TEST_F(RelaxPotentialTest, CreateComponent) {
  FakeSolver solver;

  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";
  Options options{{"units",
                   {{"seconds", 1.0},
                    {"Tesla", 1.0},
                    {"meters", 1.0},
                    {"eV", 100.},
                    {"inv_meters_cubed", 1e19}}}};

  RelaxPotential component("test", options, &solver);

  Options state = solver.getState();

  ASSERT_TRUE(state.isSet("Vort"));
  ASSERT_TRUE(state.isSet("phi1"));
}

/// Check that the RelaxPotential::transform method sets voricity and
/// phi when the component is instantiated with default parameters
TEST_F(RelaxPotentialTest, Transform) {
  FakeSolver solver;

  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";
  Options options{{"units",
                   {{"seconds", 1.0},
                    {"Tesla", 1.0},
                    {"meters", 1.0},
                    {"eV", 100.},
                    {"inv_meters_cubed", 1e19}}}};

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
  coords->Bxy = 1.0;

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
