
#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh

#include "../../include/amjuel_hyd_recombination.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

#include <bout/constants.hxx>
#include <field_factory.hxx> // For generating functions

// Reuse the "standard" fixture for FakeMesh
using HydrogenRCTest = FakeMeshFixture;

TEST_F(HydrogenRCTest, CreateComponent) {
  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}}};

  AmjuelHydRecombinationIsotope<'h'> component("test", options, nullptr);
}
