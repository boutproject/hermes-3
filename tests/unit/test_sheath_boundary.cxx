
#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh

#include "../../include/sheath_boundary.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

#include <bout/constants.hxx>
#include <field_factory.hxx>  // For generating functions

// Reuse the "standard" fixture for FakeMesh
using SheathBoundaryTest = FakeMeshFixture;

TEST_F(SheathBoundaryTest, CreateComponent) {
  Options options;
  
  SheathBoundary component("test", options, nullptr);
}

TEST_F(SheathBoundaryTest, CalculatePotential) {
  Options options;
  
  SheathBoundary component("test", options, nullptr);

  Field3D N = FieldFactory::get()->create3D("1 + y", &options, mesh);
  BoutReal Te = 2.0;
  BoutReal Ti = 3.0;
  BoutReal Zi = 1.1;
  BoutReal si = 0.5;
  
  Options state {{"species",
                  {// Electrons
                   {"e", {{"density", N},
                          {"temperature", Te},
                          {"velocity", 0.0}}},
                   // Ion species
                   {"h", {{"density", si*N},
                          {"temperature", Ti},
                          {"AA", 1.0},
                          {"charge", Zi},
                          {"velocity", 0.0}}}}}};

  component.transform(state);
  
  // Should have calculated potential
  ASSERT_TRUE(state["fields"].isSet("phi"));

  // Calculate the expected value of phi
  const BoutReal adiabatic = 5./3;
  BoutReal Vzi = sqrt(adiabatic * Ti + Zi * Te);
  BoutReal phi_ref = Te * log(sqrt(Te * SI::Mp / SI::Me / TWOPI) / (si * Zi * Vzi));

  output.write("TEST: {:e} {:e} {:e}\n", Te, si * Zi * Vzi, phi_ref);
  output.write("ION: {:e} {:e} {:e}\n", adiabatic * Ti, Zi * Te * si / (si + 1), Vzi);
  
  Field3D phi = state["fields"]["phi"];
  
  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
    for (int jz = 0; jz < mesh->LocalNz; jz++) {
      ASSERT_DOUBLE_EQ(phi_ref, phi(r.ind, mesh->yend, jz));
    }
  }
}
