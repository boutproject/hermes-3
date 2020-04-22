
#include "gtest/gtest.h"

#include "../../include/integrate.hxx"

#include "test_extras.hxx" // FakeMesh

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using CellAverageTest = FakeMeshFixture;

TEST_F(CellAverageTest, ConstantValue) {
  Field3D field{1.0};
  Field3D result = cellAverage([](BoutReal) { return 3.0; },
                                field.getRegion("RGN_NOBNDRY"))(field);

  ASSERT_TRUE(result.isAllocated());
  ASSERT_TRUE(areFieldsCompatible(field, result));
  ASSERT_TRUE(IsFieldEqual(result, 3.0, "RGN_NOBNDRY"));
}

TEST_F(CellAverageTest, ConstantField) {
  Field3D field{1.0};
  Field3D result = cellAverage([](BoutReal val) { return val * 2.0; },
                                field.getRegion("RGN_NOBNDRY"))(field);

  ASSERT_TRUE(result.isAllocated());
  ASSERT_TRUE(areFieldsCompatible(field, result));
  ASSERT_TRUE(IsFieldEqual(result, 2.0, "RGN_NOBNDRY"));
}
