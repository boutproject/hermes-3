
#include "gtest/gtest.h"

#include "../../include/integrate.hxx"

#include "test_extras.hxx" // FakeMesh
#include "fake_mesh_fixture.hxx"

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

// ─── cellAverageInto tests ───────────────────────────────────────────────────

using CellAverageIntoTest = FakeMeshFixture;

TEST_F(CellAverageIntoTest, ConstantValue) {
  Field3D field{1.0};
  Field3D result;
  cellAverageInto(result, [](BoutReal) { return 3.0; },
                  field.getRegion("RGN_NOBNDRY"), field);

  ASSERT_TRUE(result.isAllocated());
  ASSERT_TRUE(areFieldsCompatible(field, result));
  ASSERT_TRUE(IsFieldEqual(result, 3.0, "RGN_NOBNDRY"));
}

TEST_F(CellAverageIntoTest, ConstantField) {
  Field3D field{1.0};
  Field3D result;
  cellAverageInto(result, [](BoutReal val) { return val * 2.0; },
                  field.getRegion("RGN_NOBNDRY"), field);

  ASSERT_TRUE(result.isAllocated());
  ASSERT_TRUE(areFieldsCompatible(field, result));
  ASSERT_TRUE(IsFieldEqual(result, 2.0, "RGN_NOBNDRY"));
}

// Verify cellAverageInto gives bit-identical results to cellAverage.
TEST_F(CellAverageIntoTest, MatchesCellAverage) {
  Field3D Ne{1e19};
  Field3D Te{10.0};

  Field3D expected = cellAverage(
      [](BoutReal ne, BoutReal te) { return ne * te; },
      Ne.getRegion("RGN_NOBNDRY"))(Ne, Te);

  Field3D result;
  cellAverageInto(result,
      [](BoutReal ne, BoutReal te) { return ne * te; },
      Ne.getRegion("RGN_NOBNDRY"), Ne, Te);

  ASSERT_TRUE(IsFieldEqual(result, expected, "RGN_NOBNDRY", 0.0));
}

// Verify that calling cellAverageInto a second time reuses the allocation
// and updates the values correctly (workspace reuse path).
TEST_F(CellAverageIntoTest, WorkspaceReuse) {
  Field3D field_a{2.0};
  Field3D field_b{5.0};
  Field3D result;

  cellAverageInto(result, [](BoutReal v) { return v; },
                  field_a.getRegion("RGN_NOBNDRY"), field_a);
  ASSERT_TRUE(IsFieldEqual(result, 2.0, "RGN_NOBNDRY"));

  // Second call must overwrite with field_b's values, not keep field_a's.
  cellAverageInto(result, [](BoutReal v) { return v; },
                  field_b.getRegion("RGN_NOBNDRY"), field_b);
  ASSERT_TRUE(IsFieldEqual(result, 5.0, "RGN_NOBNDRY"));
}
