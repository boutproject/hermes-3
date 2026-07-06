#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx"

#include "../../include/sheath_boundary_penalty.hxx"

#include <set>

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

using namespace bout::globals;

using SheathBoundaryPenaltyTest = FakeMeshFixture;

namespace {

std::set<int> regionIndices(const Region<Ind3D>& region) {
  std::set<int> indices;
  BOUT_FOR_SERIAL(i, region) { indices.insert(i.ind); }
  return indices;
}

} // namespace

TEST_F(SheathBoundaryPenaltyTest, BuildPenaltyRegionUsesStrictThreshold) {
  auto mask = makeField<Field3D>([](const Ind3D&) { return 0.0; });

  const auto below = indexAt(mask, mesh->xstart, mesh->ystart, 0);
  const auto equal = indexAt(mask, mesh->xstart, mesh->ystart + 1, 0);
  const auto above = indexAt(mask, mesh->xstart, mesh->ystart, 1);
  const auto large = indexAt(mask, mesh->xend, mesh->yend, 2);

  mask[below] = 0.9e-5;
  mask[equal] = 1.0e-5;
  mask[above] = 1.1e-5;
  mask[large] = 1.0;

  const auto region = SheathBoundaryPenalty::buildPenaltyRegion(mask);
  const auto indices = regionIndices(region);

  EXPECT_EQ(indices.size(), 2);
  EXPECT_FALSE(indices.count(below.ind));
  EXPECT_FALSE(indices.count(equal.ind));
  EXPECT_TRUE(indices.count(above.ind));
  EXPECT_TRUE(indices.count(large.ind));
}

TEST_F(SheathBoundaryPenaltyTest, PreparePenaltyMaskReturnsThresholdedRegion) {
  auto mask = makeField<Field3D>([](const Ind3D&) { return 0.0; });

  const auto active1 = indexAt(mask, mesh->xstart, mesh->ystart, 0);
  const auto active2 = indexAt(mask, mesh->xend, mesh->yend, mesh->LocalNz - 1);

  mask[active1] = 0.2;
  mask[active2] = 0.7;

  const auto prepared = SheathBoundaryPenalty::preparePenaltyMask(mask);
  const auto indices = regionIndices(prepared.region);

  EXPECT_EQ(indices.size(), 2);
  EXPECT_TRUE(indices.count(active1.ind));
  EXPECT_TRUE(indices.count(active2.ind));
  EXPECT_DOUBLE_EQ(prepared.mask[active1], 0.2);
  EXPECT_DOUBLE_EQ(prepared.mask[active2], 0.7);
}

TEST_F(SheathBoundaryPenaltyTest, PrepareFieldAlignedPenaltyMaskCopiesYGuards) {
  auto mask_fa = makeField<Field3D>([](const Ind3D&) { return 0.0; });

  for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
    for (int jz = 0; jz < mesh->LocalNz; ++jz) {
      const auto interior = indexAt(mask_fa, r.ind, mesh->ystart, jz);
      const auto guard = interior.ym();

      mask_fa[interior] = 10.0 + r.ind + 0.1 * jz;
      mask_fa[guard] = -1.0;
    }
  }

  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
    for (int jz = 0; jz < mesh->LocalNz; ++jz) {
      const auto interior = indexAt(mask_fa, r.ind, mesh->yend, jz);
      const auto guard = interior.yp();

      mask_fa[interior] = 20.0 + r.ind + 0.1 * jz;
      mask_fa[guard] = -2.0;
    }
  }

  const auto prepared =
      SheathBoundaryPenalty::prepareFieldAlignedPenaltyMask(mask_fa, *mesh);
  const auto indices = regionIndices(prepared.region);

  for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
    for (int jz = 0; jz < mesh->LocalNz; ++jz) {
      const auto interior = indexAt(prepared.mask, r.ind, mesh->ystart, jz);
      EXPECT_DOUBLE_EQ(prepared.mask[interior.ym()], prepared.mask[interior]);
      EXPECT_TRUE(indices.count(interior.ind));
    }
  }

  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
    for (int jz = 0; jz < mesh->LocalNz; ++jz) {
      const auto interior = indexAt(prepared.mask, r.ind, mesh->yend, jz);
      EXPECT_DOUBLE_EQ(prepared.mask[interior.yp()], prepared.mask[interior]);
      EXPECT_TRUE(indices.count(interior.ind));
    }
  }
}
