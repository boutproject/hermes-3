#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx"

#include "../../include/sheath_boundary_penalty.hxx"

#include <bout/constants.hxx>

#include <cmath>
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

TEST_F(SheathBoundaryPenaltyTest, CalculateVolumetricPenaltyReturnsElectronPenaltyTerms) {
  auto mask = makeField<Field3D>([](const Ind3D&) { return 0.0; });
  const auto active = indexAt(mask, mesh->xstart, mesh->ystart, 0);
  mask[active] = 0.25;
  const auto region = SheathBoundaryPenalty::buildPenaltyRegion(mask, 0.0);

  auto Ne = makeField<Field3D>([](const Ind3D&) { return 0.5; });
  auto Te = makeField<Field3D>([](const Ind3D&) { return 2.0; });
  auto Ve = makeField<Field3D>([](const Ind3D&) { return 3.0; });
  auto density_source = makeField<Field3D>([](const Ind3D&) { return 7.0; });
  auto momentum_source = makeField<Field3D>([](const Ind3D&) { return 11.0; });
  auto energy_source = makeField<Field3D>([](const Ind3D&) { return 13.0; });
  const SheathBoundaryPenalty::PenaltyMaskData penalty_data{mask, region};

  const BoutReal Me = 0.1;
  const BoutReal gamma_e = 3.5;
  const BoutReal penalty_timescale = 4.0;

  const auto result = SheathBoundaryPenalty::calculateVolumetricPenalty(
      penalty_data, Ne, Te, Ve, Me, gamma_e, penalty_timescale, density_source,
      momentum_source, energy_source);

  EXPECT_DOUBLE_EQ(result.density[active],
                   -0.25 * 7.0 - 0.25 * (0.5 - 1e-5) / penalty_timescale);
  EXPECT_DOUBLE_EQ(result.momentum[active],
                   -0.25 * 11.0 - 0.25 * Me * 0.5 * 3.0 / penalty_timescale);
  EXPECT_DOUBLE_EQ(result.energy[active],
                   -0.25 * 13.0 - 0.25 * gamma_e * 0.5 * 2.0 / penalty_timescale);

  const auto other = indexAt(mask, mesh->xend, mesh->yend, mesh->LocalNz - 1);
  EXPECT_DOUBLE_EQ(result.density[other], 0.0);
  EXPECT_DOUBLE_EQ(result.momentum[other], 0.0);
  EXPECT_DOUBLE_EQ(result.energy[other], 0.0);
}

TEST_F(SheathBoundaryPenaltyTest, CalculateVolumetricPenaltyReturnsIonPenaltyTerms) {
  auto mask = makeField<Field3D>([](const Ind3D&) { return 0.0; });
  const auto active = indexAt(mask, mesh->xstart, mesh->ystart, 0);
  mask[active] = 0.5;
  const auto region = SheathBoundaryPenalty::buildPenaltyRegion(mask, 0.0);

  auto Ni = makeField<Field3D>([](const Ind3D&) { return 0.5; });
  auto Ti = makeField<Field3D>([](const Ind3D&) { return 2.0; });
  auto Vi = makeField<Field3D>([](const Ind3D&) { return 3.0; });
  auto density_source = makeField<Field3D>([](const Ind3D&) { return 5.0; });
  auto momentum_source = makeField<Field3D>([](const Ind3D&) { return 7.0; });
  auto energy_source = makeField<Field3D>([](const Ind3D&) { return 11.0; });
  const SheathBoundaryPenalty::PenaltyMaskData penalty_data{mask, region};

  const BoutReal Mi = 2.0;
  const BoutReal gamma_i = 4.0;
  const BoutReal penalty_timescale = 10.0;
  const BoutReal density_floor = 1.0;

  const auto result = SheathBoundaryPenalty::calculateVolumetricPenalty(
      penalty_data, Ni, Ti, Vi, Mi, gamma_i, penalty_timescale, density_source,
      momentum_source, energy_source, density_floor);

  EXPECT_DOUBLE_EQ(result.density[active], -0.5 * 5.0);
  EXPECT_DOUBLE_EQ(result.momentum[active],
                   -0.5 * 7.0 - 0.5 * Mi * density_floor * 3.0 / penalty_timescale);
  EXPECT_DOUBLE_EQ(result.energy[active],
                   -0.5 * 11.0 - 0.5 * gamma_i * density_floor * 2.0 / penalty_timescale);

  const auto other = indexAt(mask, mesh->xend, mesh->yend, mesh->LocalNz - 1);
  EXPECT_DOUBLE_EQ(result.density[other], 0.0);
  EXPECT_DOUBLE_EQ(result.momentum[other], 0.0);
  EXPECT_DOUBLE_EQ(result.energy[other], 0.0);
}

TEST_F(SheathBoundaryPenaltyTest,
       CalculateElectronSurfaceMomentumPenaltyAddsExpectedTerm) {
  auto mask_fa = makeField<Field3D>([](const Ind3D&) { return 0.0; });
  const auto active = indexAt(mask_fa, mesh->xstart, mesh->ystart, 0);
  const auto iyp = active.yp();
  const auto iym = active.ym();
  mask_fa[active] = 1.0;
  mask_fa[iyp] = 0.0;
  mask_fa[iym] = 1.0;
  const auto region = SheathBoundaryPenalty::buildPenaltyRegion(mask_fa, 0.0);

  auto Ne_fa = makeField<Field3D>([](const Ind3D&) { return 0.7; });
  auto Te_fa = makeField<Field3D>([](const Ind3D&) { return 2.0; });
  auto Ve_fa = makeField<Field3D>([](const Ind3D&) { return 0.5; });
  auto phi_fa = makeField<Field3D>([](const Ind3D&) { return 0.3; });
  const SheathBoundaryPenalty::PenaltyMaskData penalty_fa_data{mask_fa, region};

  const BoutReal Me = 0.2;
  const BoutReal penalty_timescale = 4.0;
  const BoutReal tesheath = 2.0;
  const BoutReal vesheath = 0.5;
  const BoutReal phisheath = 0.3;
  const BoutReal Cse =
      sqrt(tesheath / (TWOPI * Me)) * exp(-phisheath / BOUTMAX(tesheath, 1e-5));
  const BoutReal expected = Me * 0.7 * (-Cse - vesheath) / penalty_timescale;

  const auto result = SheathBoundaryPenalty::calculateElectronSurfaceMomentumPenalty(
      penalty_fa_data, Ne_fa, Te_fa, Ve_fa, phi_fa, Me, penalty_timescale);

  EXPECT_DOUBLE_EQ(result[active], expected);
  const auto other = indexAt(mask_fa, mesh->xend, mesh->yend, mesh->LocalNz - 1);
  EXPECT_DOUBLE_EQ(result[other], 0.0);
}

TEST_F(SheathBoundaryPenaltyTest, CalculateIonSurfaceMomentumPenaltyAddsExpectedTerm) {
  auto mask_fa = makeField<Field3D>([](const Ind3D&) { return 0.0; });
  const auto active = indexAt(mask_fa, mesh->xstart, mesh->yend, 0);
  const auto iyp = active.yp();
  const auto iym = active.ym();
  mask_fa[active] = 1.0;
  mask_fa[iyp] = 1.0;
  mask_fa[iym] = 0.0;
  const auto region = SheathBoundaryPenalty::buildPenaltyRegion(mask_fa, 0.0);

  auto Ni_fa = makeField<Field3D>([](const Ind3D&) { return 0.4; });
  auto Ti_fa = makeField<Field3D>([](const Ind3D&) { return 1.5; });
  auto Te_fa = makeField<Field3D>([](const Ind3D&) { return 2.5; });
  auto Vi_fa = makeField<Field3D>([](const Ind3D&) { return -0.2; });
  const SheathBoundaryPenalty::PenaltyMaskData penalty_fa_data{mask_fa, region};

  const BoutReal Mi = 2.0;
  const BoutReal penalty_timescale = 5.0;
  const BoutReal density_floor = 0.6;
  const BoutReal Cs = sqrt((2.5 + 1.5) / Mi);
  const BoutReal expected = Mi * density_floor * (Cs - (-0.2)) / penalty_timescale;

  const auto result = SheathBoundaryPenalty::calculateIonSurfaceMomentumPenalty(
      penalty_fa_data, Ni_fa, Ti_fa, Te_fa, Vi_fa, Mi, penalty_timescale, density_floor);

  EXPECT_DOUBLE_EQ(result[active], expected);
  const auto other = indexAt(mask_fa, mesh->xstart, mesh->ystart, mesh->LocalNz - 1);
  EXPECT_DOUBLE_EQ(result[other], 0.0);
}
