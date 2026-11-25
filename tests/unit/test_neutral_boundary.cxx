
#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/neutral_boundary.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

#include <bout/constants.hxx>

#include <cmath>
#include <vector>

// Reuse the "standard" fixture for FakeMesh
using NeutralBoundaryTest = FakeMeshFixture;

namespace {
/// Apply neutral_boundary to a neutral species "d" with uniform
/// (normalised) density and temperature, and return the energy source
Field3D neutralBoundaryEnergySource(Options& options, BoutReal density,
                                    BoutReal temperature) {
  NeutralBoundary component("d", options, nullptr);

  Options state;
  state["species"]["d"]["AA"] = 2.0;
  state["species"]["d"]["density"] = density;
  state["species"]["d"]["temperature"] = temperature;
  state["species"]["d"]["pressure"] = density * temperature;

  component.transform(state);

  return get<Field3D>(state["species"]["d"]["energy_source"]);
}

enum class Wall { target, sol, pfr };

/// Cells next to the given wall
std::vector<Ind3D> wallCells(const Field3D& f, Wall wall) {
  std::vector<Ind3D> cells;
  if (wall == Wall::target) {
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        cells.push_back(indexAt(f, r.ind, mesh->ystart, jz));
      }
    }
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        cells.push_back(indexAt(f, r.ind, mesh->yend, jz));
      }
    }
  } else {
    const int ix = (wall == Wall::sol) ? mesh->xend : mesh->xstart;
    for (int iy = 0; iy < mesh->LocalNy; iy++) {
      for (int iz = 0; iz < mesh->LocalNz; iz++) {
        cells.push_back(indexAt(f, ix, iy, iz));
      }
    }
  }
  return cells;
}

/// With only thermal reflection (fast reflection fraction 0), incident
/// neutrals carry 2 T per particle to the wall and thermally reflected
/// neutrals carry T_FC back. The energy source therefore vanishes for
/// T_FC = 2 T, is a sink below and a source above. This checks that the
/// value of thermal_refl_temperature is used at the given wall.
void checkThermalReflTemperature(Wall wall) {
  constexpr BoutReal Tnorm = 2.0; // eV
  constexpr BoutReal Nn = 1.0;
  constexpr BoutReal Tn = 1.5; // Normalised, i.e. 3 eV, so T_FC = 2 T is 6 eV

  for (const BoutReal T_FC : {5.0, 6.0, 7.0}) {
    Options options;
    options["units"]["eV"] = Tnorm;
    auto& d = options["d"];
    switch (wall) {
    case Wall::target:
      d["target_fast_refl_fraction"] = 0.0;
      break;
    case Wall::sol:
      d["neutral_boundary_lower_y"] = false;
      d["neutral_boundary_upper_y"] = false;
      d["neutral_boundary_sol"] = true;
      d["sol_fast_refl_fraction"] = 0.0;
      break;
    case Wall::pfr:
      d["neutral_boundary_lower_y"] = false;
      d["neutral_boundary_upper_y"] = false;
      d["neutral_boundary_pfr"] = true;
      d["pfr_fast_refl_fraction"] = 0.0;
      break;
    }
    d["thermal_refl_temperature"] = T_FC;

    const Field3D source = neutralBoundaryEnergySource(options, Nn, Tn);

    const auto cells = wallCells(source, wall);
    ASSERT_FALSE(cells.empty());
    for (const auto& i : cells) {
      if (T_FC < 6.0) {
        EXPECT_LT(source[i], 0.0) << "T_FC = " << T_FC;
      } else if (T_FC > 6.0) {
        EXPECT_GT(source[i], 0.0) << "T_FC = " << T_FC;
      } else {
        EXPECT_NEAR(source[i], 0.0, 1e-12) << "T_FC = " << T_FC;
      }
    }
  }
}
} // namespace

// Without the option, thermally reflected neutrals have the Franck-Condon
// temperature of 3 eV that was previously hardcoded
TEST_F(NeutralBoundaryTest, DefaultThermalReflTemperature) {
  constexpr BoutReal Tnorm = 5.0; // eV
  constexpr BoutReal Nn = 1.2;
  constexpr BoutReal Tn = 0.4; // Normalised
  constexpr BoutReal AA = 2.0;

  Options options_default;
  options_default["units"]["eV"] = Tnorm;
  const Field3D source_default = neutralBoundaryEnergySource(options_default, Nn, Tn);

  Options options_set;
  options_set["units"]["eV"] = Tnorm;
  options_set["d"]["thermal_refl_temperature"] = 3.0;
  const Field3D source_set = neutralBoundaryEnergySource(options_set, Nn, Tn);

  // Heat flux to the targets with the default reflection parameters
  // (energy reflection factor 0.75, fast reflection fraction 0.8) and
  // T_FC = 3 eV. The fake mesh has a unit metric, so the energy source
  // in the cell next to the target is minus this heat flux.
  const BoutReal v_th = 0.25 * std::sqrt(8 * Tn / (PI * AA));
  const BoutReal T_FC = 3.0 / Tnorm;
  const BoutReal q = 2 * Nn * Tn * v_th - (0.75 * 0.8) * 2 * Nn * Tn * v_th
                     - (1 - 0.8) * T_FC * Nn * v_th;

  const auto cells = wallCells(source_default, Wall::target);
  ASSERT_FALSE(cells.empty());
  for (const auto& i : cells) {
    EXPECT_NEAR(source_default[i], -q, 1e-12 * std::abs(q));
  }

  // Setting the default value explicitly gives the same result everywhere
  for (const auto& i : source_default.getRegion("RGN_ALL")) {
    EXPECT_EQ(source_default[i], source_set[i]);
  }
}

TEST_F(NeutralBoundaryTest, ThermalReflTemperatureTarget) {
  checkThermalReflTemperature(Wall::target);
}

TEST_F(NeutralBoundaryTest, ThermalReflTemperatureSOL) {
  checkThermalReflTemperature(Wall::sol);
}

TEST_F(NeutralBoundaryTest, ThermalReflTemperaturePFR) {
  // The PFR boundary is only applied where Y is not periodic
  static_cast<FakeMesh*>(mesh)->ix_separatrix = 0;
  checkThermalReflTemperature(Wall::pfr);
}
