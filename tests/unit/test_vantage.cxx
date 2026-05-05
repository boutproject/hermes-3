#include <cmath>

#include <bout/bout_types.hxx>
#include <bout/region.hxx>

#include "gtest/gtest.h"

#include "../../include/component.hxx"
#include "../../include/vantage.hxx"
#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Constructor for VantageTest which inherits from FakeMeshFixture.
// Construct a blank options and add fields to the mesh.
class VantageTest : public FakeMeshFixture {
public:
  VantageTest() : FakeMeshFixture() {}
  Options alloptions;

protected:
  // SetUp() is a gtest hook which runs something before
  // each test.
  void SetUp() override {
    auto* fake_mesh = static_cast<FakeMesh*>(mesh);

    Field2D Rxy{0.0, mesh};
    Field2D Zxy{0.0, mesh};
    Field2D Rxy_corners{0.0, mesh};
    Field2D Rxy_lower_right_corners{0.0, mesh};
    Field2D Rxy_upper_right_corners{0.0, mesh};
    Field2D Rxy_upper_left_corners{0.0, mesh};
    Field2D Zxy_corners{0.0, mesh};
    Field2D Zxy_lower_right_corners{0.0, mesh};
    Field2D Zxy_upper_right_corners{0.0, mesh};
    Field2D Zxy_upper_left_corners{0.0, mesh};

    // Create R and Z in FakeMesh
    for (int ix = mesh->xstart; ix <= mesh->xend; ++ix) {
      for (int iy = mesh->ystart; iy <= mesh->yend; ++iy) {

        // x and y are in index space
        BoutReal x = ix - mesh->xstart;
        BoutReal y = iy - mesh->ystart;
        BoutReal dx = 1.0;
        BoutReal dy = 1.0;
        Rxy(ix, iy) = x;
        Zxy(ix, iy) = y;
        Rxy_corners(ix, iy) = x - 0.5 * dx;
        Rxy_lower_right_corners(ix, iy) = x + 0.5 * dx;
        Rxy_upper_left_corners(ix, iy) = x - 0.5 * dx;
        Rxy_upper_right_corners(ix, iy) = x + 0.5 * dx;
        Zxy(ix, iy) = y;
        Zxy_corners(ix, iy) = y - 0.5 * dy;
        Zxy_lower_right_corners(ix, iy) = y - 0.5 * dy;
        Zxy_upper_left_corners(ix, iy) = y + 0.5 * dy;
        Zxy_upper_right_corners(ix, iy) = y + 0.5 * dy;
      }
    }

    // Populate FakeMesh GridDataSource with RZ data
    fake_mesh->setGridDataSource(
        new FakeGridDataSource{{{"Rxy", Rxy},
                                {"Zxy", Zxy},
                                {"Rxy_corners", Rxy_corners},
                                {"Rxy_lower_right_corners", Rxy_lower_right_corners},
                                {"Rxy_upper_right_corners", Rxy_upper_right_corners},
                                {"Rxy_upper_left_corners", Rxy_upper_left_corners},
                                {"Zxy_corners", Zxy_corners},
                                {"Zxy_lower_right_corners", Zxy_lower_right_corners},
                                {"Zxy_upper_right_corners", Zxy_upper_right_corners},
                                {"Zxy_upper_left_corners", Zxy_upper_left_corners}}});
  }
};

Options MakeOptions() {
  Options alloptions;

  alloptions["restart_files"]["path"] = ".";

  alloptions["units"]["eV"] = 1.0;
  alloptions["units"]["meters"] = 1.0;
  alloptions["units"]["seconds"] = 1.0;
  alloptions["units"]["inv_meters_cubed"] = 1e19;

  alloptions["dmplex"]["use_cxx_ivertex"] = true;

  alloptions["vantage"]["initial_neutral_density"] = 1;
  alloptions["vantage"]["npart_per_cell"] = 1;
  alloptions["vantage"]["nsteps"] = 0;
  alloptions["vantage"]["dt"] = 0.001;
  return alloptions;
}

// Create VANTAGE component
// This includes the entire kinetic neutral execution at the moment,
// including DMPlex creation, particle initialisation etc.
TEST_F(VantageTest, CreateComponent) {

  Options alloptions = MakeOptions();
  Vantage const component("vantage", alloptions, nullptr);
}