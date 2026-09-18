#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/div_ops.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
class DivOpsTest : public FakeMeshFixture_tmpl<7, 5, 3, false, 2>,
                   public testing::WithParamInterface<std::tuple<bool, bool, bool>> {};

TEST_P(DivOpsTest, Div_n_bxGrad_f_B_XPPM) {
  Field3D n = 1.0;
  Field3D f = 1.0;
  auto all = GetParam();
  Div_n_bxGrad_f_B_XPPM(n, f, std::get<0>(all), std::get<1>(all), std::get<2>(all));
}
INSTANTIATE_TEST_SUITE_P(DIVOPS, DivOpsTest,
                         testing::Combine(testing::Bool(), testing::Bool(),
                                          testing::Bool()));
