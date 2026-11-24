#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/div_ops.hxx"
#include "../../include/domain_balance.hxx"

#include <bout/fv_ops.hxx>
#include <bout/paralleltransform.hxx>

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

class ParallelBalanceTest : public FakeMeshFixture_tmpl<7, 6, 8, false, 2, 1> {};

TEST_F(ParallelBalanceTest, ShiftedFluxesIncludeBothOpenBoundaries) {
  static_cast<FakeMesh*>(mesh)->ix_separatrix = 0;
  Field2D shift = 0.0;
  for (const auto& i : shift.getRegion("RGN_ALL")) {
    shift[i] = 0.1 * i.y();
  }
  mesh->getCoordinates()->setParallelTransform(
      std::make_unique<ShiftedMetric>(*mesh, CELL_CENTRE, shift, 8.0));
  const Field3D density{2.0}, velocity{1.0}, speed{1.0};
  Field3D particle_flow, kinetic_flow;
  const Field3D divergence = FV::Div_par_mod(density, velocity, speed, particle_flow);
  FV::Div_par_fvv_heating(density, velocity, speed, kinetic_flow);
  DomainBalance ledger;
  ledger.addTargets("particles", "target", particle_flow);
  ledger.addTargets("kinetic", "target", kinetic_flow);
  const auto lower = "target_" + DomainBalance::target(mesh, true);
  const auto upper = "target_" + DomainBalance::target(mesh, false);
  for (const auto* quantity : {"particles", "kinetic"}) {
    const auto incoming = ledger.local.at({quantity, lower});
    const auto outgoing = ledger.local.at({quantity, upper});
    EXPECT_LT(incoming, 0.0);
    EXPECT_GT(outgoing, 0.0);
    EXPECT_NEAR(incoming + outgoing, 0.0, 1e-12);
  }
  EXPECT_NEAR(DomainBalance::integral(divergence),
              ledger.local.at({"particles", lower})
                  + ledger.local.at({"particles", upper}),
              1e-12);
}
