#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include "test_extras.hxx"

#include "../../include/evolve_density.hxx"
#include "../../include/evolve_energy.hxx"
#include "../../include/evolve_momentum.hxx"
#include "../../include/evolve_pressure.hxx"

#include <bout/field_factory.hxx>

#include <cmath>
#include <type_traits>

template <class Component>
class EvolveHyperZTest : public FakeMeshFixture_tmpl<7, 7, 8, false, 2, 2> {
protected:
  void SetUp() override {
    FakeMeshFixture_tmpl::SetUp();
    this->getMesh().getCoordinates()->setBxy(1.0);
    this->getMesh().getCoordinates()->setDz(2.0);
  }

  Options options() {
    return Options{
        {"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
        {"i", {{"AA", 1.0}, {"charge", 0.0}, {"numerical_viscous_heating", false}}}};
  }

  std::string variable() {
    if constexpr (std::is_same_v<Component, EvolveDensity>) {
      return "Ni";
    } else if constexpr (std::is_same_v<Component, EvolveMomentum>) {
      return "NVi";
    } else if constexpr (std::is_same_v<Component, EvolvePressure>) {
      return "Pi";
    } else {
      return "Ei";
    }
  }

  Field3D derivative(Options& opts) {
    FakeSolver solver;
    Component component("i", opts, &solver);
    const Field3D profile = makeField<Field3D>(
        [](Ind3D& i) { return 2.0 + std::cos(2.0 * PI * i.z() / 8.0); });
    solver.setState(Options{{variable(), profile}});
    Options state{
        {"fastest_wave", 0.0},
        {"species",
         {{"i", {{"AA", 1.0}, {"density", Field3D{2.0}}, {"velocity", Field3D{0.0}}}}}}};
    component.transform(state);
    // Isolate Z diffusion from parallel advection, including for momentum.
    state["species"]["i"]["velocity"] = Field3D{0.0};
    component.finally(state);
    return solver.getTimeDerivs()[variable()].template as<Field3D>();
  }

  void expectDamping(const Field3D& rhs, const Field2D& coefficient) {
    // The discrete fourth derivative has this eigenvalue for the first
    // periodic Fourier mode. The constant part of the profile is unchanged.
    const BoutReal eigenvalue = 16.0 * std::pow(std::sin(PI / 8.0), 4);
    for (const auto& i : rhs.getRegion("RGN_NOBNDRY")) {
      EXPECT_NEAR(rhs[i],
                  -coefficient(i.x(), i.y()) * eigenvalue
                      * std::cos(2.0 * PI * i.z() / 8.0),
                  1e-12)
          << "at " << i.x() << ", " << i.y() << ", " << i.z();
    }
  }
};

using HyperZComponents =
    testing::Types<EvolveDensity, EvolveMomentum, EvolvePressure, EvolveEnergy>;
TYPED_TEST_SUITE(EvolveHyperZTest, HyperZComponents);

TYPED_TEST(EvolveHyperZTest, SpatialExpressionDampsEachXYPlane) {
  auto opts = this->options();
  opts["i"]["hyper_z"] = "0.01 + 0.001*x + 0.002*y";
  const Field2D coefficient =
      FieldFactory::get()->create2D("0.01 + 0.001*x + 0.002*y", &opts, &this->getMesh());
  this->expectDamping(this->derivative(opts), coefficient);
}

TYPED_TEST(EvolveHyperZTest, SuppliedFieldUsesNegativeAndZeroRegions) {
  auto opts = this->options();
  const Field2D coefficient =
      makeField<Field2D>([](Ind2D& i) { return 0.01 * (i.x() - 3) * (i.y() - 3); });
  opts["i"]["hyper_z"] = coefficient;
  this->expectDamping(this->derivative(opts), coefficient);
}

TYPED_TEST(EvolveHyperZTest, AbsentAndScalarCoefficients) {
  auto absent = this->options();
  this->expectDamping(this->derivative(absent), Field2D{0.0});
  for (const BoutReal coefficient : {0.0, 0.01, -0.01}) {
    auto opts = this->options();
    opts["i"]["hyper_z"] = coefficient;
    this->expectDamping(this->derivative(opts), Field2D{coefficient});
  }
}
