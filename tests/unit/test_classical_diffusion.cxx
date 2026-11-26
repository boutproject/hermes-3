#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx"

#include "../../include/classical_diffusion.hxx"

#include <bout/field_factory.hxx>

using bout::globals::mesh;

class ClassicalDiffusionTest : public FakeMeshFixture_tmpl<7, 5, 8> {
protected:
  Options options() {
    return Options{
        {"units", {{"meters", 1.0}, {"seconds", 1.0}, {"inv_meters_cubed", 1e19}}}};
  }

  BoutReal integral(const Field3D& field) {
    const auto volume = mesh->getCoordinates()->cell_volume();
    BoutReal result = 0.0;
    for (const auto& i : field.getRegion("RGN_NOBNDRY")) {
      result += field[i] * volume[i];
    }
    return result;
  }
};

TEST_F(ClassicalDiffusionTest, ClosedWallsConserveParticlesMomentumAndEnergy) {
  // Nonorthogonal X-Y geometry and nonzero wall gradients exercise closure of
  // the complete face flux, including mixed-metric contributions.
  static_cast<FakeMesh*>(mesh)->ix_separatrix = 0;
  CovariantMetricTensor covariant(1.0, 1.0, 1.0, 0.2, 0.0, 0.0, mesh);
  mesh->getCoordinates()->setMetricTensor(covariant.inverse(), covariant);

  for (const BoutReal custom_D : {-1.0, 0.4}) {
    SCOPED_TRACE(custom_D);
    auto config = options();
    config["classical_diffusion"]["custom_D"] = custom_D;
    ClassicalDiffusion component("classical_diffusion", config, nullptr);
    component.declareAllSpecies({"e"});
    Options state;
    auto& electrons = state["species"]["e"];
    const auto density = FieldFactory::get()->create3D(
        "2 + 0.03*x^2 + 0.02*y^2 + 0.01*x*y + 0.04*cos(z)", &config, mesh);
    const auto temperature =
        FieldFactory::get()->create3D("2 + 0.05*x + 0.1*y + 0.03*sin(z)", &config, mesh);
    electrons["density"] = density;
    electrons["temperature"] = temperature;
    electrons["pressure"] = density * temperature;
    electrons["velocity"] =
        FieldFactory::get()->create3D("0.4 + 0.02*x + 0.01*y", &config, mesh);
    electrons["charge"] = -1.0;
    electrons["AA"] = 0.5;
    electrons["collision_frequency"] = Field3D{0.2};
    component.transform(state);

    for (const auto* source : {"density_source", "momentum_source", "energy_source"}) {
      SCOPED_TRACE(source);
      const auto rhs = get<Field3D>(electrons[source]);
      EXPECT_GT(max(abs(rhs), false, "RGN_NOBNDRY"), 0.0);
      EXPECT_NEAR(integral(rhs), 0.0, 1e-12);
    }
  }
}
