#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx"

#include "../../include/classical_diffusion.hxx"

#include <bout/constants.hxx>
#include <bout/field_factory.hxx>

using bout::globals::mesh;

class ClassicalDiffusionTest : public FakeMeshFixture_tmpl<7, 5, 8> {
protected:
  Options options(BoutReal meters = 1.0, BoutReal seconds = 1.0) {
    return Options{
        {"units",
         {{"meters", meters}, {"seconds", seconds}, {"inv_meters_cubed", 1e19}}}};
  }

  Options diagnostics(BoutReal rho_s0 = 1.0, BoutReal Omega_ci = 1.0) {
    return Options{
        {"rho_s0", rho_s0}, {"Omega_ci", Omega_ci}, {"Nnorm", 1e19}, {"Tnorm", 1.0}};
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

TEST_F(ClassicalDiffusionTest, CustomDiffusionDampsThermalEnergyAtConstantDensity) {
  auto config = options(2.0, 0.5);
  config["classical_diffusion"]["custom_D"] = 3.2; // 0.4 in normalised units
  config["classical_diffusion"]["diagnose"] = true;
  const auto pressure =
      makeField<Field3D>([](const Ind3D& i) { return 2.0 + cos(2.0 * PI * i.z() / nz); });
  const Field3D expected_energy = -1.5 * 0.4 * 4.0 * SQ(sin(PI / nz)) * (pressure - 2.0);

  // Pressure takes precedence when supplied. Temperature-only states recover
  // the same thermal-energy diffusion from P = N*T, without collisions.
  for (const bool supply_pressure : {false, true}) {
    SCOPED_TRACE(supply_pressure);
    ClassicalDiffusion component("classical_diffusion", config, nullptr);
    component.declareAllSpecies({"n"});
    Options state;
    auto& species = state["species"]["n"];
    species["density"] = Field3D{2.0};
    species["temperature"] = supply_pressure ? Field3D{99.0} : pressure / 2.0;
    if (supply_pressure) {
      species["pressure"] = pressure;
    }
    component.transform(state);
    auto output = diagnostics(2.0, 2.0);
    component.outputVars(output);
    EXPECT_DOUBLE_EQ(output["D_classical"].attributes["conversion"].as<BoutReal>(), 8.0);

    EXPECT_TRUE(
        IsFieldEqual(get<Field3D>(species["density_source"]), 0.0, "RGN_NOBNDRY"));
    EXPECT_TRUE(IsFieldEqual(get<Field3D>(species["energy_source"]), expected_energy,
                             "RGN_NOBNDRY", 1e-12));
  }
}

TEST_F(ClassicalDiffusionTest, PeriodicTransportIsConservativeAndTranslationInvariant) {
  mesh->periodicX = true;
  static_cast<FakeMesh*>(mesh)->copy_periodic_halos = true;
  const int cells = mesh->xend - mesh->xstart + 1;

  for (const BoutReal custom_D : {-1.0, 0.4}) {
    SCOPED_TRACE(custom_D);
    auto config = options();
    config["classical_diffusion"]["custom_D"] = custom_D;
    ClassicalDiffusion component("classical_diffusion", config, nullptr);
    component.declareAllSpecies({"e"});
    std::map<std::string, Field3D> reference;

    // Moving a profile through the periodic seam must give the same sources
    // as moving it through an interior face, including for variable Dn.
    for (const int shift : {0, 1}) {
      auto density = makeField<Field3D>([&](const Ind3D& i) {
        return 2.0 + cos(TWOPI * (i.x() - mesh->xstart + shift) / cells);
      });
      BOUT_FOR(i, density.getRegion("RGN_XGUARDS")) { density[i] = 99.0; }
      mesh->communicate(density);
      Options state;
      auto& electrons = state["species"]["e"];
      electrons["density"] = density;
      electrons["pressure"] = SQ(density);
      electrons["temperature"] = density;
      electrons["velocity"] = 0.5 * density;
      electrons["charge"] = -1.0;
      electrons["AA"] = 0.5;
      electrons["collision_frequency"] = Field3D{0.2};
      component.transform(state);

      for (const auto* source : {"density_source", "momentum_source", "energy_source"}) {
        SCOPED_TRACE(source);
        const auto rhs = get<Field3D>(electrons[source]);
        EXPECT_GT(max(abs(rhs), false, "RGN_NOBNDRY"), 0.0);
        EXPECT_NEAR(integral(rhs), 0.0, 1e-12);
        if (shift == 0) {
          reference[source] = rhs;
        } else {
          for (const auto& i : rhs.getRegion("RGN_NOBNDRY")) {
            const int x = i.x() == mesh->xend ? mesh->xstart : i.x() + 1;
            EXPECT_NEAR(rhs[i], reference[source](x, i.y(), i.z()), 1e-12);
          }
        }
      }
    }
  }
}
