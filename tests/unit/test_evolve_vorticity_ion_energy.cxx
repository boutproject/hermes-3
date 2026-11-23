#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include "reynolds_reference.hxx"
#include "test_extras.hxx"

#include "../../include/component.hxx"
#include "../../include/component_scheduler.hxx"
#include "../../include/evolve_pressure.hxx"
#include "../../include/evolve_vorticity_ion_energy.hxx"
#include "../../include/polarisation_drift.hxx"
#include "../../include/sound_speed.hxx"

#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/constants.hxx>
#include <bout/coordinates.hxx>
#include <bout/derivs.hxx>
#include <bout/difops.hxx>
#include <bout/field_factory.hxx>
#include <bout/globals.hxx>
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <random>
#include <string>
#include <tuple>

using bout::globals::mesh;
using bout::ion_energy::PolarisationACBoundary;

// A test harness for the same two separately registered production components.
class IonEnergyTestSystem : public EvolveVorticityIonEnergy {
  static Options& prepare(const std::string& name, Options& options) {
    options[name]["type"] = std::string{"evolve_vorticity_ion_energy"};
    options["polarisation_drift"]["type"] = std::string{"polarisation_drift"};
    options["polarisation_drift"]["implementation"] = std::string{"conservative"};
    return options;
  }
  PolarisationDrift drift;

public:
  IonEnergyTestSystem(std::string name, Options& options, Solver* solver)
      : EvolveVorticityIonEnergy(name, prepare(name, options), solver),
        drift("polarisation_drift", options, solver) {}
  void declareAllSpecies(const SpeciesInformation& species) {
    EvolveVorticityIonEnergy::declareAllSpecies(species);
    drift.declareAllSpecies(species);
  }
  void transform(Options& state) {
    EvolveVorticityIonEnergy::transform(state);
    drift.transform(state);
  }
  void outputVars(Options& state) override {
    EvolveVorticityIonEnergy::outputVars(state);
    drift.outputVars(state);
  }
};

template <int NY, int MYG = 1, int NZ = 8>
class IonEnergyFixture : public FakeMeshFixture_tmpl<10, NY, NZ, false, 2, MYG> {
protected:
  // Laplacian solver state must not outlive the fixture's mesh.
  void SetUp() override {
    Laplacian::cleanup();
    mesh->addBoundary(
        bout::boundary::NewBoundaryRegionXIn("core", mesh->ystart, mesh->yend, mesh));
    mesh->addBoundary(
        bout::boundary::NewBoundaryRegionXOut("sol", mesh->ystart, mesh->yend, mesh));
  }
  void TearDown() override { Laplacian::cleanup(); }

  Options makeOptions() {
    Options::root()["mesh"]["paralleltransform"]["type"] = "identity";
    return {{"units",
             {{"eV", 1.0},
              {"inv_meters_cubed", 1.0},
              {"seconds", 1.0},
              {"Tesla", 1.0},
              {"meters", 1.0}}},
            {"polarisation_drift", {{"implementation", "conservative"}}},
            {"test",
             {{"average_atomic_mass", 1.5},
              {"phi_boundary_velocity", 0.4},
              {"diagnose", true},
              {"phi_rtol", 1e-10},
              {"phi_max_iterations", 100}}}};
  }

  Options makeState() {
    return {{"species",
             {{"e", {{"charge", -1.0}, {"AA", 1. / 1836}, {"pressure", Field3D{1.0}}}},
              {"d+", {{"charge", 1.0}, {"AA", 2.0}, {"density", Field3D{0.5}}}},
              {"he++", {{"charge", 2.0}, {"AA", 4.0}, {"density", Field3D{0.25}}}}}}};
  }

  void declareSpecies(IonEnergyTestSystem& component) {
    component.declareAllSpecies(SpeciesInformation({"e"}, {}, {"d+", "he++"}, {}));
  }

  Options diagnostics(IonEnergyTestSystem& component) {
    Options output{{"Nnorm", 1.0}, {"Tnorm", 1.0}, {"Omega_ci", 1.0}, {"rho_s0", 1.0}};
    component.outputVars(output);
    return output;
  }

  BoutReal integral(const Field3D& field) {
    const auto* coord = mesh->getCoordinates();
    BoutReal result = 0.0;
    for (const auto& i : field.getRegion("RGN_NOBNDRY")) {
      result +=
          field[i] * coord->J()[i] * coord->dx()[i] * coord->dy()[i] * coord->dz()[i];
    }
    return result;
  }
};

using EvolveVorticityIonEnergyTest = IonEnergyFixture<5>;
using IonEnergyParallelDampingTest = IonEnergyFixture<12, 2>;

using IonEnergySingletonYTest = IonEnergyFixture<1, 0>;

using IonEnergyDCWallTest = IonEnergyFixture<5, 1, 1>;

// Density boundary values need not be valid during thermodynamic recovery.
TEST_F(IonEnergyDCWallTest, ZeroDensityGuardsDoNotCreateHotWallEnergyFlux) {
  mesh->getCoordinates()->setNon_uniform(false);
  FakeSolver solver;
  auto options = makeOptions();
  options["test"]["radial_boundary"] = "no_slip";
  IonEnergyTestSystem component("test", options, &solver);
  declareSpecies(component);
  solver.setState({{"Vort", Field3D{0.0}}, {"Eion", Field3D{3.0}}});
  auto state = makeState();
  for (const auto* name : {"d+", "he++"}) {
    Field3D density = copy(state["species"][name]["density"].as<Field3D>());
    for (const auto& i : density.getRegion("RGN_NOY")) {
      if (i.x() < mesh->xstart || i.x() > mesh->xend) {
        density[i] = 0.0;
      }
    }
    state["species"][name]["density"] = density;
  }
  // A smooth parallel current drives a nonzero radial polarisation wall flux.
  state["species"]["e"]["momentum"] =
      makeField<Field3D>([](const Ind3D& i) { return 1e-5 * i.y() / 1836; });
  component.transform(state);
  component.finally(state);
  const auto output = diagnostics(component);
  const Field3D temperature = output["Tion"].as<Field3D>();
  const Field3D phi = output["phi"].as<Field3D>();
  for (const auto* name : {"d+", "he++"}) {
    const Field3D pressure = state["species"][name]["pressure"].as<Field3D>();
    const Field3D density = state["species"][name]["density"].as<Field3D>();
    for (const auto& i : pressure.getRegion("RGN_NOY")) {
      const int interior = std::clamp(i.x(), mesh->xstart, mesh->xend);
      EXPECT_NEAR(temperature[i], 8. / 3, 1e-12);
      EXPECT_NEAR(pressure[i], pressure(interior, i.y(), i.z()), 1e-12);
      EXPECT_NEAR(phi[i], 0.0, 1e-12);
      if (i.x() != interior) {
        EXPECT_EQ(density[i], 0.0); // Recovery must preserve the density input.
      }
    }
  }
  const Field3D flux = output["polarisation_flux_x"].as<Field3D>();
  EXPECT_GT(max(abs(flux), true), 1e-8);
  EXPECT_LT(max(abs(output["DivEpol"].as<Field3D>()), true), 1e-3);
  EXPECT_LT(max(abs(solver.getTimeDerivs()["Eion"].as<Field3D>()), true), 1e-3);
}

// Recovery and finite-volume energy identities.
TEST_F(EvolveVorticityIonEnergyTest, UniformStateHasZeroRHS) {
  FakeSolver solver;
  auto options = makeOptions();
  options["components"] = "e, d+, he++, test, polarisation_drift";
  options["test"]["type"] = "evolve_vorticity_ion_energy";
  options["test"]["radial_boundary"] = "fixed";
  options["e"]["type"] = "quasineutral, evolve_pressure, zero_current";
  options["e"]["charge"] = -1.0;
  options["e"]["AA"] = 1. / 1836;
  options["e"]["thermal_conduction"] = false;
  for (const auto& name : {"d+", "he++"}) {
    options[name]["type"] = "evolve_density, evolve_momentum";
    options[name]["charge"] = name == std::string{"d+"} ? 1.0 : 2.0;
    options[name]["AA"] = name == std::string{"d+"} ? 2.0 : 4.0;
  }
  auto scheduler = ComponentScheduler::create(options, options, &solver);
  solver.setState({{"Vort", Field3D{0.0}},
                   {"Eion", Field3D{3.0}},
                   {"Pe", Field3D{1.0}},
                   {"Nd+", Field3D{0.5}},
                   {"Nhe++", Field3D{0.25}},
                   {"NVd+", Field3D{0.0}},
                   {"NVhe++", Field3D{0.0}}});
  Options state{{"time", 0.0}, {"fastest_wave", Field3D{1.0}}};
  scheduler->transform(state);
  for (const auto* ion : {"d+", "he++"}) {
    EXPECT_TRUE(IsFieldEqual(state["species"][ion]["temperature"].as<Field3D>(), 8. / 3,
                             "RGN_NOBNDRY", 1e-12));
  }
  const auto derivatives = solver.getTimeDerivs();
  for (const auto& [name, value] : derivatives.getChildren()) {
    EXPECT_TRUE(IsFieldEqual(value.as<Field3D>(), 0.0, "RGN_NOBNDRY", 1e-13)) << name;
  }
  for (const auto& [name, species] : state["species"].getChildren()) {
    for (const auto& source : {"density_source", "momentum_source", "energy_source"}) {
      if (species.isSet(source)) {
        EXPECT_TRUE(
            IsFieldEqual(species[source].as<Field3D>(), 0.0, "RGN_NOBNDRY", 1e-13))
            << name << ":" << source;
      }
    }
  }
}

TEST_F(EvolveVorticityIonEnergyTest, PotentialRecoverySatisfiesFiniteVolumeVorticity) {
  auto* coord = mesh->getCoordinates();
  coord->setContravariantMetricTensor(
      ContravariantMetricTensor{1.0, 1.0, 1.0, 0.0, 0.3, 0.0});
  // Avoid the fixture's unrelated nonuniform first-derivative correction.
  coord->setNon_uniform(false);
  const BoutReal kz = 2 * PI / nz;
  const Field3D Phi = makeField<Field3D>([&](const Ind3D& i) {
    return 0.05 * sin((0.4 * (i.x() - mesh->xstart)) + (3 * kz * i.z()));
  });
  const Field3D coefficient{Coordinates::FieldMetric{1.5 / SQ(coord->Bxy())}};
  const Field3D vort = bout::ion_energy::Div_a_Grad_perp_XZ(coefficient, Phi);
  const Field3D energy =
      makeField<Field3D>([&](const Ind3D& i) { return 3.0 + 0.1 * cos(2 * kz * i.z()); });
  for (int iterations : {1, 100}) {
    SCOPED_TRACE(iterations);
    FakeSolver solver;
    auto options = makeOptions();
    options["test"]["phi_max_iterations"].force(iterations);
    IonEnergyTestSystem component("test", options, &solver);
    declareSpecies(component);
    solver.setState({{"Vort", vort}, {"Eion", energy}});
    auto state = makeState();
    if (iterations == 1) {
      EXPECT_THROW(component.transform(state), BoutException);
      continue;
    }
    component.transform(state);
    component.finally(state);
    auto output = diagnostics(component);
    const Field3D phi = output["phi"].as<Field3D>();
    const Field3D psi = output["ion_pressure_potential"].as<Field3D>();
    const Field3D forward = bout::ion_energy::Div_a_Grad_perp_XZ(coefficient, phi + psi);
    EXPECT_TRUE(IsFieldEqual(forward, vort, "RGN_NOBNDRY", 3e-10));
    const Field3D partition = 1.5
                                  * (get<Field3D>(state["species"]["d+"]["pressure"])
                                     + get<Field3D>(state["species"]["he++"]["pressure"]))
                              + output["ion_KE_perp"].as<Field3D>();
    EXPECT_TRUE(IsFieldEqual(partition, energy, "RGN_NOBNDRY", 1e-12));
  }
}

TEST_F(EvolveVorticityIonEnergyTest, IntegrationByPartsOnNonuniformCrossMetricMesh) {
  auto* coord = mesh->getCoordinates();
  const BoutReal kx = 2 * PI / (mesh->xend - mesh->xstart + 1);
  const BoutReal kz = 2 * PI / nz;
  // Profiles, coefficients and guard cells are periodic in X and Z. In a
  // metric-3D build this also exercises metrics varying in Z.
  const auto phase = [&](const auto& i) {
    const BoutReal value = kx * (i.x() - mesh->xstart);
    if constexpr (bout::build::use_metric_3d) {
      value += kz * i.z();
    }
    return value;
  };
  using Metric = Coordinates::FieldMetric;
  coord->setDx(
      makeField<Metric>([&](const auto& i) { return 1.0 + 0.2 * cos(phase(i)); }));
  coord->setDz(
      makeField<Metric>([&](const auto& i) { return 1.0 + 0.1 * sin(phase(i)); }));
  coord->setJ(
      makeField<Metric>([&](const auto& i) { return 1.3 + 0.1 * cos(phase(i)); }));
  coord->setContravariantMetricTensor(ContravariantMetricTensor{
      makeField<Metric>([&](const auto& i) { return 1.2 + 0.1 * sin(phase(i)); }),
      Metric{1.0}, Metric{1.4}, Metric{0.0},
      makeField<Metric>([&](const auto& i) { return 0.3 + 0.1 * cos(phase(i)); }),
      Metric{0.0}});
  const Field3D a = makeField<Field3D>([&](const Ind3D& i) {
    return 1.0 + 0.2 * sin(kx * (i.x() - mesh->xstart) - kz * i.z());
  });
  const Field3D f = makeField<Field3D>([&](const Ind3D& i) {
    return cos(kx * (i.x() - mesh->xstart) + kz * i.z()) + 0.3 * sin(2 * kz * i.z());
  });
  const Field3D h = makeField<Field3D>(
      [&](const Ind3D& i) { return sin(kx * (i.x() - mesh->xstart) + 2 * kz * i.z()); });
  const Field3D Lf = bout::ion_energy::Div_a_Grad_perp_XZ(a, f);
  const Field3D Lh = bout::ion_energy::Div_a_Grad_perp_XZ(a, h);
  const Field3D ke = bout::ion_energy::kineticEnergyPerp(a, f);
  BoutReal fLf = 0.0, hLf = 0.0, fLh = 0.0, energy = 0.0, divergence = 0.0;
  for (const auto& i : f.getRegion("RGN_NOBNDRY")) {
    const BoutReal volume = coord->J()[i] * coord->dx()[i] * coord->dz()[i];
    fLf += volume * f[i] * Lf[i];
    hLf += volume * h[i] * Lf[i];
    fLh += volume * f[i] * Lh[i];
    energy += volume * ke[i];
    divergence += volume * Lf[i];
  }
  EXPECT_GT(energy, 0.0);
  EXPECT_NEAR(fLf, -2 * energy, 1e-11);
  EXPECT_NEAR(hLf, fLh, 1e-11);
  EXPECT_NEAR(divergence, 0.0, 1e-11);
}

TEST_F(EvolveVorticityIonEnergyTest, IntegrationByPartsIncludesRadialBoundaryWork) {
  auto* coord = mesh->getCoordinates();
  coord->setContravariantMetricTensor(
      ContravariantMetricTensor{1.0, 1.0, 1.0, 0.0, 0.3, 0.0});
  const Field3D f = makeField<Field3D>([&](const Ind3D& i) {
    return (0.2 * i.x()) + sin((0.4 * i.x()) + (2 * PI * i.z() / nz));
  });
  const Field3D Lf = bout::ion_energy::Div_a_Grad_perp_XZ(Field3D{1.5}, f);
  const Field3D ke = bout::ion_energy::kineticEnergyPerp(Field3D{1.5}, f);
  BoutReal fLf = 0.0, energy = 0.0, boundary_work = 0.0;
  for (const auto& i : f.getRegion("RGN_NOBNDRY")) {
    fLf += f[i] * Lf[i];
    energy += ke[i];
    if (i.x() == mesh->xend || i.x() == mesh->xstart) {
      const auto l = i.x() == mesh->xend ? i : i.xm();
      const auto r = l.xp();
      // On a uniform mesh the normal flux at a radial face uses the
      // centred transverse derivative averaged between its two cells.
      const BoutReal transverse = 0.25 * (f[l.zp()] - f[l.zm()] + f[r.zp()] - f[r.zm()]);
      const BoutReal flux = 1.5 * (f[r] - f[l] + (0.3 * transverse));
      const BoutReal outward = i.x() == mesh->xend ? 1.0 : -1.0;
      boundary_work += outward * 0.5 * (f[l] + f[r]) * flux;
    }
  }
  EXPECT_GT(abs(boundary_work), 1e-3);
  EXPECT_NEAR(fLf, -2 * energy + boundary_work, 1e-11);
}

// Energy transfers and species transport.
TEST_F(EvolveVorticityIonEnergyTest, ParallelFluxCarriesEnthalpyAndKineticEnergy) {
  static_cast<FakeMesh*>(mesh)->ix_separatrix = -1; // Open ends for the linear profile.
  FakeSolver solver;
  auto options = makeOptions();
  options["test"]["average_atomic_mass"].force(2.0);
  IonEnergyTestSystem component("test", options, &solver);
  component.declareAllSpecies(SpeciesInformation({"e"}, {}, {"d+"}, {}));
  solver.setState({{"Vort", Field3D{0.0}}, {"Eion", Field3D{1.5}}});
  Options state{{"species",
                 {{"e", {{"charge", -1.0}, {"pressure", Field3D{1.0}}}},
                  {"d+",
                   {{"charge", 1.0},
                    {"AA", 2.0},
                    {"density", Field3D{1.0}},
                    {"velocity", Field3D{0.0}}}}}}};
  component.transform(state);

  // Set the completed state independently of the inversion. With v=1 and
  // n=p=1+0.1*y, d/dy[(5/2*p + K)*v] = 0.35, including kinetic transport.
  state["species"]["d+"]["velocity"].force(Field3D{1.0});
  const Field3D profile =
      makeField<Field3D>([](const Ind3D& i) { return 1.0 + (0.1 * i.y()); });
  state["species"]["d+"]["density"].force(profile);
  state["species"]["d+"]["pressure"].force(profile);
  state["fastest_wave"] = Field3D{1.0};
  component.finally(state);
  auto output = diagnostics(component);

  EXPECT_TRUE(
      IsFieldEqual(output["DivEparallel"].as<Field3D>(), 0.35, "RGN_NOBNDRY", 1e-12));
  EXPECT_TRUE(IsFieldEqual(solver.getTimeDerivs()["Eion"].as<Field3D>(), -0.35,
                           "RGN_NOBNDRY", 1e-12));
}

TEST_F(EvolveVorticityIonEnergyTest, ParticleSourcesCarryTransferredEnergy) {
  for (const bool injection : {false, true}) {
    SCOPED_TRACE(injection);
    FakeSolver solver;
    auto options = makeOptions();
    IonEnergyTestSystem component("test", options, &solver);
    declareSpecies(component);
    solver.setState({{"Vort", Field3D{0.0}}, {"Eion", Field3D{10.0}}});
    auto state = makeState();
    auto& ion = state["species"]["d+"];
    ion["velocity"] = Field3D{injection ? 1.0 : 3.0};
    ion["density_source"] = Field3D{injection ? 1.0 : -1.0};
    ion["momentum_source"] = Field3D{injection ? 6.0 : -6.0};
    ion["energy_source"] = Field3D{injection ? 4.0 : -4.0};
    component.transform(state);
    component.finally(state);
    // Removal carries thermal energy and the existing ion kinetic energy;
    // injection at velocity 3 includes mixing heat relative to velocity 1.
    EXPECT_TRUE(IsFieldEqual(diagnostics(component)["Eion_src"].as<Field3D>(),
                             injection ? 9.0 : -13.0, "RGN_NOBNDRY", 1e-12));
  }
}

TEST_F(EvolveVorticityIonEnergyTest,
       ElectronHeatingDistinguishesCompressionAndElectricalWork) {
  for (bool curvature : {false, true}) {
    SCOPED_TRACE(curvature);
    static_cast<FakeMesh*>(mesh)->ix_separatrix = -1;
    static_cast<FakeMesh*>(mesh)->setGridDataSource(new FakeGridDataSource{
        {{"bxcvx", 0.0}, {"bxcvy", 0.0}, {"bxcvz", curvature ? 0.5 : 0.0}}});
    FakeSolver solver;
    auto options = makeOptions();
    options["e"]["poloidal_flows"] = false;
    options["e"]["bndry_flux"] = false;
    EvolvePressure electrons("e", options, &solver);
    IonEnergyTestSystem component("test", options, &solver);
    declareSpecies(component);
    electrons.declareAllSpecies(SpeciesInformation({"e"}, {}, {"d+", "he++"}, {}));
    const Field3D vort = makeField<Field3D>([&](const Ind3D& i) {
      return curvature ? 0.03 * cos(TWOPI * i.z() / nz) : 0.02 * i.y();
    });
    solver.setState({{"Vort", vort}, {"Eion", Field3D{3.0}}, {"Pe", Field3D{1.0}}});
    auto state = makeState();
    state["species"]["e"]["density"] = Field3D{1.0};
    if (!curvature) {
      state["species"]["e"]["momentum"] =
          makeField<Field3D>([](const Ind3D& i) { return -(1.0 + 0.1 * i.y()) / 1836; });
    }
    electrons.transform(state);
    component.transform(state);
    component.finally(state);
    electrons.finally(state);
    const auto output = diagnostics(component);
    // Curvature supplies thermal compression work; parallel electrical work
    // changes total energy without heating uniform electron pressure.
    const Field3D expected =
        curvature ? Field3D{-DDZ(output["phi"].as<Field3D>())} : Field3D{0.0};
    EXPECT_TRUE(IsFieldEqual(state["species"]["e"]["energy_source"].as<Field3D>(),
                             expected, "RGN_NOBNDRY", 1e-12));
    if (curvature) {
      EXPECT_GT(max(abs(expected), true), 1e-5);
    } else {
      EXPECT_GT(max(abs(output["phi_Div_Je"].as<Field3D>()), true), 1e-5);
      EXPECT_TRUE(IsFieldEqual(solver.getTimeDerivs()["Pe"].as<Field3D>(), 0.0,
                               "RGN_NOBNDRY", 1e-12));
    }
  }
}

TEST_F(EvolveVorticityIonEnergyTest, VaryingSpeciesFractionsPreserveChargeCurrent) {
  for (bool poloidal : {false, true}) {
    SCOPED_TRACE(poloidal);
    mesh->getCoordinates()->setContravariantMetricTensor(
        ContravariantMetricTensor{1.0, 1.0, 1.0, 0.0, 0.0, 0.3});
    FakeSolver solver;
    auto options = makeOptions();
    options["test"]["poloidal_flows"] = poloidal;
    IonEnergyTestSystem component("test", options, &solver);
    declareSpecies(component);
    solver.setState({{"Vort", makeField<Field3D>([](const Ind3D& i) {
                        return 0.02 * sin(0.8 * i.x() + 1.7 * i.y() + TWOPI * i.z() / nz);
                      })},
                     {"Eion", Field3D{4.0}}});
    auto state = makeState();
    state["species"]["d+"]["density"].force(makeField<Field3D>([](const Ind3D& i) {
      return 0.5 + (0.2 * sin((0.6 * i.x()) + (TWOPI * i.z() / nz)));
    }));
    state["species"]["he++"]["density"].force(makeField<Field3D>([](const Ind3D& i) {
      return 0.25 + 0.1 * cos(0.4 * i.x() - 2 * TWOPI * i.z() / nz);
    }));
    state["species"]["he++"]["AA"].force(3.0);
    for (const auto* name : {"d+", "he++"}) {
      auto& species = state["species"][name];
      const BoutReal velocity = std::string{name} == "d+" ? 2.0 : -1.0;
      species["velocity"] = Field3D{velocity};
      species["momentum"] = Field3D{species["AA"].as<BoutReal>() * velocity
                                    * species["density"].as<Field3D>()};
    }
    component.transform(state);
    component.finally(state);
    const Field3D current = diagnostics(component)["DivJpol"].as<Field3D>();
    const Field3D charge = state["species"]["d+"]["density_source"].as<Field3D>()
                           + 2 * state["species"]["he++"]["density_source"].as<Field3D>();
    EXPECT_GT(max(abs(current), true), 1e-10);
    EXPECT_TRUE(IsFieldEqual(diagnostics(component)["Eion_src"].as<Field3D>(), 0.0,
                             "RGN_NOBNDRY", 1e-12));
    for (const auto& i : current.getRegion("RGN_NOBNDRY")) {
      // FakeMesh does not communicate the species weights into Y halos.
      if (!poloidal || (i.y() > mesh->ystart && i.y() < mesh->yend)) {
        EXPECT_NEAR(charge[i], -current[i], 1e-12);
        for (const auto* name : {"d+", "he++"}) {
          const auto& species = state["species"][name];
          const BoutReal velocity = std::string{name} == "d+" ? 2.0 : -1.0;
          EXPECT_NEAR(species["momentum_source"].as<Field3D>()[i],
                      species["AA"].as<BoutReal>() * velocity
                          * species["density_source"].as<Field3D>()[i],
                      1e-12);
        }
      }
    }
  }
}

// Dissipation changes the thermal/kinetic partition and preserves energy.
TEST_F(EvolveVorticityIonEnergyTest, KinematicViscosityDampsZModeInPhysicalUnits) {
  const Field3D vort = makeField<Field3D>(
      [&](const Ind3D& i) { return 0.001 * cos(2 * PI * i.z() / nz); });
  const auto run = [&](bool viscous) {
    FakeSolver solver;
    auto options = makeOptions();
    options["units"]["meters"].force(2.0);
    options["units"]["seconds"].force(0.5);
    if (viscous) {
      options["test"]["viscosity"] = 8.0; // 8 / (2^2 * 2) = 1 normalised.
    }
    IonEnergyTestSystem component("test", options, &solver);
    declareSpecies(component);
    solver.setState({{"Vort", vort}, {"Eion", Field3D{10.0}}});
    auto state = makeState();
    component.transform(state);
    component.finally(state);
    auto output = diagnostics(component);
    output["vorticity_rhs"] = solver.getTimeDerivs()["Vort"].as<Field3D>();
    return output;
  };
  const auto baseline = run(false);
  const auto viscous = run(true);
  EXPECT_FALSE(baseline.isSet("kinematic_viscosity"));
  const Field3D expected = -4 * SQ(sin(PI / nz)) * vort;
  EXPECT_TRUE(IsFieldEqual(viscous["kinematic_viscosity"].as<Coordinates::FieldMetric>(),
                           1.0, "RGN_NOBNDRY", 1e-12));
  EXPECT_TRUE(
      IsFieldEqual(viscous["DivJviscous"].as<Field3D>(), expected, "RGN_NOBNDRY", 1e-12));
  EXPECT_TRUE(IsFieldEqual(Field3D{viscous["vorticity_rhs"].as<Field3D>()
                                   - baseline["vorticity_rhs"].as<Field3D>()},
                           expected, "RGN_NOBNDRY", 1e-12));
  for (const auto* field : {"ion_chi", "DivJpol", "DivEpol", "ddt(Eion)"}) {
    EXPECT_TRUE(IsFieldEqual(viscous[field].as<Field3D>(), baseline[field].as<Field3D>(),
                             "RGN_NOBNDRY", 1e-13))
        << field;
  }
  // Missing electron density is supported, as in the legacy component.
  EXPECT_TRUE(IsFieldEqual(viscous["electron_viscous_heating"].as<Field3D>(), 0.0,
                           "RGN_NOBNDRY", 0.0));
}

TEST_F(EvolveVorticityIonEnergyTest, ViscousHeatingPairsIonAndElectronEnergy) {
  mesh->getCoordinates()->setNon_uniform(false);
  FakeSolver solver;
  auto options = makeOptions();
  options["test"]["viscosity"] = 0.1;
  options["e"]["poloidal_flows"] = false;
  options["e"]["bndry_flux"] = false;
  EvolvePressure electrons("e", options, &solver);
  IonEnergyTestSystem component("test", options, &solver);
  declareSpecies(component);
  electrons.declareAllSpecies(SpeciesInformation({"e"}, {}, {"d+", "he++"}, {}));
  const Field3D vort = makeField<Field3D>([&](const Ind3D& i) {
    return 0.02
           * sin(2 * PI * (i.x() - mesh->xstart + 0.5) / (mesh->xend - mesh->xstart + 1));
  });
  solver.setState({{"Vort", vort}, {"Eion", Field3D{10.0}}, {"Pe", Field3D{1.0}}});
  auto state = makeState();
  state["species"]["e"]["density"] = Field3D{1.0};
  electrons.transform(state);
  component.transform(state);
  component.finally(state);
  electrons.finally(state);
  const auto output = diagnostics(component);
  const auto rhs = solver.getTimeDerivs();
  const Field3D heat_e = output["electron_viscous_heating"].as<Field3D>();
  EXPECT_GT(integral(output["kinematic_viscous_heating"].as<Field3D>()), 1e-8);
  EXPECT_TRUE(IsFieldEqual(
      heat_e, Field3D{output["kinematic_viscous_heating"].as<Field3D>() / 3673},
      "RGN_NOBNDRY", 1e-12));
  EXPECT_TRUE(
      IsFieldEqual(Field3D{1.5 * rhs["Pe"].as<Field3D>()}, heat_e, "RGN_NOBNDRY", 1e-12));
  const BoutReal wall_power = integral(output["DivEpol"].as<Field3D>());
  EXPECT_NEAR(integral(Field3D{rhs["Eion"].as<Field3D>() + heat_e}), -wall_power, 1e-12);

  // Check the finite-update energy partition without assuming monotonic Kperp:
  // warm-ion pressure-potential wall values also evolve during recovery.
  const BoutReal kinetic_before = integral(output["ion_KE_perp"].as<Field3D>());
  const BoutReal ion_thermal_before =
      integral(Field3D{1.5
                       * (state["species"]["d+"]["pressure"].as<Field3D>()
                          + state["species"]["he++"]["pressure"].as<Field3D>())});
  const BoutReal dt = 1e-4;
  solver.setState({{"Vort", Field3D{vort + dt * rhs["Vort"].as<Field3D>()}},
                   {"Eion", Field3D{10.0 + dt * rhs["Eion"].as<Field3D>()}},
                   {"Pe", Field3D{1.0 + dt * rhs["Pe"].as<Field3D>()}}});
  state = makeState();
  state["species"]["e"]["density"] = Field3D{1.0};
  electrons.transform(state);
  component.transform(state);
  const BoutReal kinetic_after =
      integral(diagnostics(component)["ion_KE_perp"].as<Field3D>());
  const BoutReal ion_thermal_after =
      integral(Field3D{1.5
                       * (state["species"]["d+"]["pressure"].as<Field3D>()
                          + state["species"]["he++"]["pressure"].as<Field3D>())});
  EXPECT_NEAR(ion_thermal_after - ion_thermal_before
                  + dt * (integral(heat_e) + wall_power),
              kinetic_before - kinetic_after, 1e-11);
}

TEST_F(IonEnergyParallelDampingTest, DampingChangesOnlyPhysicalVorticityRate) {
  mesh->getCoordinates()->setNon_uniform(false);
  for (const auto* mode : {"fixed", "free_slip", "no_slip"}) {
    const auto run = [&](bool enabled) {
      FakeSolver solver;
      auto options = makeOptions();
      options["test"]["radial_boundary"] = std::string{mode};
      options["test"]["phi_dissipation"] = enabled;
      options["sound_speed"]["electron_dynamics"] = false;
      IonEnergyTestSystem component("test", options, &solver);
      SoundSpeed sound("sound_speed", options, &solver);
      declareSpecies(component);
      sound.declareAllSpecies(SpeciesInformation({"e"}, {}, {"d+", "he++"}, {}));
      const Field3D vort = makeField<Field3D>([](const Ind3D& i) {
        return 0.001 * (1 + (0.1 * i.x())) * cos(TWOPI * i.z() / 8)
               * (i.y() % 2 ? -1 : 1);
      });
      solver.setState({{"Vort", vort}, {"Eion", Field3D{3.0}}});
      auto state = makeState();
      component.transform(state);
      sound.transform(state);
      component.finally(state);
      auto result = diagnostics(component);
      result["rhs"] = solver.getTimeDerivs()["Vort"].as<Field3D>();
      result["species"] = state["species"].copy();
      return result;
    };
    const auto off = run(false), on = run(true);
    const Field3D damping = on["ion_phi_dissipation"].as<Field3D>();
    EXPECT_GT(max(abs(damping), true), 1e-8);
    EXPECT_TRUE(
        IsFieldEqual(off["ion_phi_dissipation"].as<Field3D>(), 0.0, "RGN_ALL", 0.0));
    EXPECT_TRUE(IsFieldEqual(Field3D{on["rhs"].as<Field3D>() - off["rhs"].as<Field3D>()},
                             damping, "RGN_NOBNDRY", 1e-13));
    for (const auto* field : {"phi", "ion_chi", "DivJpol", "ddt(Eion)"}) {
      EXPECT_TRUE(IsFieldEqual(on[field].as<Field3D>(), off[field].as<Field3D>(),
                               "RGN_NOBNDRY", 0.0))
          << field;
    }
    for (const auto& [name, species] : on["species"].getChildren()) {
      for (const auto* source : {"density_source", "momentum_source", "energy_source"}) {
        if (species.isSet(source)) {
          EXPECT_TRUE(IsFieldEqual(species[source].as<Field3D>(),
                                   off["species"][name][source].as<Field3D>(),
                                   "RGN_NOBNDRY", 0.0))
              << name << ':' << source;
        }
      }
    }
    for (int x : {mesh->xstart - 1, mesh->xend + 1}) {
      for (int y = mesh->ystart; y <= mesh->yend; ++y) {
        for (int z = 0; z < mesh->LocalNz; ++z) {
          EXPECT_EQ(on["rhs"].as<Field3D>()(x, y, z), off["rhs"].as<Field3D>()(x, y, z));
          EXPECT_EQ(damping(x, y, z), 0.0);
        }
      }
    }
  }
}

TEST_F(IonEnergyParallelDampingTest, DampingLossIsRecoveredAsIonThermalEnergy) {
  mesh->periodicX = true;
  static_cast<FakeMesh*>(mesh)->copy_periodic_halos = true;
  mesh->getCoordinates()->setNon_uniform(false);
  FakeSolver solver;
  auto options = makeOptions();
  options["test"]["phi_dissipation"] = true;
  options["sound_speed"]["electron_dynamics"] = false;
  IonEnergyTestSystem component("test", options, &solver);
  SoundSpeed sound("sound_speed", options, &solver);
  declareSpecies(component);
  sound.declareAllSpecies(SpeciesInformation({"e"}, {}, {"d+", "he++"}, {}));
  const Field3D vort = makeField<Field3D>([](const Ind3D& i) {
    return 0.02 * cos(TWOPI * (i.x() - 1.5) / 6) * cos(TWOPI * i.z() / 8)
           * (i.y() % 2 ? -1 : 1);
  });
  solver.setState({{"Vort", vort}, {"Eion", Field3D{3.0}}});
  auto state = makeState();
  component.transform(state);
  sound.transform(state);
  component.finally(state);
  const auto before = diagnostics(component);
  const auto thermal = [&](const Options& value) {
    return integral(Field3D{1.5
                            * (value["species"]["d+"]["pressure"].as<Field3D>()
                               + value["species"]["he++"]["pressure"].as<Field3D>())});
  };
  const BoutReal thermal_before = thermal(state);
  // One small, isolated numerical-damping update at fixed evolved Eion. The
  // off/on test checks the full physical RHS independently; this update tests
  // the recovery partition without contamination from physical transport.
  solver.setState(
      {{"Vort", Field3D{vort + 1e-3 * before["ion_phi_dissipation"].as<Field3D>()}},
       {"Eion", Field3D{3.0}}});
  state = makeState();
  component.transform(state);
  const auto after = diagnostics(component);
  const BoutReal kinetic_loss = integral(before["ion_KE_perp"].as<Field3D>())
                                - integral(after["ion_KE_perp"].as<Field3D>());
  EXPECT_GT(kinetic_loss, 1e-9);
  EXPECT_NEAR(thermal(state) - thermal_before, kinetic_loss, 1e-11);
}

// Boundary and restart behaviour.
TEST_F(EvolveVorticityIonEnergyTest, RestartAndBoundaryHooksPreserveSolverWallState) {
  Options::root()["Vort"]["bndry_all"] = "neumann";
  FakeSolver solver;
  auto options = makeOptions();
  options["test"]["radial_boundary"] = "no_slip";
  options["test"]["phi_boundary_velocity"] = 0.4;
  options["test"]["viscosity"] = 0.3;
  options["polarisation_drift"]["polarisation_corrections"] = 6;
  options["units"]["seconds"].force(0.2); // Normalised velocity is 0.08.
  IonEnergyTestSystem component("test", options, &solver);
  declareSpecies(component);
  ASSERT_TRUE(solver.evolvesBoundary("Vort"));
  Field3D encoded = makeField<Field3D>([&](const Ind3D& i) {
    if (i.x() < mesh->xstart) {
      return 0.08;
    }
    if (i.x() > mesh->xend) {
      return -0.03;
    }
    return 0.002 * sin(2 * PI * i.z() / nz);
  });
  solver.setState({{"Vort", encoded}, {"Eion", Field3D{4.0}}});
  auto packed = solver.packState();
  auto state = makeState();
  solver.applyValueBoundaries(5.0);
  EXPECT_EQ(solver.packState(), packed);
  component.transform(state);
  component.finally(state);
  EXPECT_EQ(solver.packState(), packed);
  const auto before = copy(solver.getTimeDerivs()["Vort"].as<Field3D>());
  solver.applyDerivativeBoundaries();
  EXPECT_TRUE(
      IsFieldEqual(before, solver.getTimeDerivs()["Vort"].as<Field3D>(), "RGN_ALL", 0.0));
  const auto output = diagnostics(component);
  const Field3D phi = output["phi"].as<Field3D>();
  const Field3D physical = state["fields"]["vorticity"].as<Field3D>();
  const Field3D saved = output["Vort"].as<Field3D>();
  EXPECT_TRUE(IsFieldEqual(physical, saved, "RGN_ALL", 0.0));
  for (int y = mesh->ystart; y <= mesh->yend; ++y) {
    for (const auto side : {std::pair{mesh->xstart, mesh->xstart - 1},
                            std::pair{mesh->xend, mesh->xend + 1}}) {
      const auto [interior, guard] = side;
      BoutReal phi0 = 0.0;
      for (int z = 0; z < nz; ++z) {
        phi0 += phi(interior, y, z) / nz;
      }
      const BoutReal rate = 0.08 * 2.0 * (phi0 - encoded(guard, y, 0));
      for (int z = 0; z < nz; ++z) {
        EXPECT_NEAR(0.5 * (phi(guard, y, z) + phi(interior, y, z)), encoded(guard, y, z),
                    1e-12);
        EXPECT_NEAR(before(guard, y, z), rate, 1e-12);
        EXPECT_EQ(output["phi_wall_rate"].as<Field3D>()(guard, y, z),
                  before(guard, y, z));
        EXPECT_NEAR(physical(guard, y, z), physical(interior, y, z), 1e-12);
      }
    }
  }
  // Repeated and earlier RHS evaluations leave the integrator's state intact.
  solver.applyValueBoundaries(1.0);
  component.transform(state = makeState());
  component.finally(state);
  EXPECT_EQ(solver.packState(), packed);
  Options restart{{"phi", Field3D{0.2}}, {"ion_energy_boundary_state_version", 1}};
  component.restartVars(restart);
  EXPECT_EQ(solver.packState(), packed);
}

TEST_F(EvolveVorticityIonEnergyTest, CoreWallUsesWeightedStateAndRate) {
  for (const auto* mode : {"fixed", "free_slip", "no_slip"}) {
    auto* coord = mesh->getCoordinates();
    coord->setJ(makeField<Coordinates::FieldMetric>(
        [](const auto& i) { return 1.0 + 0.1 * i.x() + 0.3 * i.y(); }));
    coord->setDy(makeField<Coordinates::FieldMetric>(
        [](const auto& i) { return 0.5 + 0.2 * i.y(); }));
    coord->setContravariantMetricTensor(
        ContravariantMetricTensor{makeField<Coordinates::FieldMetric>(
                                      [](const auto& i) { return 1.0 + 0.2 * i.y(); }),
                                  1.0, 1.0, 0.0, 0.0, 0.0});
    FakeSolver solver;
    auto options = makeOptions();
    options["test"]["radial_boundary"] = mode;
    IonEnergyTestSystem component("test", options, &solver);
    declareSpecies(component);
    const Field3D restart_phi = makeField<Field3D>([&](const Ind3D& i) {
      return i.y() < mesh->ystart || i.y() > mesh->yend ? 999.0 : 0.01 * (1 + i.y());
    });
    Options restart{{"phi", restart_phi}};
    component.restartVars(restart);
    const int interior = mesh->xstart, guard = interior - 1;
    BoutReal expected_midpoint = 0.0, weight_sum = 0.0;
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      const BoutReal weight = 0.5 * (coord->J()(guard, y, 0) + coord->J()(interior, y, 0))
                              * coord->dy()(interior, y, 0);
      expected_midpoint += weight * restart_phi(interior, y, 0);
      weight_sum += weight;
    }
    expected_midpoint /= weight_sum;
    const Field3D initial = solver.getState()["Vort"].as<Field3D>();
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      EXPECT_NEAR(initial(guard, y, 0), expected_midpoint, 1e-15);
    }
    // A version-1 restart may contain per-Y differences. Recovery must apply
    // only their flux-surface average, without altering the solver state.
    Field3D encoded{0.0};
    for (int x = 0; x < mesh->xstart; ++x) {
      for (int y = mesh->ystart; y <= mesh->yend; ++y) {
        for (int z = 0; z < nz; ++z) {
          encoded(x, y, z) = restart_phi(interior, y, 0);
        }
      }
    }
    solver.setState({{"Vort", encoded}, {"Eion", Field3D{4.0}}});
    auto packed = solver.packState();
    auto state = makeState();
    component.transform(state);
    component.finally(state);
    EXPECT_EQ(packed, solver.packState());
    const auto output = diagnostics(component);
    const Field3D phi = output["phi"].as<Field3D>();
    const Field3D rates = output["phi_wall_rate"].as<Field3D>();
    BoutReal expected_rate = 0.0;
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      const BoutReal midpoint = 0.5 * (phi(guard, y, 0) + phi(interior, y, 0));
      EXPECT_NEAR(midpoint, expected_midpoint, 1e-12);
      const BoutReal weight = 0.5 * (coord->J()(guard, y, 0) + coord->J()(interior, y, 0))
                              * coord->dy()(interior, y, 0);
      const BoutReal phi0 = phi(interior, y, 0);
      expected_rate += weight * 0.4 * sqrt(coord->g11()(interior, y, 0)) * 2.0
                       * (phi0 - expected_midpoint);
    }
    expected_rate = std::string{mode} == "no_slip" ? expected_rate / weight_sum : 0.0;
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      EXPECT_NEAR(rates(guard, y, 0), expected_rate, 1e-12);
    }
  }
}

TEST_F(IonEnergySingletonYTest, RecoveryAndRelaxationWorkWithoutYGuards) {
  FakeSolver solver;
  auto options = makeOptions();
  options["test"]["radial_boundary"] = "no_slip";
  options["test"]["viscosity"] = 0.1;
  IonEnergyTestSystem component("test", options, &solver);
  declareSpecies(component);
  const Field3D vort = makeField<Field3D>([](const Ind3D& i) {
    return i.x() >= mesh->xstart && i.x() <= mesh->xend ? 0.001 * sin(TWOPI * i.z() / nz)
                                                        : 0.0;
  });
  solver.setState({{"Vort", vort}, {"Eion", Field3D{4.0}}});
  auto packed = solver.packState();
  auto state = makeState();
  component.transform(state);
  component.finally(state);
  EXPECT_EQ(solver.packState(), packed);
  const auto output = diagnostics(component);
  const Field3D total_potential =
      output["phi"].as<Field3D>() + output["ion_pressure_potential"].as<Field3D>();
  EXPECT_TRUE(
      IsFieldEqual(bout::ion_energy::Div_a_Grad_perp_XZ(Field3D{1.5}, total_potential),
                   vort, "RGN_NOBNDRY", 1e-10));
  EXPECT_GT(max(abs(output["DivJviscous"].as<Field3D>()), true), 1e-8);
}

TEST_F(EvolveVorticityIonEnergyTest, WallModelsPairParticleChargeAndEnergyBudgets) {
  for (const auto* mode : {"fixed", "free_slip", "no_slip"}) {
    for (const auto* ac : {"neumann", "dirichlet"}) {
      SCOPED_TRACE(std::string{mode} + "/" + ac);
      FakeSolver solver;
      auto options = makeOptions();
      options["test"]["radial_boundary"] = mode;
      options["polarisation_drift"]["polarisation_ac_boundary"] = ac;
      options["test"]["viscosity"] = 0.1;
      IonEnergyTestSystem component("test", options, &solver);
      declareSpecies(component);
      ASSERT_TRUE(solver.evolvesBoundary("Vort"));
      Field3D vort = makeField<Field3D>([&](const Ind3D& i) {
        if (i.x() < mesh->xstart || i.x() > mesh->xend) {
          return 0.01;
        }
        return 0.002 * (1.0 + sin(0.4 * i.x() + 2 * PI * i.z() / nz));
      });
      solver.setState({{"Vort", vort}, {"Eion", Field3D{6.0}}});
      auto state = makeState();
      // Neumann density guards retain the local Z-dependent wall composition.
      state["species"]["he++"]["density"] = makeField<Field3D>([&](const Ind3D& i) {
        const int x = std::clamp(i.x(), mesh->xstart, mesh->xend);
        return 0.25 + 0.03 * cos(0.6 * x + 2 * PI * i.z() / nz);
      });
      component.transform(state);
      component.finally(state);
      const auto output = diagnostics(component);
      const Field3D physical = state["fields"]["vorticity"].as<Field3D>();
      const Field3D psi = output["ion_pressure_potential"].as<Field3D>();
      const Field3D phi = output["phi"].as<Field3D>();
      for (int x : {mesh->xstart, mesh->xend}) {
        const int guard = x == mesh->xstart ? x - 1 : x + 1;
        for (int y = mesh->ystart; y <= mesh->yend; ++y) {
          for (int z = 0; z < nz; ++z) {
            EXPECT_EQ(psi(guard, y, z), psi(x, y, z));
            EXPECT_NEAR((phi(x, y, z) + phi(guard, y, z)) / 2, 0.01, 1e-12);
            EXPECT_EQ(physical(guard, y, z),
                      (std::string{mode} == "free_slip" ? -1 : 1) * physical(x, y, z));
            if (std::string{mode} != "no_slip") {
              EXPECT_EQ(solver.getTimeDerivs()["Vort"].as<Field3D>()(guard, y, z), 0.0);
            }
          }
        }
      }
      EXPECT_NEAR(integral(output["DivEExB"].as<Field3D>()), 0.0, 1e-12);
      for (const auto* species : {"d+", "he++"}) {
        const Field3D loss =
            output[std::string{"particle_pol_wall_loss_"} + species].as<Field3D>();
        EXPECT_NEAR(integral(state["species"][species]["density_source"].as<Field3D>()),
                    -integral(loss), 1e-12);
      }
      const Field3D flux = output["polarisation_flux_x"].as<Field3D>();
      BoutReal wall_power = 0.0, wall_signal = 0.0;
      for (int x : {mesh->xstart - 1, mesh->xend}) {
        const int cell = x == mesh->xend ? x : x + 1;
        for (int y = mesh->ystart; y <= mesh->yend; ++y) {
          for (int z = 0; z < nz; ++z) {
            BoutReal weight = 0.5 * (phi(x, y, z) + phi(x + 1, y, z)) + psi(cell, y, z);
            BoutReal mass_density = 0.0;
            for (const auto* name : {"d+", "he++"}) {
              const auto& species = state["species"][name];
              mass_density += species["AA"].as<BoutReal>()
                              * species["density"].as<Field3D>()(cell, y, z);
            }
            for (const auto* name : {"d+", "he++"}) {
              const auto& species = state["species"][name];
              weight += 1.5 * species["AA"].as<BoutReal>()
                        / species["charge"].as<BoutReal>()
                        * species["pressure"].as<Field3D>()(cell, y, z) / mass_density;
            }
            const BoutReal outward = x == mesh->xend ? flux(x, y, z) : -flux(x, y, z);
            wall_power += outward * weight;
            wall_signal += std::abs(outward);
          }
        }
      }
      EXPECT_GT(wall_signal, 1e-9);
      EXPECT_NEAR(integral(solver.getTimeDerivs()["Eion"].as<Field3D>()), -wall_power,
                  1e-12);
      if (std::string{mode} != "free_slip") {
        EXPECT_NEAR(integral(output["DivJviscous"].as<Field3D>()), 0.0, 1e-12);
        EXPECT_NEAR(integral(output["DivEviscous"].as<Field3D>()), 0.0, 1e-12);
      } else {
        EXPECT_GT(std::abs(integral(output["DivJviscous"].as<Field3D>())), 1e-8);
      }
    }
  }
}

TEST_F(EvolveVorticityIonEnergyTest, AuxiliaryInverseRecoversDCAndNyquist) {
  Options options{{"type", "cyclic"},
                  {"inner_boundary_flags", INVERT_AC_GRAD},
                  {"outer_boundary_flags", INVERT_AC_GRAD}};
  auto inverse = Laplacian::create(&options);
  inverse->setCoefC(Coordinates::FieldMetric{1.5});
  for (int mode : {0, 1, nz / 2}) {
    SCOPED_TRACE(mode);
    Field3D exact = makeField<Field3D>([&](const Ind3D& i) {
      const BoutReal x = i.x() - mesh->xstart + 0.5;
      return x * (x - (mesh->xend - mesh->xstart + 1)) - 0.25
             + (mode ? cos(TWOPI * mode * i.z() / nz) : 0.0);
    });
    bout::ion_energy::applyPolarisationBoundary(exact, PolarisationACBoundary::neumann);
    const Field3D rhs = bout::ion_energy::Div_a_Grad_perp_XZ(Field3D{1.5}, exact);
    const Field3D recovered = bout::ion_energy::fixedPolarisationPotential(
        Field3D{1.5}, rhs, *inverse, 64, PolarisationACBoundary::neumann);
    EXPECT_TRUE(IsFieldEqual(recovered, exact, "RGN_NOBNDRY", 1e-11));
    EXPECT_TRUE(
        IsFieldEqual(bout::ion_energy::Div_a_Grad_perp_XZ(Field3D{1.5}, recovered), rhs,
                     "RGN_NOBNDRY", 1e-11));
  }
}

TEST_F(EvolveVorticityIonEnergyTest, ConservativeComponentsRequireEachOther) {
  FakeSolver solver;
  auto options = makeOptions();
  options["test"]["type"] = "evolve_vorticity_ion_energy";
  options["hermes"]["components"] = "test";
  EXPECT_THROW(EvolveVorticityIonEnergy("test", options, &solver), BoutException);
  options["hermes"]["components"].force("polarisation_drift");
  EXPECT_THROW(PolarisationDrift("polarisation_drift", options, &solver), BoutException);
}

class IonEnergyPeriodicCurrentTest : public FakeMeshFixture_tmpl<10, 5, 8, false, 2> {
protected:
  void SetUp() override {
    Laplacian::cleanup();
    mesh->periodicX = true;
  }
  void TearDown() override { Laplacian::cleanup(); }
};

TEST_F(IonEnergyPeriodicCurrentTest, RejectsUnsolvableMeanWithoutDiscardingCurrent) {
  for (bool zero_mean : {false, true}) {
    FakeSolver solver;
    Options options{{"units",
                     {{"eV", 1.0},
                      {"inv_meters_cubed", 1.0},
                      {"seconds", 1.0},
                      {"Tesla", 1.0},
                      {"meters", 1.0}}},
                    {"test", {{"phi_rtol", 1e-10}}}};
    Options::root()["mesh"]["paralleltransform"]["type"] = "identity";
    IonEnergyTestSystem component("test", options, &solver);
    component.declareAllSpecies(SpeciesInformation({"e"}, {}, {"d+"}, {}));
    solver.setState({{"Vort", Field3D{0.0}}, {"Eion", Field3D{1.5}}});
    Options state{
        {"species",
         {{"e", {{"charge", -1.0}, {"AA", 1. / 1836}, {"pressure", Field3D{1.0}}}},
          {"d+", {{"charge", 1.0}, {"AA", 2.0}, {"density", Field3D{1.0}}}}}}};
    state["fields"]["DivJextra"] = makeField<Field3D>(
        [&](const Ind3D& i) { return zero_mean ? 0.1 * sin(TWOPI * i.z() / nz) : 0.1; });
    if (zero_mean) {
      EXPECT_NO_THROW(component.transform(state));
    } else {
      EXPECT_THROW(component.transform(state), BoutException);
    }
  }
}

// Repeated evaluations and production component scheduling.
// Note: RHS is only independent of history if phi_reset_initial_guess = true
TEST_F(EvolveVorticityIonEnergyTest, RHSDoesNotDependOnEvaluationHistory) {
  for (const bool relax : {false, true}) {
    FakeSolver solver;
    auto options = makeOptions();
    options["test"]["radial_boundary"] = relax ? "no_slip" : "fixed";
    // Reset phi each time rhs is called
    options["test"]["phi_reset_initial_guess"] = true;
    IonEnergyTestSystem component("test", options, &solver);
    declareSpecies(component);
    const auto evaluate = [&](int phase) {
      solver.setState({{"Vort", makeField<Field3D>([&](const Ind3D& i) {
                          return 0.01 * sin(2 * PI * i.z() / nz + phase);
                        })},
                       {"Eion", Field3D{3.0}}});
      auto state = makeState();
      component.transform(state);
      component.finally(state);
      return state;
    };
    auto state = evaluate(0);
    const Field3D phi = copy(state["fields"]["phi"].as<Field3D>());
    const Field3D rhs = copy(solver.getTimeDerivs()["Eion"].as<Field3D>());
    evaluate(1);
    state = evaluate(0);
    EXPECT_TRUE(
        IsFieldEqual(phi, state["fields"]["phi"].as<Field3D>(), "RGN_NOBNDRY", 0.0));
    EXPECT_TRUE(IsFieldEqual(rhs, solver.getTimeDerivs()["Eion"].as<Field3D>(),
                             "RGN_NOBNDRY", 0.0));
  }
}

TEST_F(EvolveVorticityIonEnergyTest, ConservativeStateAndSourcesDoNotDependOnDiagnose) {
  Options reference;
  for (bool diagnose : {true, false}) {
    FakeSolver solver;
    auto options = makeOptions();
    options["test"]["diagnose"].force(diagnose);
    IonEnergyTestSystem component("test", options, &solver);
    declareSpecies(component);
    solver.setState(
        {{"Vort", Field3D{0.0}}, {"Eion", makeField<Field3D>([&](const Ind3D& i) {
                                    return 3.0 + 0.2 * sin(2 * PI * i.z() / nz);
                                  })}});
    auto state = makeState();
    component.transform(state);
    component.finally(state);
    Options actual;
    actual["phi"] = state["fields"]["phi"].as<Field3D>();
    actual["ion_chi"] = state["fields"]["ion_chi"].as<Field3D>();
    for (const auto* name : {"d+", "he++"}) {
      for (const auto* field :
           {"pressure", "density_source", "momentum_source", "energy_source"}) {
        if (state["species"][name].isSet(field)) {
          actual[std::string{name} + "/" + field] =
              state["species"][name][field].as<Field3D>();
        }
      }
    }
    actual["ddt_Eion"] = solver.getTimeDerivs()["Eion"].as<Field3D>();
    actual["ddt_Vort"] = solver.getTimeDerivs()["Vort"].as<Field3D>();
    if (diagnose) {
      reference = actual.copy();
    } else {
      for (const auto& [name, value] : actual.getChildren()) {
        EXPECT_TRUE(IsFieldEqual(value.as<Field3D>(), reference[name].as<Field3D>(),
                                 "RGN_NOBNDRY", 0.0))
            << name;
      }
    }
  }
}

TEST_F(EvolveVorticityIonEnergyTest, ProductionSheathAutosortsBeforeConservativeDrift) {
  static_cast<FakeMesh*>(mesh)->ix_separatrix = -1;
  Options reference;
  for (bool reversed : {false, true}) {
    FakeSolver solver;
    auto options = makeOptions();
    options["components"] =
        reversed ? "polarisation_drift, sheath_boundary, test, he++, d+, e"
                 : "e, d+, he++, test, sheath_boundary, polarisation_drift";
    options["test"]["type"] = "evolve_vorticity_ion_energy";
    options["e"]["type"] = "quasineutral, evolve_pressure, evolve_momentum";
    options["e"]["charge"] = -1.0;
    options["e"]["AA"] = 1. / 1836;
    options["e"]["thermal_conduction"] = false;
    options["e"]["poloidal_flows"] = false;
    for (const auto* name : {"d+", "he++"}) {
      options[name]["type"] = "evolve_density, evolve_momentum";
      options[name]["charge"] = std::string{name} == "d+" ? 1.0 : 2.0;
      options[name]["AA"] = std::string{name} == "d+" ? 2.0 : 4.0;
      options[name]["poloidal_flows"] = false;
    }
    auto scheduler = ComponentScheduler::create(options, options, &solver);
    solver.setState({{"Vort", Field3D{0.02}},
                     {"Eion", Field3D{1.5}},
                     {"Pe", Field3D{1.0}},
                     {"Nd+", Field3D{0.5}},
                     {"Nhe++", Field3D{0.25}},
                     {"NVd+", Field3D{0.0}},
                     {"NVhe++", Field3D{0.0}},
                     {"NVe", Field3D{0.0}}});
    Options state{{"time", 0.0}, {"fastest_wave", Field3D{4.0}}};
    scheduler->transform(state);
    EXPECT_GT(max(abs(state["fields"]["DivJpar"].as<Field3D>()), true), 1e-4);
    const auto& fields = state["fields"];
    EXPECT_TRUE(IsFieldEqual(
        Field3D{fields["DivJpol"].as<Field3D>() + fields["DivJdia"].as<Field3D>()
                + fields["DivJpar"].as<Field3D>()},
        fields["ion_auxiliary_residual"].as<Field3D>(), "RGN_NOBNDRY", 1e-10));
    Options actual = solver.getTimeDerivs().copy();
    for (const auto* name : {"phi", "ion_chi", "DivJpol", "DivJpar", "ion_current_rhs",
                             "ion_auxiliary_residual"}) {
      actual[name] = state["fields"][name].as<Field3D>();
    }
    if (!reversed) {
      reference = actual.copy();
    } else {
      for (const auto& [name, value] : actual.getChildren()) {
        EXPECT_TRUE(IsFieldEqual(value.as<Field3D>(), reference[name].as<Field3D>(),
                                 "RGN_NOBNDRY", 0.0))
            << name;
      }
    }
  }
}

// Reynolds current: conservation, work, symmetry and continuum convergence.
TEST_F(EvolveVorticityIonEnergyTest, ReynoldsObeysGaussAndBoundaryWork) {
  using namespace bout::ion_energy;
  auto* coord = mesh->getCoordinates();
  coord->setJ(makeField<Coordinates::FieldMetric>(
      [](const auto& i) { return 1.2 + 0.05 * i.x(); }));
  coord->setContravariantMetricTensor(
      ContravariantMetricTensor{1.2, 1.0, 1.4, 0.0, 0.3, 0.0});
  Field3D phi = makeField<Field3D>(
      [](const Ind3D& i) { return sin(0.7 * i.x()) * cos(2 * PI * i.z() / nz); });
  Field3D psi = makeField<Field3D>([](const Ind3D& i) {
    return 1.0 + 0.1 * cos(0.8 * i.x()) * sin(2 * PI * i.z() / nz);
  });
  applyVorticityWallBoundary(phi, true);
  applyWallPressureBoundary(psi);
  const Field3D c = makeField<Field3D>(
      [](const Ind3D& i) { return 1.0 + 0.1 * cos(i.x() + 2 * PI * i.z() / nz); });
  const auto flux = reynoldsFaceFluxes(c, Div_a_Grad_perp_XZ(c, phi + psi), phi, psi);
  std::mt19937 generator(1729);
  std::uniform_real_distribution<BoutReal> distribution(0.1, 2.0);
  const Field3D weight =
      makeField<Field3D>([&](const Ind3D&) { return distribution(generator); });
  Field3D changed = copy(weight);
  for (const auto& i : changed.getRegion("RGN_ALL")) {
    if (i.x() < mesh->xstart || i.x() > mesh->xend) {
      changed[i] = 100 * distribution(generator);
    }
  }
  EXPECT_TRUE(IsFieldEqual(divergence(flux, &weight), divergence(flux, &changed),
                           "RGN_NOBNDRY", 0.0));
  BoutReal boundary_work = 0, weighted_outward = 0, signal = 0;
  for (int x : {mesh->xstart - 1, mesh->xend}) {
    const int interior = x == mesh->xend ? x : x + 1;
    const BoutReal sign = x == mesh->xend ? 1 : -1;
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      BoutReal charge = 0;
      for (int z = 0; z < nz; ++z) {
        const auto F = flux.x(x, y, z);
        charge += F;
        signal += std::abs(F);
        boundary_work +=
            sign * F * 0.5
            * (phi(x, y, z) + psi(x, y, z) + phi(x + 1, y, z) + psi(x + 1, y, z));
        weighted_outward += sign * F * weight(interior, y, z);
      }
      EXPECT_TRUE(std::isfinite(charge));
    }
  }
  EXPECT_GT(signal, 1e-6);
  EXPECT_NEAR(integral(divergence(flux, &weight)), weighted_outward, 1e-12);
  EXPECT_NEAR(integral(faceWork(flux, phi + psi)),
              boundary_work - integral((phi + psi) * divergence(flux)), 1e-12);
  const Field3D complement = 1.0 - weight;
  EXPECT_LT(max(abs(divergence(flux, &weight) + divergence(flux, &complement)
                    - divergence(flux)),
                true),
            1e-13);
}

TEST_F(EvolveVorticityIonEnergyTest, ReynoldsMirrorSymmetry) {
  using namespace bout::ion_energy;
  auto* coord = mesh->getCoordinates();
  coord->setNon_uniform(false);
  const Field3D phi = makeField<Field3D>(
      [](const Ind3D& i) { return 0.2 * sin(0.4 * i.x()) * cos(2 * PI * i.z() / nz); });
  const Field3D psi = makeField<Field3D>([](const Ind3D& i) {
    return 1.0 + 0.1 * cos(0.3 * i.x()) * sin(2 * PI * i.z() / nz);
  });
  const Field3D vort = Div_a_Grad_perp_XZ(Field3D{1.0}, phi + psi);
  const auto reflect = [&](const Field3D& f) {
    return makeField<Field3D>(
        [&](const Ind3D& i) { return f(i.x(), i.y(), (nz - i.z()) % nz); });
  };
  const auto original = reynoldsFaceFluxes(Field3D{1.0}, vort, phi, psi);
  const auto mirrored =
      reynoldsFaceFluxes(Field3D{1.0}, -reflect(vort), -reflect(phi), -reflect(psi));
  for (const auto& i : phi.getRegion("RGN_NOBNDRY")) {
    EXPECT_NEAR(mirrored.x[i], -original.x(i.x(), i.y(), (nz - i.z()) % nz), 1e-12);
    EXPECT_NEAR(mirrored.z[i], original.z(i.x(), i.y(), (nz - i.z() - 1) % nz), 1e-12);
  }
}

// The analytic reference is generated from covariant derivatives with
// Christoffel symbols, independently of the discrete current.
template <int NX, int NZ>
class ReynoldsConvergenceMesh : public FakeMeshFixture_tmpl<NX, 1, NZ, false, 2> {
public:
  void TestBody() override {}

  BoutReal error() {
    using namespace bout::ion_energy;
    const BoutReal h = 1.0 / NZ;
    Options::root()["h"] = h;
    Options::root()["r"] = std::string{"1+h*(x-1.5)"};
    // Parsed option references retain constants; each refinement needs its h.
    FieldFactory::get()->cleanCache();
    auto* coord = mesh->getCoordinates();
    coord->setNon_uniform(false);
    coord->setDx(h);
    coord->setDz(2 * PI / NZ);
    const auto metric = [&](const char* expression) {
      return Coordinates::FieldMetric{
          FieldFactory::get()->create2D(expression, nullptr, mesh)};
    };
    coord->setJ(metric(reynolds_reference::J));
    coord->setBxy(metric("1/r"));
    coord->setContravariantMetricTensor(ContravariantMetricTensor{
        metric(reynolds_reference::g11), Coordinates::FieldMetric{1.0},
        metric(reynolds_reference::g33), Coordinates::FieldMetric{0.0},
        metric(reynolds_reference::g13), Coordinates::FieldMetric{0.0}});
    const auto field = [&](const char* expression) {
      return FieldFactory::get()->create3D(expression, nullptr, mesh);
    };
    const Field3D phi = field(reynolds_reference::phi);
    const Field3D psi = field(reynolds_reference::psi);
    const Field3D expected = field(reynolds_reference::expected);
    const Field3D actual = divergence(reynoldsFaceFluxes(
        field(reynolds_reference::c), field(reynolds_reference::omega), phi, psi));
    BoutReal sum = 0;
    int count = 0;
    for (const auto& i : expected.getRegion("RGN_NOBNDRY")) {
      if (i.x() >= mesh->xstart + 3 && i.x() <= mesh->xend - 3) {
        sum += abs(actual[i] - expected[i]);
        ++count;
      }
    }
    return sum / count;
  }
};

TEST(IonEnergyReynoldsConvergence, LimitedCurrentConvergesInShearedAnnulus) {
  const BoutReal coarse = ReynoldsConvergenceMesh<36, 32>{}.error();
  const BoutReal fine = ReynoldsConvergenceMesh<68, 64>{}.error();
  const BoutReal finest = ReynoldsConvergenceMesh<132, 128>{}.error();
  ASSERT_GT(finest, 0.0);
  EXPECT_GE(std::log2(coarse / fine), 1.5);
  EXPECT_GE(std::log2(fine / finest), 1.5);
}
