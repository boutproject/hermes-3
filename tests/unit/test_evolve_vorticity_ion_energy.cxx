#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include "test_extras.hxx"

#include "../../include/div_ops.hxx"
#include "../../include/evolve_pressure.hxx"
#include "../../include/evolve_vorticity_ion_energy.hxx"

#include <bout/constants.hxx>
#include <bout/derivs.hxx>
#include <bout/difops.hxx>

using bout::globals::mesh;

class EvolveVorticityIonEnergyTest : public FakeMeshFixture_tmpl<10, 5, 8, false, 2> {
protected:
  // Delp2 uses a global Laplacian that must not outlive the fixture's mesh.
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
            {"test",
             {{"average_atomic_mass", 1.5},
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

  void declareSpecies(EvolveVorticityIonEnergy& component) {
    component.declareAllSpecies(SpeciesInformation({"e"}, {}, {"d+", "he++"}, {}));
  }

  Options diagnostics(EvolveVorticityIonEnergy& component) {
    Options output{{"Nnorm", 1.0}, {"Tnorm", 1.0}, {"Omega_ci", 1.0}, {"rho_s0", 1.0}};
    component.outputVars(output);
    return output;
  }
};

TEST_F(EvolveVorticityIonEnergyTest, FreshStartRestartHookPreservesRecoveredPotential) {
  FakeSolver solver;
  auto options = makeOptions();
  EvolveVorticityIonEnergy component("test", options, &solver);
  declareSpecies(component);
  // postInit saves the unsolved potential before the initial RHS evaluation.
  Options restart;
  component.restartVars(restart);
  const BoutReal k = 2 * PI / nz;
  const Field3D energy =
      makeField<Field3D>([&](const Ind3D& i) { return 3.0 + 0.2 * cos(k * i.z()); });
  solver.setState({{"Vort", Field3D{0.0}}, {"Eion", energy}});
  auto state = makeState();
  component.transform(state);
  const Field3D recovered = copy(diagnostics(component)["phi"].as<Field3D>());
  EXPECT_GT(max(abs(recovered), false), 1e-3);
  // The initial output hook must save, rather than restore, this solution.
  component.restartVars(restart);
  EXPECT_TRUE(
      IsFieldEqual(restart["phi"].as<Field3D>(), recovered, "RGN_NOBNDRY", 1e-12));
  EXPECT_TRUE(IsFieldEqual(diagnostics(component)["phi"].as<Field3D>(), recovered,
                           "RGN_NOBNDRY", 1e-12));
}

TEST_F(EvolveVorticityIonEnergyTest, RestartRestorationIsIndependentForEachInstance) {
  auto options = makeOptions();
  FakeSolver first_solver;
  EvolveVorticityIonEnergy first("test", options, &first_solver);
  Options first_restart{{"phi", Field3D{2.0}}};
  first.restartVars(first_restart);
  first_restart["phi"] = Field3D{7.0};
  first.restartVars(first_restart);
  EXPECT_TRUE(
      IsFieldEqual(first_restart["phi"].as<Field3D>(), 2.0, "RGN_NOBNDRY", 1e-12));

  FakeSolver second_solver;
  EvolveVorticityIonEnergy second("test", options, &second_solver);
  Options second_restart{{"phi", Field3D{3.0}}};
  second.restartVars(second_restart);
  EXPECT_TRUE(
      IsFieldEqual(second_restart["phi"].as<Field3D>(), 3.0, "RGN_NOBNDRY", 1e-12));
}

TEST_F(EvolveVorticityIonEnergyTest, CommonTemperatureUsesMassToChargePressurePotential) {
  FakeSolver solver;
  auto options = makeOptions();
  EvolveVorticityIonEnergy component("test", options, &solver);
  declareSpecies(component);
  solver.setState({{"Vort", Field3D{0.0}}, {"Eion", Field3D{3.0}}});
  auto state = makeState();

  component.transform(state);
  component.finally(state);
  auto output = diagnostics(component);

  // There is no flow, so Eion = 3/2 * (n_d + n_he) * Ti.
  EXPECT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["d+"]["temperature"]), 8. / 3,
                           "RGN_NOBNDRY", 1e-12));
  EXPECT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["he++"]["temperature"]), 8. / 3,
                           "RGN_NOBNDRY", 1e-12));
  EXPECT_TRUE(IsFieldEqual(output["ion_pressure_potential"].as<Field3D>(), 8. / 3,
                           "RGN_NOBNDRY", 1e-12));
  EXPECT_TRUE(IsFieldEqual(solver.getTimeDerivs()["Eion"].as<Field3D>(), 0.0,
                           "RGN_NOBNDRY", 1e-12));
}

TEST_F(EvolveVorticityIonEnergyTest, IonCurrentUsesIonMomentumWithoutElectronMomentum) {
  FakeSolver solver;
  auto options = makeOptions();
  EvolveVorticityIonEnergy component("test", options, &solver);
  declareSpecies(component);
  const Field3D vort = makeField<Field3D>([](const Ind3D& i) { return 0.02 * i.y(); });
  solver.setState({{"Vort", vort}, {"Eion", Field3D{10.0}}});
  auto state = makeState();
  state["species"]["d+"]["velocity"] = Field3D{1.0};
  state["species"]["d+"]["momentum"] = Field3D{1.0};
  state["species"]["he++"]["velocity"] = Field3D{-0.5};
  state["species"]["he++"]["momentum"] = Field3D{-0.5};

  ASSERT_NO_THROW(component.transform(state));
  component.finally(state);
  auto output = diagnostics(component);

  // Jpar_i = (1/2)*1 + (2/4)*(-1/2) = 1/4, independently of electrons.
  const Field3D expected = 0.25 * Div_par(output["phi"].as<Field3D>());
  EXPECT_GT(max(abs(expected), false), 1e-5);
  EXPECT_TRUE(
      IsFieldEqual(output["DivEelectric"].as<Field3D>(), expected, "RGN_NOBNDRY", 1e-12));
  const Field3D energy_balance =
      solver.getTimeDerivs()["Eion"].as<Field3D>() + output["DivEExB"].as<Field3D>()
      + output["DivEReynolds"].as<Field3D>() + output["DivEpol"].as<Field3D>()
      + output["DivEparallel"].as<Field3D>() + expected;
  EXPECT_TRUE(IsFieldEqual(energy_balance, 0.0, "RGN_NOBNDRY", 1e-12));
}

TEST_F(EvolveVorticityIonEnergyTest, ElectronCurrentCanReadElectronMass) {
  FakeSolver solver;
  auto options = makeOptions();
  EvolveVorticityIonEnergy component("test", options, &solver);
  declareSpecies(component);
  solver.setState({{"Vort", Field3D{0.0}}, {"Eion", Field3D{3.0}}});
  auto state = makeState();
  state["species"]["e"]["momentum"] = Field3D{0.0};

  EXPECT_NO_THROW(component.transform(state));
}

TEST_F(EvolveVorticityIonEnergyTest, ParallelFluxCarriesEnthalpyAndKineticEnergy) {
  static_cast<FakeMesh*>(mesh)->ix_separatrix = -1; // Open ends for the linear profile.
  FakeSolver solver;
  auto options = makeOptions();
  options["test"]["average_atomic_mass"].force(2.0);
  EvolveVorticityIonEnergy component("test", options, &solver);
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
      makeField<Field3D>([](const Ind3D& i) { return 1.0 + 0.1 * i.y(); });
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

TEST_F(EvolveVorticityIonEnergyTest, ParticleRemovalCarriesItsKineticEnergy) {
  FakeSolver solver;
  auto options = makeOptions();
  EvolveVorticityIonEnergy component("test", options, &solver);
  declareSpecies(component);
  solver.setState({{"Vort", Field3D{0.0}}, {"Eion", Field3D{10.0}}});
  auto state = makeState();
  auto& ion = state["species"]["d+"];
  ion["velocity"] = Field3D{3.0};
  // Remove one ion per unit time at its existing velocity and temperature.
  ion["density_source"] = Field3D{-1.0};
  ion["momentum_source"] = Field3D{-6.0};
  ion["energy_source"] = Field3D{-4.0};

  component.transform(state);
  component.finally(state);
  auto output = diagnostics(component);

  // Each removed ion carries 4 thermal + (1/2)*2*3^2 = 13 total energy.
  EXPECT_TRUE(
      IsFieldEqual(output["Eion_src"].as<Field3D>(), -13.0, "RGN_NOBNDRY", 1e-12));
}

TEST_F(EvolveVorticityIonEnergyTest, ParticleInjectionConservesTransferredEnergy) {
  FakeSolver solver;
  auto options = makeOptions();
  EvolveVorticityIonEnergy component("test", options, &solver);
  declareSpecies(component);
  solver.setState({{"Vort", Field3D{0.0}}, {"Eion", Field3D{10.0}}});
  auto state = makeState();
  auto& ion = state["species"]["d+"];
  ion["velocity"] = Field3D{1.0};
  ion["density_source"] = Field3D{1.0};
  ion["momentum_source"] = Field3D{6.0};
  // Incident velocity 3: reaction supplies (1/2)*2*(3-1)^2 mixing heat.
  ion["energy_source"] = Field3D{4.0};

  component.transform(state);
  component.finally(state);
  auto output = diagnostics(component);

  EXPECT_TRUE(IsFieldEqual(output["Eion_src"].as<Field3D>(), 9.0, "RGN_NOBNDRY", 1e-12));
}

TEST_F(EvolveVorticityIonEnergyTest, PolarisationPotentialSatisfiesVorticityDerivative) {
  static_cast<FakeMesh*>(mesh)->setGridDataSource(
      new FakeGridDataSource{{{"bxcvx", 0.0}, {"bxcvy", 0.0}, {"bxcvz", 0.5}}});
  FakeSolver solver;
  auto options = makeOptions();
  EvolveVorticityIonEnergy component("test", options, &solver);
  declareSpecies(component);
  solver.setState({{"Vort", Field3D{0.0}}, {"Eion", Field3D{3.0}}});
  auto state = makeState();
  const BoutReal k = 2 * PI / nz;
  state["species"]["e"]["pressure"] =
      makeField<Field3D>([&](const Ind3D& i) { return 1.0 + 0.2 * cos(k * i.z()); });

  component.transform(state);
  component.finally(state);
  auto output = diagnostics(component);
  const Field3D chi = output["ddt(Phi)"].as<Field3D>();
  const Field3D rhs = solver.getTimeDerivs()["Vort"].as<Field3D>();
  EXPECT_GT(max(abs(rhs), false), 1e-5);

  // The inversion converges to the finite-volume forward operator, even
  // though its approximate inverse uses spectral Z differentiation.
  const Field3D forward = bout::ion_energy::Div_a_Grad_perp_XZ(Field3D{1.5}, chi);
  for (int x = mesh->xstart + 1; x < mesh->xend; ++x) {
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      for (int z = 0; z < nz; ++z) {
        EXPECT_NEAR(forward(x, y, z), rhs(x, y, z), 2e-10);
      }
    }
  }
}

TEST_F(EvolveVorticityIonEnergyTest, PolarisationFluxConservesEnergyOnNonuniformMesh) {
  static_cast<FakeMesh*>(mesh)->setGridDataSource(
      new FakeGridDataSource{{{"bxcvx", 0.0}, {"bxcvy", 0.0}, {"bxcvz", 0.5}}});
  auto* coord = mesh->getCoordinates();
  const Field2D dx =
      makeField<Field2D>([](const Ind2D& i) { return 1.0 + 0.05 * i.x() * i.x(); });
  const Field2D dz{1.0}; // The spectral inversion requires a uniform Z length.
  coord->setDx(dx);
  coord->setDz(dz);
  FakeSolver solver;
  auto options = makeOptions();
  EvolveVorticityIonEnergy component("test", options, &solver);
  declareSpecies(component);
  const BoutReal k = 2 * PI / nz;
  const Field3D energy =
      makeField<Field3D>([&](const Ind3D& i) { return 3.0 + 0.2 * sin(k * i.z()); });
  solver.setState({{"Vort", Field3D{0.0}}, {"Eion", energy}});
  auto state = makeState();
  state["species"]["e"]["pressure"] =
      makeField<Field3D>([&](const Ind3D& i) { return 1.0 + 0.2 * sin(k * i.z()); });

  component.transform(state);
  component.finally(state);
  auto output = diagnostics(component);
  const Field3D div_flux = output["DivEpol"].as<Field3D>();
  EXPECT_GT(max(abs(div_flux), false), 1e-5);
  BoutReal integral = 0.0;
  for (const auto& i : div_flux.getRegion("RGN_NOBNDRY")) {
    integral += div_flux[i] * dx[i] * dz[i];
  }
  // Zero radial flux and periodic Z must give zero net energy transport,
  // even when radial cell widths vary with radius.
  EXPECT_NEAR(integral, 0.0, 1e-11);

  // Integrate each shared face's kinetic energy once. The result must agree
  // with the sum of the cell energies, including half-faces at radial boundaries.
  const Field3D Phi =
      output["phi"].as<Field3D>() + output["ion_pressure_potential"].as<Field3D>();
  const Field3D kinetic_energy = output["ion_KE_perp"].as<Field3D>();
  BoutReal cell_energy = 0.0;
  BoutReal face_energy = 0.0;
  for (const auto& i : kinetic_energy.getRegion("RGN_NOBNDRY")) {
    cell_energy += kinetic_energy[i] * dx[i] * dz[i];
    face_energy += 0.5 * 1.5 * dx[i] * SQ(Phi[i.zp()] - Phi[i]) / dz[i];
    const BoutReal right_weight = i.x() == mesh->xend ? 0.25 : 0.5;
    face_energy += right_weight * 1.5 * dz[i] * SQ(Phi[i.xp()] - Phi[i])
                   / (0.5 * (dx[i] + dx[i.xp()]));
    if (i.x() == mesh->xstart) {
      face_energy +=
          0.25 * 1.5 * dz[i] * SQ(Phi[i] - Phi[i.xm()]) / (0.5 * (dx[i] + dx[i.xm()]));
    }
  }
  EXPECT_NEAR(cell_energy, face_energy, 1e-9);
}

TEST_F(EvolveVorticityIonEnergyTest, CrossMetricHasExpectedMixedDerivative) {
  auto* coord = mesh->getCoordinates();
  coord->setContravariantMetricTensor(
      ContravariantMetricTensor{1.0, 1.0, 1.0, 0.0, 0.3, 0.0});
  const BoutReal kx = 2 * PI / (mesh->xend - mesh->xstart + 1);
  const BoutReal kz = 2 * PI / nz;
  const Field3D f = makeField<Field3D>(
      [&](const Ind3D& i) { return cos(kx * (i.x() - mesh->xstart) + kz * i.z()); });
  const Field3D forward = bout::ion_energy::Div_a_Grad_perp_XZ(Field3D{1.0}, f);
  const BoutReal eigenvalue =
      -4 * SQ(sin(kx / 2)) - 4 * SQ(sin(kz / 2)) - 2 * 0.3 * sin(kx) * sin(kz);
  EXPECT_TRUE(IsFieldEqual(forward, Field3D{eigenvalue * f}, "RGN_NOBNDRY", 1e-12));
  const Field3D ke = bout::ion_energy::kineticEnergyPerp(Field3D{1.0}, f);
  BoutReal integral = 0.0;
  for (const auto& i : ke.getRegion("RGN_NOBNDRY")) {
    integral += ke[i];
  }
  EXPECT_NEAR(integral, -0.25 * eigenvalue * forward.getRegion("RGN_NOBNDRY").size(),
              1e-11);
}

TEST_F(EvolveVorticityIonEnergyTest, IntegrationByPartsOnNonuniformCrossMetricMesh) {
  auto* coord = mesh->getCoordinates();
  const BoutReal kx = 2 * PI / (mesh->xend - mesh->xstart + 1);
  const BoutReal kz = 2 * PI / nz;
  // Profiles, coefficients and guard cells are periodic in X and Z. In a
  // metric-3D build this also exercises metrics varying in Z.
  const auto phase = [&](const auto& i) {
    BoutReal value = kx * (i.x() - mesh->xstart);
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

TEST_F(EvolveVorticityIonEnergyTest, PotentialRecoverySatisfiesFiniteVolumeVorticity) {
  auto* coord = mesh->getCoordinates();
  coord->setContravariantMetricTensor(
      ContravariantMetricTensor{1.0, 1.0, 1.0, 0.0, 0.3, 0.0});
  // Avoid the fixture's unrelated nonuniform first-derivative correction.
  coord->setNon_uniform(false);
  FakeSolver solver;
  auto options = makeOptions();
  EvolveVorticityIonEnergy component("test", options, &solver);
  declareSpecies(component);
  const BoutReal kz = 2 * PI / nz;
  const Field3D Phi = makeField<Field3D>([&](const Ind3D& i) {
    return 0.05 * sin(0.4 * (i.x() - mesh->xstart) + 3 * kz * i.z());
  });
  const Field3D coefficient{Coordinates::FieldMetric{1.5 / SQ(coord->Bxy())}};
  const Field3D vort = bout::ion_energy::Div_a_Grad_perp_XZ(coefficient, Phi);
  const Field3D energy =
      makeField<Field3D>([&](const Ind3D& i) { return 3.0 + 0.1 * cos(2 * kz * i.z()); });
  solver.setState({{"Vort", vort}, {"Eion", energy}});
  auto state = makeState();
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

TEST_F(EvolveVorticityIonEnergyTest, PotentialRecoveryReportsFailureToConverge) {
  FakeSolver solver;
  auto options = makeOptions();
  options["test"]["phi_max_iterations"].force(1);
  EvolveVorticityIonEnergy component("test", options, &solver);
  declareSpecies(component);
  const Field3D vort =
      makeField<Field3D>([&](const Ind3D& i) { return 0.05 * cos(6 * PI * i.z() / nz); });
  solver.setState({{"Vort", vort}, {"Eion", Field3D{3.0}}});
  auto state = makeState();
  EXPECT_THROW(component.transform(state), BoutException);
}

TEST_F(EvolveVorticityIonEnergyTest,
       ElectronThermalSourceExcludesParallelElectricalWork) {
  static_cast<FakeMesh*>(mesh)->ix_separatrix = -1;
  FakeSolver solver;
  auto options = makeOptions();
  options["e"]["poloidal_flows"] = false;
  options["e"]["bndry_flux"] = false;
  EvolvePressure electrons("e", options, &solver);
  EvolveVorticityIonEnergy component("test", options, &solver);
  declareSpecies(component);
  electrons.declareAllSpecies(SpeciesInformation({"e"}, {}, {"d+", "he++"}, {}));
  const Field3D vort = makeField<Field3D>([](const Ind3D& i) { return 0.02 * i.y(); });
  solver.setState({{"Vort", vort}, {"Eion", Field3D{3.0}}, {"Pe", Field3D{1.0}}});
  auto state = makeState();
  state["species"]["e"]["density"] = Field3D{1.0};
  state["species"]["e"]["momentum"] =
      makeField<Field3D>([](const Ind3D& i) { return -(1.0 + 0.1 * i.y()) / 1836; });
  electrons.transform(state);
  component.transform(state);
  component.finally(state);
  electrons.finally(state);
  const auto output = diagnostics(component);
  // A varying parallel current gives nonzero total energy exchange.
  // Constant Pe and zero curvature need no thermal heating.
  EXPECT_GT(max(abs(output["phi_Div_Je"].as<Field3D>()), false), 1e-5);
  EXPECT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["e"]["energy_source"]), 0.0,
                           "RGN_NOBNDRY", 1e-12));
  EXPECT_TRUE(IsFieldEqual(solver.getTimeDerivs()["Pe"].as<Field3D>(), 0.0, "RGN_NOBNDRY",
                           1e-12));
}

TEST_F(EvolveVorticityIonEnergyTest, ElectronPressureReceivesExBCompressionWork) {
  static_cast<FakeMesh*>(mesh)->setGridDataSource(
      new FakeGridDataSource{{{"bxcvx", 0.0}, {"bxcvy", 0.0}, {"bxcvz", 0.5}}});
  FakeSolver solver;
  auto options = makeOptions();
  options["e"]["poloidal_flows"] = false;
  options["e"]["bndry_flux"] = false;
  EvolvePressure electrons("e", options, &solver);
  EvolveVorticityIonEnergy component("test", options, &solver);
  declareSpecies(component);
  electrons.declareAllSpecies(SpeciesInformation({"e"}, {}, {"d+", "he++"}, {}));
  const Field3D vort =
      makeField<Field3D>([&](const Ind3D& i) { return 0.03 * cos(2 * PI * i.z() / nz); });
  solver.setState({{"Vort", vort}, {"Eion", Field3D{3.0}}, {"Pe", Field3D{1.0}}});
  auto state = makeState();
  state["species"]["e"]["density"] = Field3D{1.0};
  electrons.transform(state);
  component.transform(state);
  electrons.finally(state);
  auto output = diagnostics(component);
  // Curlb_B.z = 1 after normalisation: -Jdia_e dot Grad(phi) = -DDZ(phi).
  const Field3D expected = -DDZ(output["phi"].as<Field3D>());
  EXPECT_GT(max(abs(expected), false), 1e-5);
  EXPECT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["e"]["energy_source"]), expected,
                           "RGN_NOBNDRY", 1e-12));
  const Field3D pressure_advection = Div_n_bxGrad_f_B_XPPM(
      Field3D{1.0}, output["phi"].as<Field3D>(), false, false, true);
  EXPECT_TRUE(IsFieldEqual(
      Field3D{1.5 * (solver.getTimeDerivs()["Pe"].as<Field3D>() + pressure_advection)},
      expected, "RGN_NOBNDRY", 1e-12));
}

TEST_F(EvolveVorticityIonEnergyTest, PolarisationTransportUsesMassFractionAndTotalFlux) {
  static_cast<FakeMesh*>(mesh)->setGridDataSource(
      new FakeGridDataSource{{{"bxcvx", 0.0}, {"bxcvy", 0.0}, {"bxcvz", 0.5}}});
  FakeSolver solver;
  auto options = makeOptions();
  EvolveVorticityIonEnergy component("test", options, &solver);
  declareSpecies(component);
  // Use helium-3 so A/Z differs between species. rho=7/4 differs from
  // rho0=1.5. With v_d=2, v_he=-1, K_parallel=19/8 and total pressure=2.
  solver.setState({{"Vort", Field3D{0.0}}, {"Eion", Field3D{43. / 8}}});
  auto state = makeState();
  state["species"]["he++"]["AA"].force(3.0);
  state["species"]["d+"]["velocity"] = Field3D{2.0};
  state["species"]["d+"]["momentum"] = Field3D{2.0};
  state["species"]["he++"]["velocity"] = Field3D{-1.0};
  state["species"]["he++"]["momentum"] = Field3D{-0.75};
  state["species"]["e"]["pressure"] = makeField<Field3D>(
      [&](const Ind3D& i) { return 1.0 + 0.1 * sin(2 * PI * i.z() / nz); });
  component.transform(state);
  component.finally(state);
  const auto output = diagnostics(component);
  const Field3D dN = get<Field3D>(state["species"]["d+"]["density_source"]);
  const Field3D heN = get<Field3D>(state["species"]["he++"]["density_source"]);
  EXPECT_GT(max(abs(dN), false), 1e-5);
  // Mass fractions are 4/7 and 3/7: Jpd/Jphe=4/3, dN/heN=8/3.
  EXPECT_TRUE(IsFieldEqual(dN, Field3D{(8. / 3) * heN}, "RGN_NOBNDRY", 1e-12));
  EXPECT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["d+"]["momentum_source"]),
                           Field3D{4 * dN}, "RGN_NOBNDRY", 1e-12));
  EXPECT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["he++"]["momentum_source"]),
                           Field3D{-3 * heN}, "RGN_NOBNDRY", 1e-12));
  const Field3D chi = output["ddt(Phi)"].as<Field3D>();
  const Field3D divJpol = -bout::ion_energy::Div_a_Grad_perp_XZ(Field3D{1.5}, chi, false);
  EXPECT_TRUE(
      IsFieldEqual(Field3D{dN + 2 * heN}, Field3D{-divJpol}, "RGN_NOBNDRY", 1e-12));
  // p_d=4/3, p_he=2/3, K_d=2, K_he=3/8. The coefficient of Jpol is
  // (8/7)*(10/3+2) + (6/7)*(5/3+3/8) = 659/84, plus phi.
  const Field3D E_adv = output["phi"].as<Field3D>() + 659. / 84;
  const Field3D expected = -bout::ion_energy::Div_a_Grad_perp_XZ(1.5 * E_adv, chi, false);
  EXPECT_TRUE(
      IsFieldEqual(output["DivEpol"].as<Field3D>(), expected, "RGN_NOBNDRY", 1e-12));
  // Internal transport contributes to the flux, not an additional source.
  EXPECT_TRUE(IsFieldEqual(output["Eion_src"].as<Field3D>(), 0.0, "RGN_NOBNDRY", 1e-12));
}

TEST_F(EvolveVorticityIonEnergyTest, IntegrationByPartsIncludesRadialBoundaryWork) {
  auto* coord = mesh->getCoordinates();
  coord->setContravariantMetricTensor(
      ContravariantMetricTensor{1.0, 1.0, 1.0, 0.0, 0.3, 0.0});
  const Field3D f = makeField<Field3D>([&](const Ind3D& i) {
    return 0.2 * i.x() + sin(0.4 * i.x() + 2 * PI * i.z() / nz);
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
      const BoutReal flux = 1.5 * (f[r] - f[l] + 0.3 * transverse);
      const BoutReal outward = i.x() == mesh->xend ? 1.0 : -1.0;
      boundary_work += outward * 0.5 * (f[l] + f[r]) * flux;
    }
  }
  EXPECT_GT(abs(boundary_work), 1e-3);
  EXPECT_NEAR(fLf, -2 * energy + boundary_work, 1e-11);
}

TEST_F(EvolveVorticityIonEnergyTest, ClosedRadialFluxIncludesCrossMetric) {
  auto* coord = mesh->getCoordinates();
  coord->setContravariantMetricTensor(
      ContravariantMetricTensor{1.0, 1.0, 1.0, 0.0, 0.3, 0.0});
  const Field3D a = makeField<Field3D>(
      [&](const Ind3D& i) { return 1.0 + 0.1 * i.x() + 0.2 * cos(2 * PI * i.z() / nz); });
  // Zero X derivative, but nonzero g13 contribution to the normal flux.
  const Field3D f =
      makeField<Field3D>([&](const Ind3D& i) { return sin(2 * PI * i.z() / nz); });
  const Field3D divergence = bout::ion_energy::Div_a_Grad_perp_XZ(a, f, false);
  BoutReal integral = 0.0;
  for (const auto& i : f.getRegion("RGN_NOBNDRY")) {
    integral += divergence[i] * coord->J()[i] * coord->dx()[i] * coord->dz()[i];
  }
  EXPECT_NEAR(integral, 0.0, 1e-12);
}

TEST_F(EvolveVorticityIonEnergyTest, VaryingSpeciesFractionsPreservePolarisationCurrent) {
  static_cast<FakeMesh*>(mesh)->setGridDataSource(
      new FakeGridDataSource{{{"bxcvx", 0.0}, {"bxcvy", 0.0}, {"bxcvz", 0.5}}});
  FakeSolver solver;
  auto options = makeOptions();
  EvolveVorticityIonEnergy component("test", options, &solver);
  declareSpecies(component);
  const BoutReal kx = PI / (mesh->xend - mesh->xstart + 1);
  const BoutReal kz = 2 * PI / nz;
  const Field3D vort = makeField<Field3D>(
      [&](const Ind3D& i) { return 0.02 * sin(kx * i.x() + 3 * kz * i.z()); });
  solver.setState({{"Vort", vort}, {"Eion", Field3D{4.0}}});
  auto state = makeState();
  state["species"]["d+"]["density"].force(makeField<Field3D>(
      [&](const Ind3D& i) { return 0.5 + 0.2 * sin(2 * kx * i.x() + kz * i.z()); }));
  state["species"]["he++"]["density"].force(makeField<Field3D>(
      [&](const Ind3D& i) { return 0.25 + 0.1 * cos(kx * i.x() - 2 * kz * i.z()); }));
  state["species"]["he++"]["AA"].force(3.0);
  component.transform(state);
  component.finally(state);
  auto output = diagnostics(component);
  const Field3D divJpol = output["DivJpol"].as<Field3D>();
  const Field3D charge_source =
      get<Field3D>(state["species"]["d+"]["density_source"])
      + 2 * get<Field3D>(state["species"]["he++"]["density_source"]);
  EXPECT_GT(max(abs(divJpol), false), 1e-5);
  // The species currents must sum to the common current even when their
  // spatially varying fractions would activate different slope limiters.
  EXPECT_TRUE(IsFieldEqual(charge_source, Field3D{-divJpol}, "RGN_NOBNDRY", 1e-12));
}
