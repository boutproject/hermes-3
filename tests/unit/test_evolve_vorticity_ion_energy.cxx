#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include "test_extras.hxx"

#include "../../include/evolve_vorticity_ion_energy.hxx"

#include <bout/constants.hxx>
#include <bout/difops.hxx>

using bout::globals::mesh;

class EvolveVorticityIonEnergyTest : public FakeMeshFixture_tmpl<10, 5, 8, false, 2> {
protected:
  // Delp2 uses a global Laplacian that must not outlive the fixture's mesh.
  void SetUp() override { Laplacian::cleanup(); }
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
              {"phi_max_iterations", 30}}}};
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

  // Away from radial boundaries, apply the inversion's uniform X finite
  // difference and full spectral Z Laplacian, including any DC boundary value.
  for (int y = mesh->ystart; y <= mesh->yend; ++y) {
    for (int z = 0; z < nz; ++z) {
      BoutReal d2z = 0.0;
      for (int mode = 0; mode < nz; ++mode) {
        const BoutReal wave_number = 2 * PI * (mode <= nz / 2 ? mode : mode - nz) / nz;
        for (int neighbour = 0; neighbour < nz; ++neighbour) {
          d2z -= SQ(wave_number) * cos(wave_number * (z - neighbour))
                 * chi(4, y, neighbour) / nz;
        }
      }
      const BoutReal laplacian = chi(5, y, z) - 2 * chi(4, y, z) + chi(3, y, z) + d2z;
      EXPECT_NEAR(1.5 * laplacian, rhs(4, y, z), 1e-12);
    }
  }
}

TEST_F(EvolveVorticityIonEnergyTest, PolarisationFluxConservesEnergyOnNonuniformMesh) {
  mesh->addBoundary(
      bout::boundary::NewBoundaryRegionXIn("core", mesh->ystart, mesh->yend, mesh));
  mesh->addBoundary(
      bout::boundary::NewBoundaryRegionXOut("sol", mesh->ystart, mesh->yend, mesh));
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
