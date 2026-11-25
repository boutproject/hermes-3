#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include "test_extras.hxx"

#include "../../include/adaptive_index_diffusion.hxx"
#include "../../include/component_scheduler.hxx"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <type_traits>

using bout::globals::mesh;

class IndexDiffusionMesh : public FakeMesh {
public:
  using FakeMesh::FakeMesh;
  int x_calls{0}, y_calls{0}, all_calls{0};
  bool first_x{true}, last_x{true};
  bool firstX() const override { return first_x; }
  bool lastX() const override { return last_x; }
  comm_handle send(FieldGroup& group) override {
    ++all_calls;
    return FakeMesh::send(group);
  }
  comm_handle sendX(FieldGroup& group, comm_handle = nullptr, bool = false) override {
    ++x_calls;
    return FakeMesh::send(group);
  }
  comm_handle sendY(FieldGroup&, comm_handle = nullptr) override {
    ++y_calls;
    return nullptr;
  }
};

class IndexDiffusionTest
    : public FakeMeshFixture_tmpl<9, 7, 8, false, 2, 2, IndexDiffusionMesh> {};

TEST_F(IndexDiffusionTest, ConservesUnequalVolumesAndBoundsLocalRates) {
  const auto volumes = makeField<Field2D>(
      [](Ind2D& i) { return std::pow(10.0, 3.0 * (i.x() - 4)) * (1.0 + 0.1 * i.y()); });
  mesh->getCoordinates()->setJ(volumes);
  const auto rx = makeField<Field2D>([](Ind2D& i) { return 0.3 * (i.x() - 2); });
  const Field2D rz{0.7};
  Field3D f = makeField<Field3D>([](Ind3D& i) {
    return 2.0 + 0.1 * i.x() * i.x() + std::cos(2.0 * PI * i.z() / 8.0);
  });
  Field3D weight = makeField<Field3D>(
      [](Ind3D& i) { return 0.25 + 0.1 * (1.0 + std::cos(2.0 * PI * i.z() / 8.0)); });
  for (const auto& i : f.getRegion("RGN_GUARDS")) {
    f[i] = std::numeric_limits<BoutReal>::quiet_NaN();
    weight[i] = std::numeric_limits<BoutReal>::quiet_NaN();
  }
  getMesh().x_calls = 0;
  getMesh().y_calls = 0;
  getMesh().all_calls = 0;
  const auto faces = indexDiffusionFaces(rx, rz, weight);
  const auto rhs = applyIndexDiffusion(faces, f);
  const auto loss = indexDiffusionLossRate(faces);
  EXPECT_EQ(getMesh().x_calls, 0);
  EXPECT_EQ(getMesh().y_calls, 0);
  EXPECT_EQ(getMesh().all_calls, 0);
  BoutReal dissipation = 0.0, absolute_flux = 0.0;
  for (const auto& i : rhs.getRegion("RGN_NOBNDRY")) {
    const auto v = volumes[i];
    EXPECT_TRUE(std::isfinite(rhs[i]));
    EXPECT_GE(loss[i], 0.0);
    EXPECT_LE(loss[i], rx[i] + rz[i] + 1e-14);
    dissipation += v * f[i] * rhs[i];
    absolute_flux += v * std::abs(rhs[i]);
    if (i.x() == mesh->xstart) {
      EXPECT_EQ(faces.x[i.xm()], 0.0);
      EXPECT_EQ(faces.x[i], 0.0); // zero local ceiling closes its interior face
    }
    if (i.x() == mesh->xend) {
      EXPECT_EQ(faces.x[i], 0.0);
    }
  }
  EXPECT_NEAR(DomainBalance::integral(rhs), 0.0, 1e-14 * absolute_flux);
  EXPECT_LT(dissipation, 0.0);
  EXPECT_TRUE(IsFieldEqual(applyIndexDiffusion(faces, Field3D{4.0}), 0.0));
  for (const auto& i : rhs.getRegion("RGN_GUARDS")) {
    EXPECT_EQ(rhs[i], 0.0);
    EXPECT_TRUE(std::isnan(f[i]));
  }
}

TEST_F(IndexDiffusionTest, PeriodicXAndZFourierModesHaveIndexRates) {
  mesh->periodicX = true;
  getMesh().copy_periodic_halos = true;
  mesh->getCoordinates()->setDx(0.01);
  mesh->getCoordinates()->setDz(3.0);
  auto f = makeField<Field3D>([](Ind3D& i) {
    return std::cos(2.0 * PI * (i.x() - 2) / 5.0) + std::cos(2.0 * PI * i.z() / 8.0);
  });
  // The caller supplies valid periodic guards, just as equation owners do.
  mesh->communicateXZ(f);
  const auto faces = indexDiffusionFaces(Field2D{2.0}, Field2D{3.0}, Field3D{1.0});
  const auto rhs = applyIndexDiffusion(faces, f);
  for (const auto& i : rhs.getRegion("RGN_NOBNDRY")) {
    const auto expected =
        -2.0 * std::pow(std::sin(PI / 5.0), 2) * std::cos(2.0 * PI * (i.x() - 2) / 5.0)
        - 3.0 * std::pow(std::sin(PI / 8.0), 2) * std::cos(2.0 * PI * i.z() / 8.0);
    EXPECT_NEAR(rhs[i], expected, 1e-12);
  }
  EXPECT_NEAR(DomainBalance::integral(rhs), 0.0, 1e-12);
}

TEST_F(IndexDiffusionTest, ReadsExistingProcessorGuardsWithoutCommunication) {
  getMesh().first_x = false;
  getMesh().last_x = false;
  const auto nan = std::numeric_limits<BoutReal>::quiet_NaN();
  Field3D f{nan}, activation{nan};
  for (const auto& i : indexDiffusionInputRegion(f)) {
    activation[i] = 1.0;
    f[i] = i.x() < mesh->xstart ? 2.0 : (i.x() > mesh->xend ? 3.0 : 0.0);
  }
  getMesh().x_calls = getMesh().y_calls = getMesh().all_calls = 0;
  const auto faces = indexDiffusionFaces(Field2D{4.0}, Field2D{0.0}, activation);
  const auto rhs = applyIndexDiffusion(faces, f);
  for (const auto& i : rhs.getRegion("RGN_NOBNDRY")) {
    const auto expected = i.x() == mesh->xstart ? 2.0 : (i.x() == mesh->xend ? 3.0 : 0.0);
    EXPECT_DOUBLE_EQ(rhs[i], expected);
  }
  EXPECT_EQ(f(mesh->xstart - 1, mesh->ystart, 0), 2.0);
  EXPECT_EQ(f(mesh->xend + 1, mesh->ystart, 0), 3.0);
  EXPECT_TRUE(std::isnan(f(0, mesh->ystart, 0)));
  EXPECT_TRUE(std::isnan(f(mesh->xstart, 0, 0)));
  EXPECT_EQ(getMesh().x_calls, 0);
  EXPECT_EQ(getMesh().y_calls, 0);
  EXPECT_EQ(getMesh().all_calls, 0);
}

template <class Size>
class IndexDiffusionSmallZTest
    : public FakeMeshFixture_tmpl<7, 5, Size::value, false, 2, 1> {};
using ZSizes =
    testing::Types<std::integral_constant<int, 1>, std::integral_constant<int, 2>>;
TYPED_TEST_SUITE(IndexDiffusionSmallZTest, ZSizes);
TYPED_TEST(IndexDiffusionSmallZTest, PeriodicStencilAndMaximumPrinciple) {
  const auto f = makeField<Field3D>([](Ind3D& i) { return i.z() == 0 ? 1.0 : 0.0; });
  const auto faces = indexDiffusionFaces(Field2D{0.0}, Field2D{2.0}, Field3D{1.0});
  const auto rhs = applyIndexDiffusion(faces, f);
  for (const auto& i : rhs.getRegion("RGN_NOBNDRY")) {
    const BoutReal expected = TypeParam::value == 1 ? 0.0 : (i.z() == 0 ? -1.0 : 1.0);
    EXPECT_DOUBLE_EQ(rhs[i], expected);
    EXPECT_GE(f[i] + 0.5 * rhs[i], 0.0);
    EXPECT_LE(f[i] + 0.5 * rhs[i], 1.0);
  }
}

class AdaptiveIndexDiffusionTest : public IndexDiffusionTest {
protected:
  Options options() {
    return {{"units", {{"seconds", 0.5}, {"inv_meters_cubed", 10.0}, {"eV", 2.0}}},
            {"hermes", {{"components", "adaptive_index_diffusion"}}},
            {"adaptive_index_diffusion",
             {{"rate_x", 4.0},
              {"rate_z", 6.0},
              {"density_threshold", 20.0},
              {"diagnose", true}}}};
  }
  Options state() {
    return {{"species",
             {{"e", {{"charge", -1.0}, {"density", Field3D{2.0}}}},
              {"d+",
               {{"charge", 1.0},
                {"density", Field3D{2.0}},
                {"index_diffusion_density", Field3D{2.0}}}}}}};
  }
  Options diagnostics(AdaptiveIndexDiffusion& component) {
    Options output{{"Nnorm", 10.0}, {"Tnorm", 2.0}, {"Omega_ci", 2.0}, {"rho_s0", 1.0}};
    component.outputVars(output);
    return output;
  }
};

TEST_F(AdaptiveIndexDiffusionTest, DetectsEnabledAliasesWithoutCreatingMissingSections) {
  const Options absent;
  EXPECT_FALSE(usesAdaptiveIndexDiffusion(absent));
  const Options scheduled{{"hermes", {{"components", "e, sheath_boundary"}}}};
  EXPECT_FALSE(usesAdaptiveIndexDiffusion(scheduled));
  EXPECT_EQ(scheduled.getChildren().count("sheath_boundary"), 0);
  const Options aliased{
      {"components", "transport, e"},
      {"transport", {{"type", "adaptive_index_diffusion"}, {"rate_x", 1.0}}}};
  EXPECT_TRUE(usesAdaptiveIndexDiffusion(aliased));
  const Options omitted{{"components", "adaptive_index_diffusion, e"}};
  EXPECT_FALSE(usesAdaptiveIndexDiffusion(omitted));
  EXPECT_EQ(omitted.getChildren().count("adaptive_index_diffusion"), 0);
}

TEST_F(AdaptiveIndexDiffusionTest, ZeroRatesHaveNoDependenciesOrStateEffects) {
  auto opts = options();
  // Only physical cells determine whether the component is enabled.
  auto rx = Field2D{1.0};
  for (const auto& i : rx.getRegion("RGN_NOBNDRY")) {
    rx[i] = 0.0;
  }
  opts["adaptive_index_diffusion"]["rate_x"].force(rx);
  opts["adaptive_index_diffusion"]["rate_z"].force(0.0);
  opts["adaptive_index_diffusion"]["temperature_threshold"] = 1.0;
  EXPECT_FALSE(usesAdaptiveIndexDiffusion(opts));
  AdaptiveIndexDiffusion component("adaptive_index_diffusion", opts, nullptr);
  component.declareAllSpecies({"e", "d+"});
  EXPECT_TRUE(component.getPermissions()
                  .getVariablesWithMinimumPermission(PermissionTypes::ReadIfSet)
                  .empty());
  // No species, temperature, or geometry is needed by a disabled component.
  Options empty;
  component.transform(empty);
  EXPECT_TRUE(empty.getChildren().empty());
  EXPECT_EQ(diagnostics(component).getChildren().count("index_diffusion_activation"), 0);
}

TEST_F(AdaptiveIndexDiffusionTest, DerivedFieldsUseExistingProcessorGuards) {
  getMesh().first_x = false;
  getMesh().last_x = false;
  auto opts = options();
  opts["adaptive_index_diffusion"]["rate_z"].force(0.0);
  opts["adaptive_index_diffusion"]["temperature_threshold"] = 4.0;
  AdaptiveIndexDiffusion component("adaptive_index_diffusion", opts, nullptr);
  component.declareAllSpecies({"e", "d+"});
  auto s = state();
  const auto n = makeField<Field3D>([](Ind3D& i) {
    return i.x() < mesh->xstart ? 1.0 : (i.x() > mesh->xend ? 4.0 : 2.0);
  });
  const auto temperature = makeField<Field3D>([](Ind3D& i) {
    return i.x() < mesh->xstart ? 0.5 : (i.x() > mesh->xend ? 4.0 : 2.0);
  });
  for (const auto* name : {"e", "d+"}) {
    auto& sp = s["species"][name];
    const BoutReal mass = std::string{name} == "e" ? 0.1 : 2.0;
    sp["density"].force(n);
    if (std::string{name} != "e") {
      sp["index_diffusion_density"].force(n);
      sp["temperature"] = temperature;
    }
    sp["AA"] = mass;
    sp["index_diffusion_density_floor"] = 0.01;
    sp["index_diffusion_momentum"] = Field3D{0.02};
    sp["velocity"] = Field3D{0.02 / (mass * softFloor(n, 0.01))};
  }
  getMesh().x_calls = getMesh().y_calls = getMesh().all_calls = 0;
  component.transform(s);
  const auto activation =
      diagnostics(component)["index_diffusion_activation"].as<Field3D>();
  for (int x : {mesh->xstart - 1, mesh->xend + 1}) {
    const auto i = indexAt(n, x, mesh->ystart, 0);
    const auto wn = 4.0 / (4.0 + n[i] * n[i]);
    const auto wt = 4.0 / (4.0 + temperature[i] * temperature[i]);
    EXPECT_NEAR(activation[i], wn + (1.0 - wn) * wt, 1e-14);
    for (const auto* name : {"e", "d+"}) {
      const auto& sp = s["species"][name];
      const auto v = sp["velocity"].as<Field3D>()[i];
      const auto mass = sp["AA"].as<BoutReal>();
      EXPECT_NEAR(sp["index_diffusion_kinetic_energy"].as<Field3D>()[i],
                  0.5 * mass * n[i] * v * v, 1e-14);
    }
  }
  const auto sn = s["species"]["d+"]["index_diffusion_density_rhs"].as<Field3D>();
  EXPECT_LT(sn(mesh->xstart, mesh->ystart, 0), 0.0);
  EXPECT_GT(sn(mesh->xend, mesh->ystart, 0), 0.0);
  EXPECT_TRUE(IsFieldEqual(sn,
                           s["species"]["e"]["index_diffusion_density_rhs"].as<Field3D>(),
                           "RGN_NOBNDRY", 1e-14));
  EXPECT_EQ(getMesh().x_calls, 0);
  EXPECT_EQ(getMesh().y_calls, 0);
  EXPECT_EQ(getMesh().all_calls, 0);
}

TEST_F(AdaptiveIndexDiffusionTest, SIUnitsAndEitherTriggerActivateSharedFaces) {
  for (const auto& triggers :
       {std::string{"density"}, std::string{"temperature"}, std::string{"both"}}) {
    auto opts = options();
    if (triggers == "temperature") {
      // Reconstruct options to actually omit the density trigger.
      Options temperature_only{{"adaptive_index_diffusion",
                                {{"rate_x", 4.0}, {"rate_z", 6.0}, {"diagnose", true}}}};
      temperature_only["units"] = opts["units"].copy();
      temperature_only["hermes"] = opts["hermes"].copy();
      opts = std::move(temperature_only);
    }
    if (triggers != "density") {
      opts["adaptive_index_diffusion"]["temperature_threshold"] = 4.0;
    }
    AdaptiveIndexDiffusion component("adaptive_index_diffusion", opts, nullptr);
    component.declareAllSpecies({"e", "d+"});
    auto s = state();
    if (triggers != "density") {
      s["species"]["d+"]["temperature"] = Field3D{2.0};
    }
    component.transform(s);
    auto out = diagnostics(component);
    const auto activation = out["index_diffusion_activation"].as<Field3D>();
    const auto expected = triggers == "both" ? 0.75 : 0.5;
    EXPECT_TRUE(IsFieldEqual(activation, expected, "RGN_NOBNDRY"));
    const auto loss = out["index_diffusion_loss_rate"].as<Field3D>();
    EXPECT_NEAR(loss(4, 3, 0), 0.5 * (2.0 + 3.0) * expected, 1e-14);
    EXPECT_EQ(out["index_diffusion_loss_rate"].attributes["conversion"].as<BoutReal>(),
              2.0);
    EXPECT_FALSE(s["species"]["d+"].isSet("density_source"));
    EXPECT_FALSE(s["species"]["d+"].isSet("momentum_source"));
    EXPECT_FALSE(s["species"]["d+"].isSet("energy_source"));
  }
}

TEST_F(AdaptiveIndexDiffusionTest, ActivationChangesWithoutHistoryOrDivisionFloors) {
  auto opts = options();
  AdaptiveIndexDiffusion component("adaptive_index_diffusion", opts, nullptr);
  component.declareAllSpecies({"e", "d+"});
  for (const BoutReal n : {0.0, 1e200, -1.0, 2.0}) {
    auto s = state();
    s["species"]["d+"]["density"] = Field3D{n};
    component.transform(s);
    const auto expected = n <= 0.0 ? 1.0 : (n > 1e100 ? 0.0 : 0.5);
    EXPECT_TRUE(
        IsFieldEqual(diagnostics(component)["index_diffusion_activation"].as<Field3D>(),
                     expected, "RGN_NOBNDRY"));
  }
}

TEST_F(AdaptiveIndexDiffusionTest, ValidatesOptionsAndSpecies) {
  for (const auto* key :
       {"rate_x", "rate_z", "density_threshold", "temperature_threshold"}) {
    for (const auto value : {-1.0, std::numeric_limits<BoutReal>::infinity(),
                             std::numeric_limits<BoutReal>::quiet_NaN()}) {
      auto opts = options();
      opts["adaptive_index_diffusion"][key].force(value);
      EXPECT_THROW(AdaptiveIndexDiffusion("adaptive_index_diffusion", opts, nullptr),
                   BoutException);
    }
  }
  auto opts = options();
  opts["adaptive_index_diffusion"]["temperature_threshold"] = 4.0;
  opts["adaptive_index_diffusion"]["temperature_species"] = "missing+";
  AdaptiveIndexDiffusion component("adaptive_index_diffusion", opts, nullptr);
  component.declareAllSpecies({"e", "d+"});
  auto s = state();
  EXPECT_THROW(component.transform(s), BoutException);
  auto zero_threshold = options();
  zero_threshold["adaptive_index_diffusion"]["density_threshold"].force(0.0);
  EXPECT_THROW(
      AdaptiveIndexDiffusion("adaptive_index_diffusion", zero_threshold, nullptr),
      BoutException);
  auto duplicate = options();
  duplicate["hermes"]["components"].force("adaptive_index_diffusion, other");
  duplicate["other"]["type"] = "adaptive_index_diffusion";
  EXPECT_THROW(AdaptiveIndexDiffusion("adaptive_index_diffusion", duplicate, nullptr),
               BoutException);
  auto conflict = options();
  conflict["hermes"]["components"].force("adaptive_index_diffusion, d+");
  conflict["d+"]["charge"] = 1.0;
  conflict["d+"]["low_n_diffuse_perp"] = true;
  EXPECT_THROW(AdaptiveIndexDiffusion("adaptive_index_diffusion", conflict, nullptr),
               BoutException);
}

TEST_F(AdaptiveIndexDiffusionTest, KineticRateDifferentiatesSoftFloorAtLowDensity) {
  auto opts = options();
  opts["adaptive_index_diffusion"]["density_threshold"].force(1.0);
  AdaptiveIndexDiffusion component("adaptive_index_diffusion", opts, nullptr);
  component.declareAllSpecies({"e", "d+"});
  auto s = state();
  auto& ion = s["species"]["d+"];
  const BoutReal nf = 0.01, mass = 2.0;
  const auto n = makeField<Field3D>(
      [](Ind3D& i) { return 0.004 + 0.002 * std::cos(2.0 * PI * i.z() / 8.0); });
  const auto m = makeField<Field3D>(
      [](Ind3D& i) { return 0.02 + 0.01 * std::sin(2.0 * PI * i.z() / 8.0); });
  ion["density"] = n;
  ion["index_diffusion_density"] = n;
  ion["index_diffusion_momentum"] = m;
  ion["index_diffusion_density_floor"] = nf;
  ion["AA"] = mass;
  ion["velocity"] = Field3D{m / (mass * softFloor(n, nf))};
  component.transform(s);
  const auto sn = ion["index_diffusion_density_rhs"].as<Field3D>();
  const auto sm = ion["index_diffusion_momentum_rhs"].as<Field3D>();
  const auto sk = ion["index_diffusion_kinetic_rhs"].as<Field3D>();
  const auto kinetic = [&](BoutReal nvalue, BoutReal mvalue) {
    const auto v = mvalue / (mass * softFloor(nvalue, nf));
    return 0.5 * mass * nvalue * v * v;
  };
  const BoutReal eps = 1e-5;
  for (const auto& i : n.getRegion("RGN_NOBNDRY")) {
    const auto difference = (kinetic(n[i] + eps * sn[i], m[i] + eps * sm[i])
                             - kinetic(n[i] - eps * sn[i], m[i] - eps * sm[i]))
                            / (2.0 * eps);
    EXPECT_NEAR(sk[i], difference, 1e-10);
  }
  EXPECT_GT(max(abs(sk), false, "RGN_NOBNDRY"), 0.0);
}

TEST_F(AdaptiveIndexDiffusionTest, QuasineutralKineticRateIncludesIndividualIonClipping) {
  auto opts = options();
  AdaptiveIndexDiffusion component("adaptive_index_diffusion", opts, nullptr);
  component.declareAllSpecies({"e", "d+", "he++"});
  auto s = state();
  const auto trace_raw = makeField<Field3D>(
      [](Ind3D& i) { return -0.05 + 0.08 * std::cos(2.0 * PI * i.z() / 8.0); });
  auto& trace = s["species"]["he++"];
  trace["charge"] = 2.0;
  trace["density"] = floor(trace_raw, 0.0);
  trace["index_diffusion_density"] = trace_raw;
  auto& e = s["species"]["e"];
  const Field3D ne = 2.0 + 2.0 * floor(trace_raw, 0.0);
  const Field3D me{0.02};
  e["density"] = ne;
  e["AA"] = 0.1;
  e["index_diffusion_momentum"] = me;
  e["index_diffusion_density_floor"] = 0.01;
  e["velocity"] = Field3D{me / (0.1 * softFloor(ne, 0.01))};
  component.transform(s);
  const auto sntrace = trace["index_diffusion_density_rhs"].as<Field3D>();
  const auto sme = e["index_diffusion_momentum_rhs"].as<Field3D>();
  const auto ske = e["index_diffusion_kinetic_rhs"].as<Field3D>();
  const auto kinetic = [&](BoutReal trace_value, BoutReal momentum) {
    const auto n = 2.0 + 2.0 * std::max(trace_value, 0.0);
    const auto v = momentum / (0.1 * softFloor(n, 0.01));
    return 0.05 * n * v * v;
  };
  const BoutReal eps = 1e-5;
  for (const auto& i : ske.getRegion("RGN_NOBNDRY")) {
    const auto expected =
        (kinetic(trace_raw[i] + eps * sntrace[i], me[i] + eps * sme[i])
         - kinetic(trace_raw[i] - eps * sntrace[i], me[i] - eps * sme[i]))
        / (2.0 * eps);
    EXPECT_NEAR(ske[i], expected, 1e-12);
  }
}

class AdaptiveIndexDiffusionEquationsTest : public IndexDiffusionTest {
protected:
  struct Evaluation {
    Options state, rhs;
  };
  Field3D density(bool low) {
    return makeField<Field3D>([=](Ind3D& i) {
      return (low ? 0.002 : 2.0)
             * (1.0 + 0.1 * i.x() + 0.2 * std::cos(2.0 * PI * i.z() / 8.0));
    });
  }
  Field3D momentum() {
    return makeField<Field3D>(
        [](Ind3D& i) { return 0.02 + 0.01 * std::sin(2.0 * PI * i.z() / 8.0); });
  }
  Field3D pressure() {
    return makeField<Field3D>([](Ind3D& i) {
      return 3.0 + 0.1 * i.x() + 0.2 * std::cos(2.0 * PI * i.z() / 8.0);
    });
  }
  Evaluation evaluate(bool enabled, bool low = false, bool energy = false,
                      bool logarithmic = false, BoutReal rate = 4.0) {
    Options opts{{"units", {{"seconds", 1.0}, {"inv_meters_cubed", 1.0}, {"eV", 1.0}}},
                 {"hermes",
                  {{"components", enabled ? "adaptive_index_diffusion, e, d+" : "e, d+"},
                   {"restarting", true}}},
                 {"adaptive_index_diffusion",
                  {{"rate_x", rate}, {"rate_z", rate}, {"density_threshold", 2.0}}},
                 {"e",
                  {{"type", "quasineutral, evolve_pressure, evolve_momentum"},
                   {"charge", -1.0},
                   {"AA", 0.1},
                   {"density_floor", 0.01},
                   {"numerical_viscous_heating", false}}},
                 {"d+",
                  {{"type", energy ? "evolve_density, evolve_momentum, evolve_energy"
                                   : "evolve_density, evolve_momentum, evolve_pressure"},
                   {"charge", 1.0},
                   {"AA", 2.0},
                   {"density_floor", 0.01},
                   {"numerical_viscous_heating", false},
                   {"evolve_log", logarithmic}}}};
    FakeSolver solver;
    auto scheduler = ComponentScheduler::create(opts["hermes"], opts, &solver);
    const auto n = density(low), m = momentum(), p = pressure();
    Options variables{{"Nd+", n},
                      {"NVd+", m},
                      {"Pd+", p},
                      {"NVe", Field3D{0.2 * m}},
                      {"Pe", Field3D{0.3 * p}},
                      {"Ed+", Field3D{1.5 * p}},
                      {"logNd+", log(n)},
                      {"logPd+", log(p)}};
    solver.setState(variables);
    Options state{{"fastest_wave", Field3D{0.0}}, {"time", 0.0}};
    scheduler->transform(state);
    return {std::move(state), solver.getTimeDerivs()};
  }
};

TEST_F(AdaptiveIndexDiffusionEquationsTest, SharedMomentsAndPressureConserveTotalEnergy) {
  for (const bool low : {false, true}) {
    const auto off = evaluate(false, low), on = evaluate(true, low);
    const IndexDiffusionFaces faces{
        on.state["fields"]["index_diffusion_conductance_x"].as<Field3D>(),
        on.state["fields"]["index_diffusion_conductance_z"].as<Field3D>()};
    const Field3D sn = on.rhs["Nd+"].as<Field3D>() - off.rhs["Nd+"].as<Field3D>();
    const Field3D sm = on.rhs["NVd+"].as<Field3D>() - off.rhs["NVd+"].as<Field3D>();
    EXPECT_TRUE(
        IsFieldEqual(sn, applyIndexDiffusion(faces, density(low)), "RGN_NOBNDRY", 1e-12));
    EXPECT_TRUE(
        IsFieldEqual(sm, applyIndexDiffusion(faces, momentum()), "RGN_NOBNDRY", 1e-12));
    EXPECT_NEAR(DomainBalance::integral(sn), 0.0, 1e-12);
    EXPECT_NEAR(DomainBalance::integral(sm), 0.0, 1e-12);
    const auto inherited =
        on.state["species"]["e"]["index_diffusion_density_rhs"].as<Field3D>();
    EXPECT_TRUE(IsFieldEqual(sn, inherited, "RGN_NOBNDRY", 1e-12));
    for (const auto& name : {std::string{"d+"}, std::string{"e"}}) {
      const auto& sp = on.state["species"][name];
      const Field3D dp =
          on.rhs["P" + name].as<Field3D>() - off.rhs["P" + name].as<Field3D>();
      const auto kinetic = sp["index_diffusion_kinetic_energy"].as<Field3D>();
      const auto sk = sp["index_diffusion_kinetic_rhs"].as<Field3D>();
      const Field3D total = 1.5 * dp + sk;
      const Field3D initial_pressure =
          name == "d+" ? pressure() : Field3D{0.3 * pressure()};
      EXPECT_TRUE(IsFieldEqual(
          total, applyIndexDiffusion(faces, Field3D{1.5 * initial_pressure + kinetic}),
          "RGN_NOBNDRY", 1e-11));
      EXPECT_NEAR(DomainBalance::integral(total), 0.0, 1e-11);
      // Numerical mixing is separate from physical source conversion.
      if (sp.isSet("density_source")) {
        EXPECT_EQ(max(abs(sp["density_source"].as<Field3D>()), false), 0.0);
      }
      EXPECT_FALSE(sp.isSet("momentum_source"));
      EXPECT_FALSE(sp.isSet("energy_source"));
    }
  }
}

TEST_F(AdaptiveIndexDiffusionEquationsTest,
       ZeroRatesPreserveBaselineAndLogarithmicRates) {
  for (const bool logs : {false, true}) {
    const auto off = evaluate(false, false, false, logs);
    const auto zero = evaluate(true, false, false, logs, 0.0);
    EXPECT_FALSE(zero.state.isSection("fields"));
    for (const auto* name : {"d+", "e"}) {
      const auto& sp = zero.state["species"][name];
      EXPECT_FALSE(sp.isSet("index_diffusion_density"));
      EXPECT_FALSE(sp.isSet("index_diffusion_momentum"));
      EXPECT_FALSE(sp.isSet("index_diffusion_kinetic_rhs"));
    }
    for (const auto& [name, value] : off.rhs.getChildren()) {
      const auto expected = value.as<Field3D>();
      const auto actual = zero.rhs[name].as<Field3D>();
      for (const auto& i : expected.getRegion("RGN_NOBNDRY")) {
        EXPECT_DOUBLE_EQ(actual[i], expected[i]) << name << " at " << i;
      }
    }
    const auto on = evaluate(true, false, false, logs);
    const auto& ion = on.state["species"]["d+"];
    if (logs) {
      const Field3D dn =
          (on.rhs["logNd+"].as<Field3D>() - off.rhs["logNd+"].as<Field3D>())
          * density(false);
      EXPECT_TRUE(IsFieldEqual(dn, ion["index_diffusion_density_rhs"].as<Field3D>(),
                               "RGN_NOBNDRY", 1e-12));
    }
  }
}

TEST_F(AdaptiveIndexDiffusionEquationsTest, AlternativeTotalEnergyUsesDirectDiffusion) {
  const auto off = evaluate(false, false, true), on = evaluate(true, false, true);
  const IndexDiffusionFaces faces{
      on.state["fields"]["index_diffusion_conductance_x"].as<Field3D>(),
      on.state["fields"]["index_diffusion_conductance_z"].as<Field3D>()};
  const Field3D de = on.rhs["Ed+"].as<Field3D>() - off.rhs["Ed+"].as<Field3D>();
  EXPECT_TRUE(IsFieldEqual(de, applyIndexDiffusion(faces, Field3D{1.5 * pressure()}),
                           "RGN_NOBNDRY", 1e-12));
  EXPECT_NEAR(DomainBalance::integral(de), 0.0, 1e-12);
}
