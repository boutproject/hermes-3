#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/evolve_momentum.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using EvolveMomentumTest = FakeMeshFixture;

TEST_F(EvolveMomentumTest, CreateComponent) {
  FakeSolver solver;

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"i", {{"AA", 2.0}, {"charge", 1.0}}}};
  EvolveMomentum component("i", options, &solver);

  Options state = solver.getState();

  ASSERT_TRUE(state.isSet("NVi"));
}

TEST_F(EvolveMomentumTest, Transform) {
  FakeSolver solver;

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"i", {{"AA", 2.0}, {"charge", 1.0}}}};
  EvolveMomentum component("i", options, &solver);

  Options state;
  // Exception if density is not set
  EXPECT_THROW(component.transform(state), BoutException);
  state = {{"species", {{"i", {{"density", 1.0}, {"AA", 1.0}}}}}};
  component.transform(state);

  Options& species = state["species"]["i"];
  ASSERT_TRUE(species.isSet("velocity"));
  ASSERT_TRUE(species.isSet("momentum"));
}

TEST_F(EvolveMomentumTest, Finally) {
  FakeSolver solver;

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"i", {{"AA", 2.0}, {"charge", 1.0}}}};
  EvolveMomentum component("i", options, &solver);

  const Options state = {{"fastest_wave", 0.0},
                         {"species",
                          {{"i",
                            {{"AA", 1.0},
                             {"density", 1.0},
                             {"momentum", 0.0},
                             {"velocity", 0.0},
                             {"momentum_source", 0.5}}}}}};
  component.finally(state);

  Options ddt = solver.getTimeDerivs();

  ASSERT_TRUE(ddt.isSet("NVi"));
  Field3D ddt_NVi = ddt["NVi"].as<Field3D>();

  BOUT_FOR_SERIAL(i, ddt_NVi.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(0.5, ddt_NVi[i]);
  }
}

using EvolveMomentumBudgetTest = FakeMeshFixture_tmpl<10, 7, 8, false, 2, 2>;

TEST_F(EvolveMomentumBudgetTest, ForceAndWorkChannelsCloseTheMomentumRate) {
  static_cast<FakeMesh*>(mesh)->ix_separatrix = -1;
  mesh->getCoordinates()->setNon_uniform(false);
  for (bool diagnose : {false, true}) {
    SCOPED_TRACE(diagnose);
    FakeSolver solver;
    Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                    {"i",
                     {{"AA", 2.0},
                      {"charge", 1.0},
                      {"diagnose", diagnose},
                      {"poloidal_flows", false}}}};
    EvolveMomentum component("i", options, &solver);
    component.declareAllSpecies(SpeciesInformation({}, {}, {"i"}, {}));
    solver.setState({{"NVi", Field3D{0.8}}});
    Options state{{"fastest_wave", Field3D{1.0}},
                  {"species",
                   {{"i",
                     {{"AA", 2.0},
                      {"charge", 1.0},
                      {"density", Field3D{1.0}},
                      {"temperature", Field3D{1.0}},
                      {"momentum_source", Field3D{0.5}}}}}}};
    component.transform(state);
    state["species"]["i"]["velocity"] = Field3D{0.4};
    state["species"]["i"]["momentum"] = Field3D{0.8};
    state["species"]["i"]["pressure"] =
        makeField<Field3D>([](const Ind3D& i) { return 1.0 + 0.1 * i.y(); });
    state["fields"]["phi"] =
        makeField<Field3D>([](const Ind3D& i) { return 0.2 * i.y(); });
    component.finally(state);
    Options output{
        {"Nnorm", 1.0}, {"Tnorm", 1.0}, {"Omega_ci", 1.0}, {"Cs0", 1.0}, {"rho_s0", 1.0}};
    component.outputVars(output);
    EXPECT_TRUE(IsFieldEqual(solver.getTimeDerivs()["NVi"].as<Field3D>(), 0.2,
                             "RGN_NOBNDRY", 1e-12));
    EXPECT_EQ(output.isSet("momentum_electric_i"), diagnose);
    if (!diagnose) {
      continue;
    }
    EXPECT_TRUE(IsFieldEqual(output["momentum_electric_i"].as<Field3D>(), -0.2,
                             "RGN_NOBNDRY", 1e-12));
    EXPECT_TRUE(IsFieldEqual(output["momentum_pressure_i"].as<Field3D>(), -0.1,
                             "RGN_NOBNDRY", 1e-12));
    EXPECT_TRUE(IsFieldEqual(output["momentum_work_electric_i"].as<Field3D>(), -0.08,
                             "RGN_NOBNDRY", 1e-12));
    EXPECT_TRUE(IsFieldEqual(output["momentum_work_pressure_i"].as<Field3D>(), -0.04,
                             "RGN_NOBNDRY", 1e-12));
    Field3D sum{0.0}, work{0.0};
    for (const auto* channel :
         {"exb", "electric", "parallel_advection", "pressure", "sources", "other"}) {
      sum += output[std::string{"momentum_"} + channel + "_i"].as<Field3D>();
      work += output[std::string{"momentum_work_"} + channel + "_i"].as<Field3D>();
    }
    EXPECT_TRUE(IsFieldEqual(sum, 0.2, "RGN_NOBNDRY", 1e-12));
    EXPECT_TRUE(
        IsFieldEqual(sum, output["ddt(NVi)"].as<Field3D>(), "RGN_NOBNDRY", 1e-12));
    EXPECT_TRUE(IsFieldEqual(work, 0.08, "RGN_NOBNDRY", 1e-12));
  }
}
