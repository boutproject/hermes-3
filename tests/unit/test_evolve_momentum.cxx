#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/div_ops.hxx"
#include "../../include/evolve_momentum.hxx"
#include "hermes_build_config.hxx" // hermes::Limiter, generated in the build tree

#include <bout/field_factory.hxx> // For generating functions
#include <bout/fv_ops.hxx>

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

TEST_F(EvolveMomentumTest, FinallyScalesMomentumTimeDerivative) {
  FakeSolver solver;

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"i", {{"AA", 2.0}, {"charge", 1.0}, {"scale_momentum_ddt", 2.0}}}};
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
    ASSERT_DOUBLE_EQ(1.0, ddt_NVi[i]);
  }
}

namespace {
// wider in x since the ExB operator reaches two cells out
using EvolveMomentumChannelsTest = FakeMeshFixture_tmpl<6, 7, 8, false, 2>;

// every term in ddt(NV) switched on, so each channel has something in it
Options channelOptions(bool diagnose) {
  return Options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                 {"i",
                  {{"AA", 2.0},
                   {"charge", 1.0},
                   {"diagnose", diagnose},
                   // the poloidal ExB part reads DDY(phi) in corner guard cells that
                   // are never set (0 * junk on this metric, but junk can be NaN)
                   {"poloidal_flows", false},
                   {"low_n_diffuse_perp", true},
                   {"low_p_diffuse_perp", true},
                   {"hyper_z", 0.1}}}};
}

Options channelState(bool low_n_coeff = true) {
  Options field_options;
  auto make = [&](const std::string& expr) {
    return FieldFactory::get()->create3D(expr, &field_options, bout::globals::mesh);
  };
  Options state{
      {"species",
       {{"i",
         {{"AA", 2.0},
          {"charge", 1.0},
          {"density", make("1 + 0.1*x + 0.2*y + 0.05*sin(z)")},
          {"momentum", make("0.6 + 0.2*y + 0.1*x*cos(z)")},
          {"velocity", make("0.3 + 0.1*y + 0.05*x*cos(z)")},
          {"pressure", make("2 + 0.3*x*y + 0.1*cos(z)")},
          {"temperature", make("1.5 + 0.2*y + 0.1*sin(2*z)")},
          {"momentum_source", make("0.2 + 0.1*x*sin(z)")}}}}},
      {"fields", {{"phi", make("0.5*x*x + 0.3*sin(z) + 0.1*y*cos(z)")}}}};
  if (low_n_coeff) {
    state["species"]["i"]["low_n_coeff"] = make("0.01*(1 + y)");
  }
  return state;
}

Options channelOutputs() {
  return Options{{"Nnorm", 1.0}, {"Tnorm", 1.0}, {"Omega_ci", 1.0},
                 {"Cs0", 1.0},   {"rho_s0", 1.0}};
}

const std::vector<std::string> v_channels{"Ei_Div_VExB_KE", "Ei_Div_Vpar_KE",
                                          "Ei_V_Gradpar_P", "Ei_j_Gradpar_phi",
                                          "Ei_numerical_V"};
} // namespace

TEST_F(EvolveMomentumChannelsTest, ChannelsCloseMomentum) {
  FakeSolver solver;
  Options options = channelOptions(true);
  EvolveMomentum component("i", options, &solver);
  const Options state = channelState();
  component.finally(state);
  Options outputs = channelOutputs();
  component.outputVars(outputs);

  // V ddt(NV) against the sum of the V weighted channels and V SNV
  const Field3D V = get<Field3D>(state["species"]["i"]["velocity"]);
  const Field3D lhs = V * solver.getTimeDerivs()["NVi"].as<Field3D>();
  Field3D rhs = V * get<Field3D>(outputs["SNVi"]);
  for (const auto& ch : v_channels) {
    const Field3D E = get<Field3D>(outputs[ch]);
    EXPECT_GT(max(abs(E), false, "RGN_NOBNDRY"), 0.0) << ch;
    EXPECT_EQ("i", outputs[ch].attributes["species"].as<std::string>());
    rhs += E;
  }
  const BoutReal scale = max(abs(lhs), false, "RGN_NOBNDRY");
  ASSERT_GT(scale, 0.0);
  EXPECT_TRUE(IsFieldEqual(lhs, rhs, "RGN_NOBNDRY", 1e-12 * scale));
}

TEST_F(EvolveMomentumChannelsTest, NumericalChannelIsEachTerm) {
  // switching one numerical term off should move E_numerical_V by what it took out of
  // V ddt(NV), for every term in turn
  auto run = [](Options& options, bool low_n_coeff) {
    FakeSolver solver;
    EvolveMomentum component("i", options, &solver);
    const Options state = channelState(low_n_coeff);
    component.finally(state);
    Options outputs = channelOutputs();
    component.outputVars(outputs);
    const Field3D V = get<Field3D>(state["species"]["i"]["velocity"]);
    return std::make_pair(Field3D{V * solver.getTimeDerivs()["NVi"].as<Field3D>()},
                          get<Field3D>(outputs["Ei_numerical_V"]));
  };
  Options options_all = channelOptions(true);
  const auto [ddt_all, num_all] = run(options_all, true);
  const BoutReal scale = max(abs(num_all), false, "RGN_NOBNDRY");

  for (const std::string option :
       {"low_n_diffuse_perp", "low_p_diffuse_perp", "hyper_z", "low_n_coeff"}) {
    Options options = channelOptions(true);
    if (option == "hyper_z") {
      options["i"][option].force(-1.0); // off when not positive
    } else if (option != "low_n_coeff") {
      options["i"][option].force(false);
    }
    const auto [ddt_one, num_one] = run(options, option != "low_n_coeff");
    const Field3D d_ddt = ddt_all - ddt_one;
    const Field3D d_num = num_all - num_one;
    EXPECT_GT(max(abs(d_num), false, "RGN_NOBNDRY"), 0.0) << option;
    EXPECT_TRUE(IsFieldEqual(d_num, d_ddt, "RGN_NOBNDRY", 1e-12 * scale)) << option;
  }
}

TEST_F(EvolveMomentumChannelsTest, FluxFormsAreTheOldDefinitions) {
  FakeSolver solver;
  Options options = channelOptions(true);
  EvolveMomentum component("i", options, &solver);
  const Options state = channelState();
  component.finally(state);
  Options outputs = channelOutputs();
  component.outputVars(outputs);

  // the old branch's operators on 1/2 NV V, signed so positive is a source
  const Options& species = state["species"]["i"];
  const Field3D NV = get<Field3D>(species["momentum"]);
  const Field3D V = get<Field3D>(species["velocity"]);
  const Field3D phi = get<Field3D>(state["fields"]["phi"]);
  const Field3D fastest_wave = sqrt(get<Field3D>(species["temperature"]) / 2.0);
  Field3D flow_ylow;
  const Field3D exb = -Div_n_bxGrad_f_B_XPPM(0.5 * NV * V, phi, true, false, true);
  const Field3D par =
      -FV::Div_par_mod<hermes::Limiter>(0.5 * NV * V, V, fastest_wave, flow_ylow);

  const Field3D E_exb = get<Field3D>(outputs["Ei_Div_VExB_KE_flux"]);
  const Field3D E_par = get<Field3D>(outputs["Ei_Div_Vpar_KE_flux"]);
  EXPECT_GT(max(abs(E_exb), false, "RGN_NOBNDRY"), 0.0);
  EXPECT_GT(max(abs(E_par), false, "RGN_NOBNDRY"), 0.0);
  EXPECT_TRUE(IsFieldEqual(E_exb, exb, "RGN_NOBNDRY", 0.0));
  EXPECT_TRUE(IsFieldEqual(E_par, par, "RGN_NOBNDRY", 0.0));
}

TEST_F(EvolveMomentumChannelsTest, DiagnoseDoesNotChangeDdt) {
  FakeSolver solver_off;
  Options options_off = channelOptions(false);
  EvolveMomentum component_off("i", options_off, &solver_off);
  component_off.finally(channelState());
  Options outputs_off = channelOutputs();
  component_off.outputVars(outputs_off);

  FakeSolver solver_on;
  Options options_on = channelOptions(true);
  EvolveMomentum component_on("i", options_on, &solver_on);
  component_on.finally(channelState());

  // nothing new saved when it's off
  for (const std::string ch : {"Ei_Div_VExB_KE", "Ei_Div_Vpar_KE", "Ei_V_Gradpar_P",
                               "Ei_j_Gradpar_phi", "Ei_numerical_V",
                               "Ei_Div_VExB_KE_flux", "Ei_Div_Vpar_KE_flux"}) {
    EXPECT_FALSE(outputs_off.isSet(ch)) << ch;
  }

  // interior only, the operators leave the guard cells of ddt unset
  const Field3D off = solver_off.getTimeDerivs()["NVi"].as<Field3D>();
  const Field3D on = solver_on.getTimeDerivs()["NVi"].as<Field3D>();
  EXPECT_TRUE(IsFieldEqual(on, off, "RGN_NOBNDRY", 0.0));
}
