#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/evolve_pressure.hxx"

#include <bout/field_factory.hxx> // For generating functions

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using EvolvePressureTest = FakeMeshFixture;

TEST_F(EvolvePressureTest, CreateComponent) {
  FakeSolver solver;

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"i", {{"AA", 2.0}, {"charge", 1.0}}}};
  EvolvePressure component("i", options, &solver);

  Options state = solver.getState();

  ASSERT_TRUE(state.isSet("Pi"));
}

TEST_F(EvolvePressureTest, Transform) {
  FakeSolver solver;

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"i", {{"AA", 2.0}, {"charge", 1.0}}}};
  EvolvePressure component("i", options, &solver);

  Options state;
  // Exception if density is not set
  EXPECT_THROW(component.transform(state), BoutException);

  state = {{"species", {{"i", {{"density", 1.0}}}}}};
  component.transform(state);

  Options& species = state["species"]["i"];
  ASSERT_TRUE(species.isSet("pressure"));
  ASSERT_TRUE(species.isSet("temperature"));
}

TEST_F(EvolvePressureTest, Finally) {
  FakeSolver solver;

  Options options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"i", {{"AA", 2.0}, {"charge", 1.0}, {"thermal_conduction", false}}}};
  EvolvePressure component("i", options, &solver);

  const Options state = {{"species",
                          {{"i",
                            {{"density", 1.0},
                             {"pressure", 1.0},
                             {"temperature", 1.0},
                             {"pressure_source", 0.5}}}}}};
#if CHECKLEVEL >= 1
  // Throws exception due to pressure_source
  EXPECT_THROW(component.finally(state), BoutException);
#endif

  const Options state2 = {{"species",
                           {{"i",
                             {{"density", 1.0},
                              {"pressure", 1.0},
                              {"temperature", 1.0},
                              {"energy_source", 0.5}}}}}};
  component.finally(state2);

  Options ddt = solver.getTimeDerivs();

  ASSERT_TRUE(ddt.isSet("Pi"));
  Field3D ddt_Pi = ddt["Pi"].as<Field3D>();

  BOUT_FOR_SERIAL(i, ddt_Pi.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ((2. / 3) * 0.5, ddt_Pi[i]);
  }
}

namespace {
// wider in x since the ExB operator reaches two cells out
using EvolvePressureChannelsTest = FakeMeshFixture_tmpl<6, 7, 8, false, 2>;

// every term in ddt(P) switched on, so each channel has something in it
Options channelOptions(bool diagnose, bool p_div_v) {
  return Options{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                 {"i",
                  {{"AA", 2.0},
                   {"charge", 1.0},
                   {"thermal_conduction", false},
                   {"diagnose", diagnose},
                   {"p_div_v", p_div_v},
                   // the poloidal ExB part reads DDY(phi) in corner guard cells that
                   // are never set (0 * junk on this metric, but junk can be NaN)
                   {"poloidal_flows", false},
                   {"temperature_floor", 3.0}, // so low_T_diffuse_perp is not zero
                   {"low_n_diffuse_perp", true},
                   {"low_T_diffuse_perp", true},
                   {"low_p_diffuse_perp", true},
                   {"hyper_z", 0.1},
                   {"hyper_z_T", 0.2},
                   {"damp_p_nt", true},
                   {"numerical_viscous_heating", true}}}};
}

// P is not N * T on purpose, so damp_p_nt does something
Options channelState() {
  Options field_options;
  auto make = [&](const std::string& expr) {
    return FieldFactory::get()->create3D(expr, &field_options, bout::globals::mesh);
  };
  return Options{
      {"species",
       {{"i",
         {{"AA", 2.0},
          {"charge", 1.0},
          {"density", make("1 + 0.1*x + 0.2*y + 0.05*sin(z)")},
          {"pressure", make("2 + 0.3*x*y + 0.1*cos(z)")},
          {"temperature", make("1.5 + 0.2*y + 0.1*sin(2*z)")},
          {"velocity", make("0.3 + 0.1*y + 0.05*x*cos(z)")},
          {"energy_source", make("0.2 + 0.1*x*sin(z)")},
          {"low_n_coeff", make("0.01*(1 + y)")}}}}},
      {"fields", {{"phi", make("0.5*x*x + 0.3*sin(z) + 0.1*y*cos(z)")}}}};
}

Options channelOutputs() {
  return Options{{"Nnorm", 1.0}, {"Tnorm", 1.0}, {"Omega_ci", 1.0}, {"rho_s0", 1.0}};
}

// 3/2 ddt(P) against the sum of what's saved, relative to the biggest 3/2 ddt(P)
void checkPressureClosure(bool p_div_v) {
  FakeSolver solver;
  Options options = channelOptions(true, p_div_v);
  EvolvePressure component("i", options, &solver);

  component.finally(channelState());
  Options outputs = channelOutputs();
  component.outputVars(outputs);

  const Field3D lhs = 1.5 * solver.getTimeDerivs()["Pi"].as<Field3D>();
  const Field3D work = get<Field3D>(outputs[p_div_v ? "Ei_PdivV" : "Ei_VgradP"]);
  const Field3D rhs = 1.5 * get<Field3D>(outputs["SPi"])
                      + get<Field3D>(outputs["Ei_Div_VExB_32P"])
                      + get<Field3D>(outputs["Ei_Div_Vpar_32P"]) + work
                      + get<Field3D>(outputs["Ei_numerical_P"])
                      + get<Field3D>(outputs["Ei_nvh"]);

  const BoutReal scale = max(abs(lhs), false, "RGN_NOBNDRY");
  ASSERT_GT(scale, 0.0);
  EXPECT_TRUE(IsFieldEqual(lhs, rhs, "RGN_NOBNDRY", 1e-12 * scale));

  // and each channel actually has something in it
  for (const std::string ch : {"Ei_Div_VExB_32P", "Ei_Div_Vpar_32P", "Ei_numerical_P"}) {
    EXPECT_GT(max(abs(get<Field3D>(outputs[ch])), false, "RGN_NOBNDRY"), 0.0) << ch;
  }
}
} // namespace

TEST_F(EvolvePressureChannelsTest, ChannelsCloseVgradP) { checkPressureClosure(false); }

TEST_F(EvolvePressureChannelsTest, ChannelsClosePdivV) { checkPressureClosure(true); }

TEST_F(EvolvePressureChannelsTest, NumericalChannelIsEachTerm) {
  // switching one numerical term off should move E_numerical_P by what it took out of
  // 3/2 ddt(P), for every term in turn
  auto run = [](Options& options) {
    FakeSolver solver;
    EvolvePressure component("i", options, &solver);
    component.finally(channelState());
    Options outputs = channelOutputs();
    component.outputVars(outputs);
    return std::make_pair(Field3D{1.5 * solver.getTimeDerivs()["Pi"].as<Field3D>()},
                          get<Field3D>(outputs["Ei_numerical_P"]));
  };
  Options options_all = channelOptions(true, false);
  const auto [ddt_all, num_all] = run(options_all);
  const BoutReal scale = max(abs(num_all), false, "RGN_NOBNDRY");

  for (const std::string option : {"low_n_diffuse_perp", "low_T_diffuse_perp",
                                   "low_p_diffuse_perp", "hyper_z", "hyper_z_T",
                                   "damp_p_nt"}) {
    Options options = channelOptions(true, false);
    if (option.substr(0, 5) == "hyper") {
      options["i"][option].force(-1.0); // off when not positive
    } else {
      options["i"][option].force(false);
    }
    const auto [ddt_one, num_one] = run(options);
    const Field3D d_ddt = ddt_all - ddt_one;
    const Field3D d_num = num_all - num_one;
    EXPECT_GT(max(abs(d_num), false, "RGN_NOBNDRY"), 0.0) << option;
    EXPECT_TRUE(IsFieldEqual(d_num, d_ddt, "RGN_NOBNDRY", 1e-12 * scale)) << option;
  }
}

TEST_F(EvolvePressureChannelsTest, DiagnoseDoesNotChangeDdt) {
  for (const bool p_div_v : {false, true}) {
    FakeSolver solver_off;
    Options options_off = channelOptions(false, p_div_v);
    EvolvePressure component_off("i", options_off, &solver_off);
    component_off.finally(channelState());
    Options outputs_off = channelOutputs();
    component_off.outputVars(outputs_off);

    FakeSolver solver_on;
    Options options_on = channelOptions(true, p_div_v);
    EvolvePressure component_on("i", options_on, &solver_on);
    component_on.finally(channelState());

    // nothing new saved when it's off
    EXPECT_FALSE(outputs_off.isSet("Ei_Div_VExB_32P"));
    EXPECT_FALSE(outputs_off.isSet("Ei_Div_Vpar_32P"));
    EXPECT_FALSE(outputs_off.isSet("Ei_numerical_P"));

    // interior only, the operators leave the guard cells of ddt unset
    const Field3D off = solver_off.getTimeDerivs()["Pi"].as<Field3D>();
    const Field3D on = solver_on.getTimeDerivs()["Pi"].as<Field3D>();
    EXPECT_TRUE(IsFieldEqual(on, off, "RGN_NOBNDRY", 0.0)) << "p_div_v = " << p_div_v;
  }
}
