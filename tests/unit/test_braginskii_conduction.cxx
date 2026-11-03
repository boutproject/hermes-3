#include <string>
#include <utility>

#include <bout/bout_types.hxx>
#include <bout/region.hxx>
#include <bout/solver.hxx>
#include <bout/utils.hxx>
#include "gtest/gtest.h"

#include "../../include/braginskii_conduction.hxx"
#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
class BraginskiiConductionTest : public FakeMeshFixture {
public:
  BraginskiiConductionTest()
      : FakeMeshFixture(), options({
                               {"units",
                                {{"eV", 1.0},
                                 {"meters", 1.0},
                                 {"seconds", 1.0},
                                 {"inv_meters_cubed", 1e19}}},
                               {"test+",
                                {{"type", "evolve_energy"},
                                 {"thermal_conduction", true},
                                 {"diagnose", true}}},
                           }),
        output_base({{"Nnorm", 1e19}, {"Tnorm", 1.}, {"Omega_ci", 1.}, {"rho_s0", 1.}}),
        state_base(
            {{"fastest_wave", 0.},
             {"species",
              {{"test+",
                {{"velocity", 0.}, {"density", 1.}, {"pressure", 1.}, {"AA", 2.}}}}}}) {
    Options::root()["Ptest+"]["function"] = "1.0";
    temp1 = linearGradient(1., 0., 0., 1., 1., 0., 0.);
    temp2 = linearGradient(0., 0., 0., 1., 2., 0., 0.);
    temp_perp_variation = linearGradient(0., 1., 2., 0., 0., 0., 0.);
  }
  Options options, output_base, state_base;
  Field3D temp1, temp2, temp_perp_variation;

  static Field3D linearGradient(BoutReal a, BoutReal b1_x, BoutReal b2_x, BoutReal b1_y,
                                BoutReal b2_y, BoutReal b1_z, BoutReal b2_z) {
    Field3D result(a);
    BOUT_FOR_SERIAL(i, result.getRegion("RGN_ALL")) {
      const int x = i.x();
      const int y = i.y();
      const int z = i.z();
      result[i] +=
          b1_x * x + b2_x * SQ(x) + b1_y * y + b2_y * SQ(y) + b1_z * z + b2_z * SQ(z);
    }
    return result;
  }

  BraginskiiConduction makeComponent(BoutReal kappa, std::string mode = "braginskii") {
    Options opts = options.copy();
    opts["test+"]["kappa_coefficient"] = kappa;
    opts["test+"]["conduction_collisions_mode"] = std::move(mode);
    return BraginskiiConduction("test+", opts, nullptr);
  }

  Field3D getDeriv(const Options& state) {
    return state["species"]["test+"]["energy_source"];
  }
};

TEST_F(BraginskiiConductionTest, ConductionGradientScaling) {
  BraginskiiConduction component = this->makeComponent(1.);
  Options state0;
  Options state1 = this->state_base.copy();
  Options state2;
  state1["species"]["test+"]["charge"] = 1;
  state1["species"]["test+"]["collision_frequencies"]["test+_test+_coll"] = 1.;
  state1["species"]["test+"]["collision_frequency"] = 1.;

  state0 = state1.copy();
  state2 = state1.copy();

  state0["species"]["test+"]["temperature"] = this->temp_perp_variation;
  state1["species"]["test+"]["temperature"] = this->temp1;
  state2["species"]["test+"]["temperature"] = this->temp2;

  component.declareAllSpecies({"test+"});
  component.transform(state0);
  component.transform(state1);
  component.transform(state2);

  Field3D conduction0 = this->getDeriv(state0);
  Field3D conduction1 = this->getDeriv(state1);
  Field3D conduction2 = this->getDeriv(state2);
  BOUT_FOR_SERIAL(i, conduction1.getRegion("RGN_NOBNDRY")) {
    // Conduction is proportional to the parallel second derivative of temperature
    EXPECT_NE(conduction1[i], 0.);
    EXPECT_DOUBLE_EQ(conduction0[i], 0.);
    EXPECT_DOUBLE_EQ(conduction2[i], 2 * conduction1[i]);
  }
}

TEST_F(BraginskiiConductionTest, ConductionKappaScaling) {
  BraginskiiConduction component1 = this->makeComponent(1.);
  BraginskiiConduction component2 = this->makeComponent(2.);
  Options state1 = this->state_base.copy();
  state1["species"]["test+"]["charge"] = 1;
  state1["species"]["test+"]["collision_frequencies"]["test+_test+_coll"] = 1.;
  state1["species"]["test+"]["collision_frequency"] = 1.;
  state1["species"]["test+"]["temperature"] = this->temp1;
  Options state2 = state1.copy();

  component1.declareAllSpecies({"test+"});
  component1.transform(state1);
  component2.declareAllSpecies({"test+"});
  component2.transform(state2);

  Field3D conduction1 = this->getDeriv(state1);
  Field3D conduction2 = this->getDeriv(state2);
  BOUT_FOR_SERIAL(i, conduction1.getRegion("RGN_NOBNDRY")) {
    // Conduction is proportional to the parallel diffusivity
    EXPECT_NE(conduction1[i], 0.);
    EXPECT_DOUBLE_EQ(conduction2[i], 2 * conduction1[i]);
  }
}

TEST_F(BraginskiiConductionTest, ConductionCollisionScaling) {
  BraginskiiConduction component = this->makeComponent(1.);
  Options state0;
  Options state1 = this->state_base.copy();
  Options state2;
  state1["species"]["test+"]["charge"] = 1;
  state1["species"]["test+"]["temperature"] = this->temp1;

  state0 = state1.copy();
  state2 = state1.copy();

  state0["species"]["test+"]["collision_frequencies"]["test+_test+_coll"] = 0.5;
  state0["species"]["test+"]["collision_frequency"] = 0.5;
  state1["species"]["test+"]["collision_frequencies"]["test+_test+_coll"] = 1.;
  state1["species"]["test+"]["collision_frequency"] = 1.;
  state2["species"]["test+"]["collision_frequencies"]["test+_test+_coll"] = 2.;
  state2["species"]["test+"]["collision_frequency"] = 2.;

  component.declareAllSpecies({"test+"});
  component.transform(state0);
  component.transform(state1);
  component.transform(state2);

  Field3D conduction0 = this->getDeriv(state0);
  Field3D conduction1 = this->getDeriv(state1);
  Field3D conduction2 = this->getDeriv(state2);
  BOUT_FOR_SERIAL(i, conduction1.getRegion("RGN_NOBNDRY")) {
    // Conduction is inversely proportional to the collision frequency
    EXPECT_NE(conduction1[i], 0.);
    EXPECT_DOUBLE_EQ(conduction0[i], 2 * conduction1[i]);
    EXPECT_DOUBLE_EQ(2 * conduction2[i], conduction1[i]);
  }
}

TEST_F(BraginskiiConductionTest, ConductionCollisionsMode) {
  BraginskiiConduction component_braginskii = this->makeComponent(1., "braginskii");
  BraginskiiConduction component_multispecies = this->makeComponent(1., "multispecies");
  Options state_brag = this->state_base.copy();
  state_brag["species"]["test+"]["charge"] = 1;
  state_brag["species"]["test+"]["temperature"] = this->temp1;
  state_brag["species"]["test+"]["collision_frequencies"]["test+_test+_coll"] = 0.5;
  state_brag["species"]["test+"]["collision_frequencies"]["test+_e_coll"] = 0.5;
  state_brag["species"]["test+"]["collision_frequency"] = 1.0;
  Options state_multi = state_brag.copy();

  component_braginskii.declareAllSpecies({"test+"});
  component_braginskii.transform(state_brag);
  component_multispecies.declareAllSpecies({"test+", "e"});
  component_multispecies.transform(state_multi);

  Field3D conduction_brag = this->getDeriv(state_brag);
  Field3D conduction_multi = this->getDeriv(state_multi);
  BOUT_FOR_SERIAL(i, conduction_brag.getRegion("RGN_NOBNDRY")) {
    EXPECT_NE(conduction_brag[i], 0.);
    EXPECT_DOUBLE_EQ(conduction_brag[i], 2 * conduction_multi[i]);
  }
}

namespace {
// state for the energy channel tests, with an energy source already there or not
Options conductionChannelState(const Options& base, const Field3D& T,
                               bool with_source) {
  Options state = base.copy();
  state["species"]["test+"]["charge"] = 1;
  state["species"]["test+"]["collision_frequencies"]["test+_test+_coll"] = 1.;
  state["species"]["test+"]["collision_frequency"] = 1.;
  state["species"]["test+"]["temperature"] = T;
  if (with_source) {
    state["species"]["test+"]["energy_source"] = 0.3 + 0.1 * T;
  }
  return state;
}
} // namespace

TEST_F(BraginskiiConductionTest, ConductionChannelIsEnergySourceChange) {
  BraginskiiConduction component = this->makeComponent(1.);
  component.declareAllSpecies({"test+"});

  // no source before, so the channel is the whole energy_source
  Options state = conductionChannelState(this->state_base, this->temp1, false);
  component.transform(state);
  Options outputs = this->output_base.copy();
  component.outputVars(outputs);
  const Field3D E_cond = get<Field3D>(outputs["Etest+_conduction"]);
  EXPECT_TRUE(IsFieldEqual(E_cond, this->getDeriv(state), "RGN_NOBNDRY", 0.0));
  EXPECT_GT(max(abs(E_cond), false, "RGN_NOBNDRY"), 0.0);
  ASSERT_EQ("test+", outputs["Etest+_conduction"].attributes["species"].as<std::string>());

  // with a source already there it's the change (before from a fresh state, a copy
  // would share data)
  Options state2 = conductionChannelState(this->state_base, this->temp1, true);
  const Field3D before = conductionChannelState(this->state_base, this->temp1, true)
                             ["species"]["test+"]["energy_source"]
                                 .as<Field3D>();
  component.transform(state2);
  Options outputs2 = this->output_base.copy();
  component.outputVars(outputs2);
  const Field3D change = this->getDeriv(state2) - before;
  EXPECT_TRUE(IsFieldEqual(get<Field3D>(outputs2["Etest+_conduction"]), change,
                           "RGN_NOBNDRY", 1e-12));
}

TEST_F(BraginskiiConductionTest, ConductionChannelIsACopy) {
  BraginskiiConduction component = this->makeComponent(1.);
  component.declareAllSpecies({"test+"});

  Options state = conductionChannelState(this->state_base, this->temp1, false);
  component.transform(state);
  const Field3D expected = copy(this->getDeriv(state));

  // a later component adding in place (like the sheath) mustn't change the channel
  Field3D energy_source = state["species"]["test+"]["energy_source"].as<Field3D>();
  BOUT_FOR_SERIAL(i, energy_source.getRegion("RGN_ALL")) { energy_source[i] += 1.0; }

  Options outputs = this->output_base.copy();
  component.outputVars(outputs);
  EXPECT_TRUE(IsFieldEqual(get<Field3D>(outputs["Etest+_conduction"]), expected,
                           "RGN_NOBNDRY", 0.0));
}

// the flutter operator communicates fluxes, so it needs two x guard cells
using BraginskiiConductionFlutterTest = FakeMeshFixture_tmpl<6, 5, 7, false, 2>;

TEST_F(BraginskiiConductionFlutterTest, ConductionChannelWithFlutter) {
  Options opts{{"units",
                {{"eV", 1.0}, {"meters", 1.0}, {"seconds", 1.0}, {"inv_meters_cubed", 1e19}}},
               {"test+",
                {{"type", "evolve_energy"},
                 {"thermal_conduction", true},
                 {"diagnose", true},
                 {"conduction_collisions_mode", "braginskii"}}}};
  const Options base{
      {"fastest_wave", 0.},
      {"species",
       {{"test+", {{"velocity", 0.}, {"density", 1.}, {"pressure", 1.}, {"AA", 2.}}}}}};
  const Field3D T = BraginskiiConductionTest::linearGradient(1., 0.1, 0., 1., 1., 0.2, 0.);

  BraginskiiConduction component("test+", opts, nullptr);
  component.declareAllSpecies({"test+"});
  Options state = conductionChannelState(base, T, false);
  Options state_no_flutter = state.copy();
  state["fields"]["Apar_flutter"] =
      BraginskiiConductionTest::linearGradient(0., 0.1, 0., 0., 0.05, 0.2, 0.);
  component.transform(state);
  Options outputs{{"Nnorm", 1e19}, {"Tnorm", 1.}, {"Omega_ci", 1.}, {"rho_s0", 1.}};
  component.outputVars(outputs);

  // both adds end up in the channel
  const Field3D with_flutter = state["species"]["test+"]["energy_source"].as<Field3D>();
  EXPECT_TRUE(IsFieldEqual(get<Field3D>(outputs["Etest+_conduction"]), with_flutter,
                           "RGN_NOBNDRY", 0.0));

  // and the flutter part isn't zero, or this tests nothing
  BraginskiiConduction component2("test+", opts, nullptr);
  component2.declareAllSpecies({"test+"});
  component2.transform(state_no_flutter);
  const Field3D without =
      state_no_flutter["species"]["test+"]["energy_source"].as<Field3D>();
  EXPECT_GT(max(abs(with_flutter - without), false, "RGN_NOBNDRY"), 0.0);
}

TEST_F(BraginskiiConductionTest, ConductionDiagnoseDoesNotChangeState) {
  Options opts_off = options.copy();
  opts_off["test+"]["diagnose"].force(false);
  opts_off["test+"]["kappa_coefficient"] = 1.;
  opts_off["test+"]["conduction_collisions_mode"] = "braginskii";
  BraginskiiConduction component_off("test+", opts_off, nullptr);
  component_off.declareAllSpecies({"test+"});
  BraginskiiConduction component_on = this->makeComponent(1.);
  component_on.declareAllSpecies({"test+"});

  Options state_off = conductionChannelState(this->state_base, this->temp1, true);
  Options state_on = conductionChannelState(this->state_base, this->temp1, true);
  component_off.transform(state_off);
  component_on.transform(state_on);

  Options outputs_off = this->output_base.copy();
  component_off.outputVars(outputs_off);
  EXPECT_FALSE(outputs_off.isSet("Etest+_conduction"));

  for (const std::string var : {"energy_source", "energy_flow_ylow", "kappa_par"}) {
    const Field3D on = state_on["species"]["test+"][var].as<Field3D>();
    const Field3D off = state_off["species"]["test+"][var].as<Field3D>();
    EXPECT_TRUE(IsFieldEqual(on, off, "RGN_NOBNDRY", 0.0)) << var;
  }
}

// TODO: Add tests for more types of collisions (e.g., electron-ion, charge-exchante, ion
// neutral...)
