#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/fixed_fraction_radiation.hxx"

#include <bout/output_bout_types.hxx>

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// Reuse the "standard" fixture for FakeMesh
template <typename CoolingCurve>
class FixedFractionRadiationTest : public FakeMeshFixture {
public:
  Options options = {
      {"units",
       {{"eV", 1.0}, {"meters", 1.0}, {"seconds", 1.0}, {"inv_meters_cubed", 1.0}}}};
  using FixedFracRadType = FixedFractionRadiation<CoolingCurve>;
};

using CoolingCurves =
    ::testing::Types<HutchinsonCarbon, Argon_adas, Neon_adas, Nitrogen_adas, Carbon_adas,
                     Argon_simplified1, Argon_simplified2, Argon_simplified3,
                     Krypton_adas, Xenon_adas, Tungsten_adas>;

TYPED_TEST_SUITE(FixedFractionRadiationTest, CoolingCurves);

TYPED_TEST(FixedFractionRadiationTest, DefaultParameters) {
  FixedFractionRadiation<TypeParam> component("test", this->options, nullptr);
  Options state{{"species", {{"e", {{"density", 1.}, {"temperature", 10.}}}}}};

  component.transform(state);
  auto energy_source = get<Field3D>(state["species"]["e"]["energy_source"]);
  BOUT_FOR_SERIAL(i, energy_source.getRegion("RGN_NOBNDRY")) {
    EXPECT_EQ(energy_source[i], 0.);
  }

  Options output;
  component.outputVars(output);
  EXPECT_FALSE(output.isSet("Rtest"));
}

TYPED_TEST(FixedFractionRadiationTest, ScaleFraction) {
  Options options1 = this->options.copy();
  options1["test"]["fraction"] = 0.05;
  options1["test"]["no_core_radiation"] = false;
  Options options2 = this->options.copy();
  options2["test"]["fraction"] = 0.1;
  options2["test"]["no_core_radiation"] = false;
  FixedFractionRadiation<TypeParam> component1("test", options1, nullptr);
  FixedFractionRadiation<TypeParam> component2("test", options2, nullptr);

  Options state1{{"species", {{"e", {{"density", 1.}, {"temperature", 10.}}}}}};
  Options state2 = state1.copy();

  component1.transform(state1);
  component2.transform(state2);
  auto e1 = get<Field3D>(state1["species"]["e"]["energy_source"]);
  auto e2 = get<Field3D>(state2["species"]["e"]["energy_source"]);
  // Radiation scales linearly with fraction
  BOUT_FOR_SERIAL(i, e1.getRegion("RGN_NOBNDRY")) {
    EXPECT_NE(e1[i], 0.);
    EXPECT_DOUBLE_EQ(2 * e1[i], e2[i]);
  }
}

TYPED_TEST(FixedFractionRadiationTest, ScaleDensity) {
  this->options["test"]["fraction"] = 0.05;
  this->options["test"]["no_core_radiation"] = false;
  FixedFractionRadiation<TypeParam> component("test", this->options, nullptr);

  Options state1{{"species", {{"e", {{"density", 1.}, {"temperature", 10.}}}}}};
  Options state2{{"species", {{"e", {{"density", 2.}, {"temperature", 10.}}}}}};

  component.transform(state1);
  component.transform(state2);
  auto e1 = get<Field3D>(state1["species"]["e"]["energy_source"]);
  auto e2 = get<Field3D>(state2["species"]["e"]["energy_source"]);
  // Radiation scales quadratically with electron density
  BOUT_FOR_SERIAL(i, e1.getRegion("RGN_NOBNDRY")) {
    EXPECT_NE(e1[i], 0.);
    EXPECT_DOUBLE_EQ(4 * e1[i], e2[i]);
  }
}

TYPED_TEST(FixedFractionRadiationTest, ScaleTemperature) {
  this->options["test"]["fraction"] = 0.05;
  this->options["test"]["no_core_radiation"] = false;
  FixedFractionRadiation<TypeParam> component("test", this->options, nullptr);

  Options state1{{"species", {{"e", {{"density", 1.}, {"temperature", 0.1}}}}}};
  Options state2{{"species", {{"e", {{"density", 1.}, {"temperature", 2000.}}}}}};

  component.transform(state1);
  component.transform(state2);
  auto e1 = get<Field3D>(state1["species"]["e"]["energy_source"]);
  auto e2 = get<Field3D>(state2["species"]["e"]["energy_source"]);
  BOUT_FOR_SERIAL(i, e1.getRegion("RGN_NOBNDRY")) {
    if constexpr (TypeParam::type == "argon_simplified2"
                  or TypeParam::type == "argon_simplified3") {
      // At these extreme temperatures, argon_simplified2 and
      // argon_simplified3 set radiation to zero.
      EXPECT_EQ(e1[i], 0);
      EXPECT_EQ(e2[i], 0);
    } else {
      // Radiation is nonlinear with temperature and varies depending on
      // species, but should grow. As it represents an energy loss,
      // energy_source < 0. Therfore, more radiation equates to a more
      // negative value.
      EXPECT_LT(e2[i], e1[i]);
    }
  }
}

TYPED_TEST(FixedFractionRadiationTest, ScaleMultiplier) {
  this->options["test"]["fraction"] = 0.05;
  this->options["test"]["no_core_radiation"] = false;
  Options options1 = this->options.copy();
  options1["test"]["R_multiplier"] = 0.5;
  Options options2 = this->options.copy();
  options2["test"]["R_multiplier"] = 2.;
  FixedFractionRadiation<TypeParam> component1("test", options1, nullptr);
  FixedFractionRadiation<TypeParam> component2("test", options2, nullptr);

  Options state1{{"species", {{"e", {{"density", 1.}, {"temperature", 10.}}}}}};
  Options state2 = state1.copy();

  component1.transform(state1);
  component2.transform(state2);
  auto e1 = get<Field3D>(state1["species"]["e"]["energy_source"]);
  auto e2 = get<Field3D>(state2["species"]["e"]["energy_source"]);
  // Radiation scales linearly with R_multiplier
  BOUT_FOR_SERIAL(i, e1.getRegion("RGN_NOBNDRY")) {
    EXPECT_NE(e1[i], 0.);
    EXPECT_DOUBLE_EQ(4 * e1[i], e2[i]);
  }
}

TYPED_TEST(FixedFractionRadiationTest, Output) {
  this->options["test"]["fraction"] = 0.05;
  this->options["test"]["diagnose"] = true;
  this->options["test"]["no_core_radiation"] = false;
  FixedFractionRadiation<TypeParam> component("test", this->options, nullptr);

  Options state{{"species", {{"e", {{"density", 1.}, {"temperature", 10.}}}}}};
  component.transform(state);
  auto e = get<Field3D>(state["species"]["e"]["energy_source"]);

  Options output;
  component.outputVars(output);
  auto rad = get<Field3D>(output["Rtest"]);
  BOUT_FOR_SERIAL(i, e.getRegion("RGN_NOBNDRY")) {
    EXPECT_NE(e[i], 0.);
    EXPECT_DOUBLE_EQ(e[i], rad[i]);
  }
}

TYPED_TEST(FixedFractionRadiationTest, NoCoreRad) {
  this->options["test"]["fraction"] = 0.05;
  this->options["test"]["no_core_radiation"] = true;
  FixedFractionRadiation<TypeParam> component("test", this->options, nullptr);

  Options state{{"species", {{"e", {{"density", 1.}, {"temperature", 10.}}}}}};
  component.transform(state);
  auto e = get<Field3D>(state["species"]["e"]["energy_source"]);
  BOUT_FOR_SERIAL(i, e.getRegion("RGN_NOBNDRY")) { EXPECT_EQ(e[i], 0.); }
}

// Things to test:
// radiation is negative energy source
// scales linearly with electron density, fraction, rad_multiplier, Nnorm, cooling curve
// scales inversely with Tnorm, FreqNorm
// If no_core_radiation then 0 in core (i.e., where mesh is periodic); don't think that's testable with FakeMesh
// Sets Rspecies if diagnose
// Check max/min of radiation curves?

//428
