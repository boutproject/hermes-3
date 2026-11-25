// Tests for the optional Arrhenius continuation of OpenADAS coefficients
// below the table's minimum temperature (option
// [reactions] extrapolate_below_tmin), using the neon tables shipped in
// json_database.

#include "../../include/adas_reaction.hxx"
#include "../../include/reaction_data.hxx"

#include <bout/options.hxx>
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <string>

namespace {

std::string table(const std::string& name) {
  Options options;
  return (get_json_db_dir(options) / name).string();
}

// Density inside the table range [m^-3]
constexpr BoutReal n_eval = 1e19;

} // namespace

TEST(AdasExtrapolateBelowTmin, OptionDefaultIsFalse) {
  Options options;
  EXPECT_FALSE(adas_extrapolate_below_tmin(options));

  Options options_on{{"reactions", {{"extrapolate_below_tmin", true}}}};
  EXPECT_TRUE(adas_extrapolate_below_tmin(options_on));
}

TEST(AdasExtrapolateBelowTmin, DefaultClampsAtTmin) {
  OpenADASRateCoefficient scd(table("scd96_ne.json"), 0);
  EXPECT_FALSE(scd.extrapolate_below_tmin);
  EXPECT_DOUBLE_EQ(scd.evaluate(0.5 * scd.Tmin, n_eval), scd.evaluate(scd.Tmin, n_eval));
  EXPECT_DOUBLE_EQ(scd.evaluate(0.0, n_eval), scd.evaluate(scd.Tmin, n_eval));
}

TEST(AdasExtrapolateBelowTmin, UnchangedAtAndAboveTmin) {
  for (const auto* name : {"scd96_ne.json", "plt96_ne.json", "acd96_ne.json"}) {
    OpenADASRateCoefficient off(table(name), 0, false);
    OpenADASRateCoefficient on(table(name), 0, true);
    for (const BoutReal T : {off.Tmin, 1.0, 10.0, 100.0, off.Tmax, 2.0 * off.Tmax}) {
      for (const BoutReal n : {0.1 * off.nmin, 1e17, 1e19, 1e21, 10.0 * off.nmax}) {
        EXPECT_EQ(on.evaluate(T, n), off.evaluate(T, n))
            << name << " T=" << T << " n=" << n;
      }
    }
  }
}

TEST(AdasExtrapolateBelowTmin, EffFromLowestTwoTablePoints) {
  OpenADASRateCoefficient scd(table("scd96_ne.json"), 0, true);
  OpenADASRateCoefficient plt(table("plt96_ne.json"), 0, true);
  // ln(k1 / k0) / (1 / T0 - 1 / T1) at the lowest density, from the JSON
  // tables (T0 = 10^-0.699 = 0.2 eV, T1 = 10^-0.523 = 0.3 eV)
  EXPECT_NEAR(scd.E_eff, 22.1829, 1e-3);
  EXPECT_NEAR(plt.E_eff, 19.1423, 1e-3);
}

TEST(AdasExtrapolateBelowTmin, ArrheniusFactorBelowTmin) {
  OpenADASRateCoefficient scd(table("scd96_ne.json"), 0, true);
  ASSERT_GT(scd.E_eff, 0.0);
  for (const BoutReal frac : {0.9, 0.5, 0.25}) {
    const BoutReal T = frac * scd.Tmin;
    const BoutReal ratio = scd.evaluate(T, n_eval) / scd.evaluate(scd.Tmin, n_eval);
    const BoutReal expected = std::exp(-scd.E_eff * (1. / T - 1. / scd.Tmin));
    EXPECT_NEAR(ratio / expected, 1.0, 1e-12) << "T=" << T;
    EXPECT_LT(ratio, 1.0);
  }
}

TEST(AdasExtrapolateBelowTmin, ContinuousAtTmin) {
  OpenADASRateCoefficient scd(table("scd96_ne.json"), 0, true);
  const BoutReal at_tmin = scd.evaluate(scd.Tmin, n_eval);
  const BoutReal just_below = scd.evaluate(scd.Tmin * (1.0 - 1e-9), n_eval);
  EXPECT_NEAR(just_below / at_tmin, 1.0, 1e-6);
}

TEST(AdasExtrapolateBelowTmin, VanishesForNonPositiveT) {
  OpenADASRateCoefficient scd(table("scd96_ne.json"), 0, true);
  EXPECT_EQ(scd.evaluate(0.0, n_eval), 0.0);
  EXPECT_EQ(scd.evaluate(-1.0, n_eval), 0.0);
  // Far below Tmin the coefficient underflows to zero
  EXPECT_LT(scd.evaluate(1e-3, n_eval), 1e-300);
}

TEST(AdasExtrapolateBelowTmin, DensityStaysClamped) {
  OpenADASRateCoefficient scd(table("scd96_ne.json"), 0, true);
  const BoutReal T = 0.5 * scd.Tmin;
  EXPECT_DOUBLE_EQ(scd.evaluate(T, 1e-3 * scd.nmin), scd.evaluate(T, scd.nmin));
  EXPECT_DOUBLE_EQ(scd.evaluate(T, 1e3 * scd.nmax), scd.evaluate(T, scd.nmax));
}

TEST(AdasExtrapolateBelowTmin, NoContinuationIfEffNotPositive) {
  // Recombination coefficients increase towards low T, so the fitted
  // activation energy is negative; these stay clamped at Tmin.
  OpenADASRateCoefficient acd(table("acd96_ne.json"), 0, true);
  ASSERT_LT(acd.E_eff, 0.0);
  EXPECT_DOUBLE_EQ(acd.evaluate(0.5 * acd.Tmin, n_eval), acd.evaluate(acd.Tmin, n_eval));
  EXPECT_DOUBLE_EQ(acd.evaluate(0.0, n_eval), acd.evaluate(acd.Tmin, n_eval));
}

TEST(AdasExtrapolateBelowTmin, ChargeExchangeStaysClamped) {
  // OpenADASChargeExchange constructs its coefficient without the option, so
  // the charge exchange coefficient stays at its Tmin value although its
  // fitted activation energy is positive.
  OpenADASRateCoefficient ccd(table("ccd89_ne.json"), 0);
  EXPECT_FALSE(ccd.extrapolate_below_tmin);
  ASSERT_GT(ccd.E_eff, 0.0);
  EXPECT_DOUBLE_EQ(ccd.evaluate(0.5 * ccd.Tmin, n_eval), ccd.evaluate(ccd.Tmin, n_eval));
  EXPECT_DOUBLE_EQ(ccd.evaluate(0.0, n_eval), ccd.evaluate(ccd.Tmin, n_eval));
}
