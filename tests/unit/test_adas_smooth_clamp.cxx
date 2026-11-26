// Tests for the smooth clamp of the arguments of OpenADAS tables to the
// table range (smoothClamp in hermes_utils.hxx, option
// [reactions] adas_smooth_clamp_width), using the tables in json_database.

#include "../../include/adas_reaction.hxx"
#include "../../include/hermes_utils.hxx"
#include "../../include/reaction_data.hxx"

#include <bout/boutexception.hxx>
#include <bout/options.hxx>
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <string>
#include <vector>

namespace {

std::string table(const std::string& name) {
  Options options;
  return (get_json_db_dir(options) / name).string();
}

// Tables of ionisation, recombination, radiation and charge exchange
// coefficients used by the carbon, lithium and neon reactions
const std::vector<std::string> tables{"scd96_ne.json", "acd96_ne.json", "plt96_ne.json",
                                      "prb96_ne.json", "ccd89_ne.json", "scd96_li.json",
                                      "ccd89_li.json", "acd96_c.json",  "ccd96_c.json"};

constexpr BoutReal width = 0.1;     // Default width [decades]
constexpr BoutReal T_inside = 10.0; // Inside the T range of all tables [eV]
constexpr BoutReal n_inside = 1e19; // Inside the n range of all tables [m^-3]

/// log10 of the coefficient as a function of u = log10(T) or u = log10(n)
BoutReal log10k(OpenADASRateCoefficient& coef, BoutReal u, bool vary_T) {
  const BoutReal value = std::pow(10., u);
  return std::log10(vary_T ? coef.evaluate(value, n_inside)
                           : coef.evaluate(T_inside, value));
}

} // namespace

TEST(SmoothClamp, ZeroWidthIsClamp) {
  for (const BoutReal v : {-3.0, -0.5, -0.4999, 0.0, 0.3, 1.4999, 1.5, 5.0}) {
    EXPECT_EQ(smoothClamp(v, -0.5, 1.5, 0.0), std::clamp(v, -0.5, 1.5)) << v;
  }
}

TEST(SmoothClamp, ExactOutsideTransitions) {
  const BoutReal lo = -0.5;
  const BoutReal hi = 1.5;
  const BoutReal w = 0.2;
  for (const BoutReal v : {-10.0, -0.6}) {
    EXPECT_EQ(smoothClamp(v, lo, hi, w), lo) << v;
  }
  for (const BoutReal v : {1.6, 10.0}) {
    EXPECT_EQ(smoothClamp(v, lo, hi, w), hi) << v;
  }
  for (const BoutReal v : {-0.4, 0.0, 0.7, 1.4}) {
    EXPECT_EQ(smoothClamp(v, lo, hi, w), v) << v;
  }
  // At the edges the result is moved inwards by width / 8
  EXPECT_NEAR(smoothClamp(lo, lo, hi, w), lo + w / 8, 1e-15);
  EXPECT_NEAR(smoothClamp(hi, lo, hi, w), hi - w / 8, 1e-15);
}

TEST(SmoothClamp, ContinuouslyDifferentiable) {
  const BoutReal lo = -0.5;
  const BoutReal hi = 1.5;
  const BoutReal w = 0.2;
  const BoutReal delta = 1e-7;
  // At the edges and at both ends of each transition. The second derivative
  // is at most 1 / w, so one-sided differences differ by at most delta / w
  // (5e-7); a kink, as in the hard clamp, gives a difference of 1.
  for (const BoutReal x : {lo - w / 2, lo, lo + w / 2, hi - w / 2, hi, hi + w / 2}) {
    const BoutReal f0 = smoothClamp(x, lo, hi, w);
    const BoutReal right = (smoothClamp(x + delta, lo, hi, w) - f0) / delta;
    const BoutReal left = (f0 - smoothClamp(x - delta, lo, hi, w)) / delta;
    EXPECT_NEAR(right, left, 1e-5) << "x=" << x;
  }
  // The hard clamp has a kink at the edges
  const BoutReal right = (smoothClamp(lo + delta, lo, hi, 0.0) - lo) / delta;
  EXPECT_NEAR(right, 1.0, 1e-6);
  EXPECT_EQ(smoothClamp(lo - delta, lo, hi, 0.0), lo);
}

TEST(SmoothClamp, MonotonicAndInRange) {
  const BoutReal lo = -0.5;
  const BoutReal hi = 1.5;
  const BoutReal w = 0.2;
  BoutReal previous = smoothClamp(lo - w, lo, hi, w);
  for (int i = 1; i <= 2400; ++i) {
    const BoutReal x = lo - w + i * 1e-3;
    const BoutReal f = smoothClamp(x, lo, hi, w);
    EXPECT_GE(f, previous) << "x=" << x;
    EXPECT_GE(f, lo);
    EXPECT_LE(f, hi);
    previous = f;
  }
}

TEST(AdasSmoothClamp, OptionDefault) {
  Options options;
  EXPECT_DOUBLE_EQ(adas_smooth_clamp_width(options), width);

  Options options_off{{"reactions", {{"adas_smooth_clamp_width", 0.0}}}};
  EXPECT_EQ(adas_smooth_clamp_width(options_off), 0.0);
}

TEST(AdasSmoothClamp, ZeroWidthIsHardClamp) {
  OpenADASRateCoefficient coef(table("scd96_ne.json"), 0);
  EXPECT_EQ(coef.smooth_clamp_width, 0.0);
  const BoutReal at_tmin = coef.evaluate(coef.Tmin, n_inside);
  for (const BoutReal T : {-1.0, 0.0, 0.5 * coef.Tmin, 0.99 * coef.Tmin}) {
    EXPECT_EQ(coef.evaluate(T, n_inside), at_tmin) << "T=" << T;
  }
  const BoutReal at_nmax = coef.evaluate(T_inside, coef.nmax);
  EXPECT_EQ(coef.evaluate(T_inside, 1.01 * coef.nmax), at_nmax);
}

TEST(AdasSmoothClamp, UnchangedAwayFromEdges) {
  // More than width / 2 decades away from the table edges, inside or
  // outside the table, the smooth clamp gives the same values as the hard clamp
  const BoutReal f = std::pow(10., 0.51 * width);
  for (const auto& name : tables) {
    OpenADASRateCoefficient hard(table(name), 0, 0.0);
    OpenADASRateCoefficient smooth(table(name), 0, width);
    for (const BoutReal T : {-1.0, 0.0, 0.1 * hard.Tmin, hard.Tmin / f, hard.Tmin * f,
                             T_inside, hard.Tmax / f, hard.Tmax * f, 10 * hard.Tmax}) {
      for (const BoutReal n :
           {0.0, hard.nmin / f, hard.nmin * f, n_inside, hard.nmax / f, hard.nmax * f}) {
        EXPECT_EQ(smooth.evaluate(T, n), hard.evaluate(T, n))
            << name << " T=" << T << " n=" << n;
      }
    }
  }
}

TEST(AdasSmoothClamp, ValueAtEdges) {
  // At the table edges log10(T) and log10(n) are moved inwards by width / 8
  const BoutReal f = std::pow(10., width / 8);
  for (const auto& name : tables) {
    OpenADASRateCoefficient hard(table(name), 0, 0.0);
    OpenADASRateCoefficient smooth(table(name), 0, width);
    EXPECT_NEAR(smooth.evaluate(hard.Tmin, n_inside)
                    / hard.evaluate(hard.Tmin * f, n_inside),
                1.0, 1e-12)
        << name;
    EXPECT_NEAR(smooth.evaluate(hard.Tmax, n_inside)
                    / hard.evaluate(hard.Tmax / f, n_inside),
                1.0, 1e-12)
        << name;
    EXPECT_NEAR(smooth.evaluate(T_inside, hard.nmin)
                    / hard.evaluate(T_inside, hard.nmin * f),
                1.0, 1e-12)
        << name;
    EXPECT_NEAR(smooth.evaluate(T_inside, hard.nmax)
                    / hard.evaluate(T_inside, hard.nmax / f),
                1.0, 1e-12)
        << name;
  }
}

TEST(AdasSmoothClamp, ContinuouslyDifferentiableAtEdges) {
  // One-sided finite differences of log10(k) in u = log10(T) or log10(n), at
  // the four table edges and at both ends of each transition. Within a
  // transition the second derivative is s / width, where s is the slope of
  // the table in its edge cell, so for a continuously differentiable k the
  // one-sided differences differ by at most |s| * delta / width (1e-5 |s|).
  // The hard clamp has a jump |s| at the edges.
  const BoutReal delta = 1e-6;
  for (const auto& name : tables) {
    OpenADASRateCoefficient hard(table(name), 0, 0.0);
    OpenADASRateCoefficient smooth(table(name), 0, width);
    for (const bool vary_T : {true, false}) {
      const BoutReal lower = std::log10(vary_T ? hard.Tmin : hard.nmin);
      const BoutReal upper = std::log10(vary_T ? hard.Tmax : hard.nmax);
      for (const BoutReal edge : {lower, upper}) {
        // Slope of the table in the edge cell, from the hard clamp
        const BoutReal inwards = (edge == lower) ? 1.0 : -1.0;
        const BoutReal s = (log10k(hard, edge + inwards * width, vary_T)
                            - log10k(hard, edge + inwards * 0.75 * width, vary_T))
                           / (inwards * 0.25 * width);
        const BoutReal tolerance = 2 * std::abs(s) * delta / width + 1e-7;

        for (const BoutReal u : {edge - width / 2, edge, edge + width / 2}) {
          const BoutReal f0 = log10k(smooth, u, vary_T);
          const BoutReal right = (log10k(smooth, u + delta, vary_T) - f0) / delta;
          const BoutReal left = (f0 - log10k(smooth, u - delta, vary_T)) / delta;
          EXPECT_NEAR(right, left, tolerance)
              << name << (vary_T ? " T" : " n") << " u=" << u << " s=" << s;
          EXPECT_LE(std::abs(right), std::abs(s) * (1 + 1e-6) + 1e-7);
        }

        // The test detects the kink of the hard clamp, where the table has a slope
        if (std::abs(s) > 1e-3) {
          const BoutReal f0 = log10k(hard, edge, vary_T);
          const BoutReal right = (log10k(hard, edge + delta, vary_T) - f0) / delta;
          const BoutReal left = (f0 - log10k(hard, edge - delta, vary_T)) / delta;
          EXPECT_GT(std::abs(right - left), 0.5 * std::abs(s))
              << name << (vary_T ? " T" : " n") << " edge=" << edge;
        }
      }
    }
  }
}

TEST(AdasSmoothClamp, MonotonicAcrossEdges) {
  // In each transition the coefficient changes monotonically between its
  // values at the ends of the transition, in the direction of the table's
  // slope in the edge cell (no overshoot)
  for (const auto& name : tables) {
    OpenADASRateCoefficient hard(table(name), 0, 0.0);
    OpenADASRateCoefficient smooth(table(name), 0, width);
    for (const bool vary_T : {true, false}) {
      for (const BoutReal edge : {std::log10(vary_T ? hard.Tmin : hard.nmin),
                                  std::log10(vary_T ? hard.Tmax : hard.nmax)}) {
        const BoutReal first = log10k(smooth, edge - width / 2, vary_T);
        const BoutReal last = log10k(smooth, edge + width / 2, vary_T);
        const BoutReal sign = (last >= first) ? 1.0 : -1.0;
        BoutReal previous = first;
        for (int i = 1; i <= 100; ++i) {
          const BoutReal u = edge - width / 2 + i * width / 100;
          const BoutReal value = log10k(smooth, u, vary_T);
          EXPECT_GE(sign * (value - previous), -1e-12)
              << name << (vary_T ? " T" : " n") << " u=" << u;
          previous = value;
        }
      }
    }
  }
}

TEST(AdasSmoothClamp, InvalidWidthThrows) {
  EXPECT_THROW(OpenADASRateCoefficient(table("scd96_ne.json"), 0, -0.1), BoutException);
  EXPECT_THROW(OpenADASRateCoefficient(table("scd96_ne.json"), 0, 6.0), BoutException);
  EXPECT_THROW(OpenADASRateCoefficient(table("scd96_ne.json"), 0,
                                       std::numeric_limits<BoutReal>::quiet_NaN()),
               BoutException);
}
