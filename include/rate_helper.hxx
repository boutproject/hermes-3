#pragma once
#ifndef RATE_HELPER_H
#define RATE_HELPER_H

#include <functional>

#include <bout/bout_types.hxx>
#include <bout/region.hxx>

#include "component.hxx"
#include "hermes_build_config.hxx"
#include "integrate.hxx"

typedef std::function<BoutReal(BoutReal, BoutReal, BoutReal)> RateFunctionType;
template <typename LimiterType = hermes::Limiter, typename IdxType = Ind3D>
struct RateHelper {
  /**
   * @brief Construct a new RateHelper, extracting and storing some fields from the state
   * for use later in the rate calculation.
   *
   * @tparam LimiterType
   * @tparam RegionType
   * @param state
   * @param reactant_names vector of reactant names
   * @param rate_calc_func function with which to compute the rate from the mass action
   * factor, n_e and T_e
   * @param region the region in which to calculate the rate
   */
  RateHelper(const Options& state, const std::vector<std::string>& reactant_names,
             RateFunctionType rate_calc_func, const Region<IdxType> region)
      : rate_calc_func(rate_calc_func), region(region) {

    // Extract electron properties from state
    const Options& electron = state["species"]["e"];
    n_e = get<Field3D>(electron["density"]);
    T_e = get<Field3D>(electron["temperature"]);

    // Extract and store reactant densities
    std::transform(reactant_names.begin(), reactant_names.end(),
                   std::back_inserter(n_reactants),
                   [&](const std::string& reactant_name) {
                     return get<Field3D>(state["species"][reactant_name]["density"]);
                   });
  }

  /**
   * @brief Compute the cell-averaged reaction rate, accounting for the mass action
   * factor (product of reactant densities)
   *
   * @tparam LimiterType
   * @tparam IdxType
   * @return Field3D the cell-averaged reaction rate
   */
  Field3D calc_rate() {
    Field3D reaction_rate{emptyFrom(n_e)};
    auto J = reaction_rate.getCoordinates()->J;
    BOUT_FOR(i, region) {

      auto yp = i.yp();
      auto ym = i.ym();
      auto Ji = J[i];

      reaction_rate[i] =
          4. / 6 * rate_calc_func(mass_action(i), n_e[i], T_e[i])
          + (Ji + J[ym]) / (12. * Ji)
                * rate_calc_func(mass_action_left(i, ym, yp),
                                 cellLeft<LimiterType>(n_e[i], n_e[ym], n_e[yp]),
                                 cellLeft<LimiterType>(T_e[i], T_e[ym], T_e[yp]))
          + (Ji + J[yp]) / (12. * Ji)
                * rate_calc_func(mass_action_right(i, ym, yp),
                                 cellRight<LimiterType>(n_e[i], n_e[ym], n_e[yp]),
                                 cellRight<LimiterType>(T_e[i], T_e[ym], T_e[yp]));
    }
    return reaction_rate;
  }

private:
  /// region in which to calculate the rate
  const Region<IdxType> region;
  /// Function to calculate reaction rate as a function of n_e, T_e
  RateFunctionType rate_calc_func;
  /// Electron density and temperature
  Field3D n_e;
  Field3D T_e;
  /// Reactant densities
  std::vector<Field3D> n_reactants;

  BoutReal mass_action(IdxType i) {
    BoutReal result = 1;
    for (const auto& n : n_reactants) {
      result *= n[i];
    }
    return result;
  }

  BoutReal mass_action_left(IdxType i, IdxType ym, IdxType yp) {
    BoutReal result = 1;
    for (const auto& n : n_reactants) {
      result *= cellLeft<LimiterType>(n[i], n[ym], n[yp]);
    }
    return result;
  }

  BoutReal mass_action_right(IdxType i, IdxType ym, IdxType yp) {
    BoutReal result = 1;
    for (const auto& n : n_reactants) {
      result *= cellRight<LimiterType>(n[i], n[ym], n[yp]);
    }
    return result;
  }
};

#endif