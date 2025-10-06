#pragma once
#ifndef RATE_HELPER_H
#define RATE_HELPER_H

#include <functional>
#include <variant>

#include <bout/bout_types.hxx>
#include <bout/region.hxx>

#include "component.hxx"
#include "hermes_build_config.hxx"
#include "hermes_utils.hxx"
#include "integrate.hxx"

BOUT_ENUM_CLASS(RateParamsTypes, T, ET, nT)

/// Signatures for different rate calculations.
/// N.B. one extra arg required for the mass action factor.
using OneDRateFunc = std::function<BoutReal(BoutReal, BoutReal)>;
using TwoDRateFunc = std::function<BoutReal(BoutReal, BoutReal, BoutReal)>;
using RateFuncVariant = std::variant<OneDRateFunc, TwoDRateFunc>;

// This is a workaround before CWG2518/P2593R1, taken from cppreference.com
template <RateParamsTypes>
constexpr bool dependent_false = false;

// Name used to store effective temperature (RateParamsTypes::T)
static const std::string Teff_name = "Teff";

template <RateParamsTypes RateParamsType>
struct RateHelper {
  /**
   * @brief Construct a new RateHelper, extracting and storing some fields from the state
   * for use later in the rate calculation.
   *
   * @param state
   * @param reactant_names vector of reactant names
   * @param rate_calc_func function with which to compute the rate from the mass action
   * factor, n_e and T_e
   * @param region the region in which to calculate the rate
   */
  RateHelper(const Options& state, const std::vector<std::string>& reactant_names,
             const Region<Ind3D> region)
      : region(region) {

    // Compute / extract fields that are required as parameters for the rate calculations
    if constexpr (RateParamsType == RateParamsTypes::ET) {
      static_assert(dependent_false<RateParamsType>,
                    "RateParamsTypes::ET not implemented");
      // add_rate_param("e:density", get<Field3D>(state["species"]["e:density"]));
      // Field3D energy = ;
      // add_rate_param("e:energy", energy);
    } else if constexpr (RateParamsType == RateParamsTypes::nT) {
      for (auto field_id : {"e:density", "e:temperature"}) {
        add_rate_param(field_id, get<Field3D>(state["species"][field_id]));
      }
    } else if constexpr (RateParamsType == RateParamsTypes::T) {
      Field3D Teff;
      calc_Teff(state, reactant_names, Teff);
      add_rate_param("Teff", Teff);
    } else {
      // Compile-time error if any other RateParamsType enum exists
      static_assert(dependent_false<RateParamsType>, "Unhandled RateParamsType");
    }

    // Extract and store reactant densities
    std::transform(reactant_names.begin(), reactant_names.end(),
                   std::back_inserter(this->n_reactants),
                   [&](const std::string& reactant_name) {
                     return get<Field3D>(state["species"][reactant_name]["density"]);
                   });
  }

  void add_rate_param(const std::string& field_id, const Field3D& fld) {
    this->rate_params.insert(std::make_pair(field_id, fld));
  }

  /**
   * @brief Compute the cell-averaged reaction rate, accounting for the mass action
   * factor (product of reactant densities)
   *
   * @return Field3D the cell-averaged reaction rate
   */
  Field3D calc_rate(const RateFuncVariant& rate_calc_func_variant) {
    std::string first_key = str_keys(this->rate_params)[0];
    Field3D reaction_rate{emptyFrom(this->rate_params[first_key])};
    std::visit(
        [this, &reaction_rate](auto&& rate_calc_func) {
          auto J = reaction_rate.getCoordinates()->J;
          BOUT_FOR(i, region) {

            auto yp = i.yp();
            auto ym = i.ym();
            auto Ji = J[i];

            // Calc rate_central, rate_left, rate_right at this index
            BoutReal rate_central, rate_left, rate_right;

            using RateFuncType = std::decay_t<decltype(rate_calc_func)>;
            if constexpr (std::is_same_v<RateFuncType, OneDRateFunc>) {
              if constexpr (RateParamsType == RateParamsTypes::T) {
                rate_central =
                    rate_calc_func(mass_action(i), get_rate_param(Teff_name, i));
                rate_left = rate_calc_func(mass_action_left(i, ym, yp),
                                           get_rate_param_left(Teff_name, i, ym, yp));
                rate_right = rate_calc_func(mass_action_right(i, ym, yp),
                                            get_rate_param_right(Teff_name, i, ym, yp));
              } else {
                throw BoutException(
                    "Unhandled RateParamsType (1D rate function being passed)");
              }
            } else if constexpr (std::is_same_v<RateFuncType, TwoDRateFunc>) {
              if constexpr (RateParamsType == RateParamsTypes::nT) {
                rate_central =
                    rate_calc_func(mass_action(i), get_rate_param("e:density", i),
                                   get_rate_param("e:temperature", i));
                rate_left =
                    rate_calc_func(mass_action_left(i, ym, yp),
                                   get_rate_param_left("e:density", i, ym, yp),
                                   get_rate_param_left("e:temperature", i, ym, yp));
                rate_right =
                    rate_calc_func(mass_action_right(i, ym, yp),
                                   get_rate_param_right("e:density", i, ym, yp),
                                   get_rate_param_right("e:temperature", i, ym, yp));
              } else if constexpr (RateParamsType == RateParamsTypes::ET) {
                static_assert(
                    dependent_false<RateParamsType>,
                    "RateHelper::calc_rate not set up for RateParamsTypes::ET yet");
              } else {
                throw BoutException(
                    "Unhandled RateParamsType (2D rate function being passed)");
              }
            }

            // Overall rate at this index
            reaction_rate[i] = 4. / 6 * rate_central
                               + (Ji + J[ym]) / (12. * Ji) * rate_left
                               + (Ji + J[yp]) / (12. * Ji) * rate_right;
          }
        },
        rate_calc_func_variant);
    return reaction_rate;
  }

  /**
   * @brief Compute the effective temperature (in eV) of heavy reactants.
   *
   * @details Used to scale different isotope masses and finite neutral particle
   temperatures by using the effective temperature (Amjuel p43) T_eff = (M/M_1)T_1 +
   (M/M_2)T_2
   *
   * @param[in] state
   * @param[in] reactant_names names of all reactant species
   * @param[inout] Teff Field3D object in which to store the result
   *
   * @todo read clamp values from json?
   */
  void calc_Teff(const Options& state, const std::vector<std::string>& reactant_names,
                 Field3D& Teff) {

    std::vector<std::string> heavy_reactant_names;
    std::copy_if(reactant_names.begin(), reactant_names.end(),
                 std::back_inserter(heavy_reactant_names),
                 [](const std::string s) { return s.compare("e") != 0; });
    Teff = 0.0;
    BoutReal Tnorm = get<BoutReal>(state["units"]["eV"]);
    for (auto& sp : heavy_reactant_names) {
      Teff += (get<Field3D>(state["species"][sp]["temperature"])
               / get<Field3D>(state["species"][sp]["AA"]))
              * Tnorm;
    }

    // Clamp values
    constexpr BoutReal Teff_min = 0.01;
    constexpr BoutReal Teff_max = 10000;
    for (const auto& i : Teff.getRegion("RGN_NOBNDRY")) {
      Teff[i] = std::clamp(Teff[i], Teff_min, Teff_max);
    }
  }

  /**
   * @brief Extract the (cell-centre) value of a rate parameter.
   *
   * @param name name of the parameter (label in state["species"])
   * @param i central index
   * @return BoutReal
   */
  BoutReal get_rate_param(const std::string& name, Ind3D i) {
    return this->rate_params[name][i];
  }

  /**
   * @brief Extract the value of a rate parameter at the left of a cell.
   *
   * @param name name of the parameter (label in state["species"])
   * @param i central index
   * @param ym neighbour 1 index
   * @param yp neighbour 2 index
   * @return BoutReal
   */
  BoutReal get_rate_param_left(const std::string& name, Ind3D i, Ind3D ym, Ind3D yp) {
    return cellLeft<hermes::Limiter>(this->rate_params[name][i],
                                     this->rate_params[name][ym],
                                     this->rate_params[name][yp]);
  }

  /**
   * @brief Extract the value of a rate parameter at the right of a cell.
   * @param name name of the parameter (label in state["species"])
   * @param i central index
   * @param ym neighbour 1 index
   * @param yp neighbour 2 index
   * @return BoutReal
   */
  BoutReal get_rate_param_right(const std::string& name, Ind3D i, Ind3D ym, Ind3D yp) {
    return cellRight<hermes::Limiter>(this->rate_params[name][i],
                                      this->rate_params[name][ym],
                                      this->rate_params[name][yp]);
  }

private:
  /// region in which to calculate the rate
  const Region<Ind3D> region;
  /// Function to calculate reaction rate as a function of n_e, T_e
  RateFuncVariant rate_calc_func;

  /// Reactant densities
  std::vector<Field3D> n_reactants;

  // Rate parameter fields
  std::map<std::string, Field3D> rate_params;

  BoutReal mass_action(Ind3D i) {
    BoutReal result = 1;
    for (const auto& n : n_reactants) {
      result *= n[i];
    }
    return result;
  }

  BoutReal mass_action_left(Ind3D i, Ind3D ym, Ind3D yp) {
    BoutReal result = 1;
    for (const auto& n : n_reactants) {
      result *= cellLeft<hermes::Limiter>(n[i], n[ym], n[yp]);
    }
    return result;
  }

  BoutReal mass_action_right(Ind3D i, Ind3D ym, Ind3D yp) {
    BoutReal result = 1;
    for (const auto& n : n_reactants) {
      result *= cellRight<hermes::Limiter>(n[i], n[ym], n[yp]);
    }
    return result;
  }
};

#endif