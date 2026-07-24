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
#include "reaction_data.hxx"

namespace hermes {

/// Struct to hold values used in cell averaging
struct DataTriple {
  DataTriple() : DataTriple(0, 0, 0) {}
  DataTriple(BoutReal c, BoutReal l, BoutReal r) : centre(c), left(l), right(r) {}
  BoutReal centre, left, right;
};

// Temporary storage for central,left,right vals of various properties at each cell index
struct CellProps {
  CellProps(size_t num_reactants)
      : densities(num_reactants), collision_freqs(num_reactants) {}

  DataTriple mass_action = {1, 1, 1};
  DataTriple rate;

  std::vector<DataTriple> densities;
  std::vector<DataTriple> collision_freqs;
};

/// Signatures for different rate calculations.
/// N.B. one extra arg required for the mass action factor.
using OneDRateFunc = std::function<BoutReal(BoutReal, BoutReal)>;
using TwoDRateFunc = std::function<BoutReal(BoutReal, BoutReal, BoutReal)>;
using RateFuncVariant = std::variant<OneDRateFunc, TwoDRateFunc>;

/// Struct to hold final (possibly averaged) reaction rate and collision frequency fields that get returned by calc_rates()
struct RateData {
  /// The reaction rate Field3D
  Field3D rate;
  /// Collision frequencies keyed by reactant name
  std::map<std::string, Field3D> collision_frequencies;

  /**
   * @brief Extract collision frequency for a reactant.
   *
   * @param reactant_name Name of the reactant
   * @return const Field3D& The collision frequency field
   * @throws BoutException if reactant name not found
   */
  const Field3D& coll_freq(const std::string& reactant_name) const {
    auto it = this->collision_frequencies.find(reactant_name);
    if (it == this->collision_frequencies.end()) {
      throw BoutException(
          fmt::format("Collision frequency not found for reactant '{}'", reactant_name));
    }
    return it->second;
  }
};

// This is a workaround before CWG2518/P2593R1, taken from cppreference.com
template <RateParamsTypes>
constexpr bool dependent_false = false;

// Name used to store effective temperature (RateParamsTypes::T)
static const std::string Teff_name = "Teff";

/**
 * @brief Struct to encapsulate reaction rate and collision frequency calculations for a
 * number of different parameterisations.
 *
 * @tparam RateParamsType type identifying the reaction rate function parameters.
 */
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
  RateHelper(const GuardedOptions state, const Options& units,
             const std::vector<std::string>& reactant_names, const Region<Ind3D> region)
      : num_reactants(reactant_names.size()), reactant_names(reactant_names),
        region(region) {

    // Extract and store reactant densities
    for (const auto& reactant : reactant_names) {
      // Hold the wrapper in a local so the referenced field lifetime is obvious to GCC.
      const GuardedOptions reactant_density = state["species"][reactant]["density"];
      this->reactant_densities.push_back(&reactant_density.GetRef<Field3D>());
    }

    // Compute / extract fields that are required as parameters for the rate calculations
    if constexpr (RateParamsType == RateParamsTypes::ET) {
      static_assert(dependent_false<RateParamsType>,
                    "RateParamsTypes::ET not implemented");
    } else if constexpr (RateParamsType == RateParamsTypes::nT) {

      const GuardedOptions n_e = state["species"]["e:density"];
      const GuardedOptions T_e = state["species"]["e:temperature"];
      this->rate_params.push_back(&n_e.GetRef<Field3D>());
      this->rate_params.push_back(&T_e.GetRef<Field3D>());
    } else if constexpr (RateParamsType == RateParamsTypes::T) {
      calc_Teff(state, units, reactant_names, Teff_storage);
      this->rate_params.push_back(&Teff_storage);
    } else {
      // Compile-time error if any other RateParamsType enum exists
      static_assert(dependent_false<RateParamsType>, "Unhandled RateParamsType");
    }
  }

  /**
   * @brief Compute the cell-averaged reaction rate and collision frequencies, accounting
   * for reactant densities.
   *
   * @param rate_calc_func_variant a function that calculates the rate. Typed as
   * std::variant to easily switch between different rate parameterisations.
   * @param do_averaging whether to perform cell averaging
   * @return RateData containing the calculated rates and collision frequencies
   */
  RateData calc_rates(const RateFuncVariant& rate_calc_func_variant,
                      bool do_averaging = true) {

    // Initialize result data structure
    RateData result;
    result.rate = emptyFrom(*this->rate_params[0]);
    for (const std::string& reactant : this->reactant_names) {
      result.collision_frequencies[reactant] = emptyFrom(*this->rate_params[0]);
    }

    // Temporary storage for densities, rate and collision frequencies associated with a cell
    CellProps cell_props(this->num_reactants);

    // Populate cell_props differently according to the rate function type
    std::visit(
        [this, &cell_props, &do_averaging, &result](auto&& rate_calc_func) {
          auto J = result.rate.getCoordinates()->J;
          BOUT_FOR(i, region) {
            // Get densities and mass action-factor for this cell
            add_densities_and_mass_action(i, cell_props, do_averaging);

            using RateFuncType = std::decay_t<decltype(rate_calc_func)>;
            if constexpr (std::is_same_v<RateFuncType, OneDRateFunc>) {
              // 1D rate function
              if constexpr (RateParamsType == RateParamsTypes::T) {
                // Get rate params
                DataTriple T_data = compute_rate_param(0, i, do_averaging);

                // Compute rate(s)
                cell_props.rate.centre =
                    rate_calc_func(cell_props.mass_action.centre, T_data.centre);
                if (do_averaging) {
                  cell_props.rate.left =
                      rate_calc_func(cell_props.mass_action.left, T_data.left);
                  cell_props.rate.right =
                      rate_calc_func(cell_props.mass_action.right, T_data.right);
                }
              } else {
                throw BoutException(
                    "Unhandled RateParamsType (1D rate function being passed)");
              }
            } else if constexpr (std::is_same_v<RateFuncType, TwoDRateFunc>) {
              // 2D rate function
              if constexpr (RateParamsType == RateParamsTypes::nT) {
                // Get rate params
                DataTriple ne_data = compute_rate_param(0, i, do_averaging);
                DataTriple Te_data = compute_rate_param(1, i, do_averaging);

                // Compute rate(s)
                cell_props.rate.centre = rate_calc_func(cell_props.mass_action.centre,
                                                        ne_data.centre, Te_data.centre);
                if (do_averaging) {
                  cell_props.rate.left = rate_calc_func(cell_props.mass_action.left,
                                                        ne_data.left, Te_data.left);
                  cell_props.rate.right = rate_calc_func(cell_props.mass_action.right,
                                                         ne_data.right, Te_data.right);
                }

              } else if constexpr (RateParamsType == RateParamsTypes::ET) {
                static_assert(
                    dependent_false<RateParamsType>,
                    "RateHelper::calc_rate not set up for RateParamsTypes::ET yet");
              } else {
                throw BoutException(
                    "Unhandled RateParamsType (2D rate function being passed)");
              }
            }

            // Add the collision freqs (rate(s) divided by each reactant density)
            add_collision_freqs(cell_props, do_averaging);

            if (do_averaging) {
              // Compute averaged rate, collision frequencies and store in result
              auto Ji = J[i];
              auto Jm = J.yup()[i.ym()];
              auto Jp = J.ydown()[i.yp()];

              result.rate[i] = 4. / 6 * cell_props.rate.centre
                               + (Ji + Jm) / (12. * Ji) * cell_props.rate.left
                               + (Ji + Jp) / (12. * Ji) * cell_props.rate.right;
              for (std::size_t reactant_idx = 0; reactant_idx < this->num_reactants;
                   reactant_idx++) {

                std::string reactant_name = this->reactant_names[reactant_idx];
                result.collision_frequencies[reactant_name][i] =
                    4. / 6 * cell_props.collision_freqs[reactant_idx].centre
                    + (Ji + Jm) / (12. * Ji)
                          * cell_props.collision_freqs[reactant_idx].left
                    + (Ji + Jp) / (12. * Ji)
                          * cell_props.collision_freqs[reactant_idx].right;
              }
            } else {
              // Extract rate, collision frequencies at cell centre and store in result
              result.rate[i] = cell_props.rate.centre;
              for (std::size_t reactant_idx = 0; reactant_idx < this->num_reactants;
                   reactant_idx++) {
                std::string reactant_name = this->reactant_names[reactant_idx];
                result.collision_frequencies[reactant_name][i] =
                    cell_props.collision_freqs[reactant_idx].centre;
              }
            }
          }
        },
        rate_calc_func_variant);
    return result;
  }

private:
  /// Size of reactant_names, cached to avoid repeated .size() calls
  size_t num_reactants;

  /// Function to calculate reaction rate
  RateFuncVariant rate_calc_func;

  /// Reactant densities, keyed by species name (stored as pointers to avoid copying)
  std::vector<const Field3D*> reactant_densities;

  /// Reactant names in the order provided to the constructor
  std::vector<std::string> reactant_names;

  /// region in which to calculate the rate
  const Region<Ind3D> region;

  // Rate parameter fields (stored as pointers to avoid copying)
  std::vector<const Field3D*> rate_params;

  // Storage for Teff when RateParamsType == T (needed to keep the object alive)
  Field3D Teff_storage;

  /**
   * @brief Compute the effective temperature (in eV) of heavy reactants.
   *
   * @details Used to scale different isotope masses and finite neutral particle
   * temperatures by using the effective temperature (Amjuel p43)
   *   T_eff = (M/M_1)T_1 + (M/M_2)T_2
   *
   * @param[in] state
   * @param[in] reactant_names names of all reactant species
   * @param[inout] Teff Field3D object in which to store the result
   *
   * @todo read clamp values from json
   */
  void calc_Teff(const GuardedOptions state, const Options& units,
                 const std::vector<std::string>& reactant_names, Field3D& Teff) {

    std::vector<std::string> heavy_reactant_names;
    std::copy_if(reactant_names.begin(), reactant_names.end(),
                 std::back_inserter(heavy_reactant_names),
                 [](const std::string s) { return s.compare("e") != 0; });
    Teff = 0.0;
    BoutReal Tnorm = get<BoutReal>(units["eV"]);
    for (auto& sp : heavy_reactant_names) {
      // Keep the child option alive before taking a reference to its stored field.
      const GuardedOptions temperature_opt = state["species"][sp]["temperature"];
      const Field3D& temperature = temperature_opt.template GetRef<Field3D>();
      BoutReal AA = get<BoutReal>(state["species"][sp]["AA"]);
      Teff += (temperature / AA) * Tnorm;
    }

    // Clamp values
    constexpr BoutReal Teff_min = 0.01;
    constexpr BoutReal Teff_max = 10000;
    for (const auto& i : Teff.getRegion("RGN_NOBNDRY")) {
      Teff[i] = std::clamp(Teff[i], Teff_min, Teff_max);
    }
  }

  /**
   * @brief For each reactant, divide the reaction rate by reactant density to compute the associated collision frequency
   *
   * @param cell_props[inout] temporary struct already containing reaction rates and densities
   * @param do_averaging[in] whether to compute left and right values for averaging
   */
  inline void add_collision_freqs(CellProps& cell_props, bool do_averaging = true) {
    for (std::size_t reactant_idx = 0; reactant_idx < this->num_reactants;
         reactant_idx++) {
      cell_props.collision_freqs[reactant_idx].centre =
          cell_props.rate.centre / cell_props.densities[reactant_idx].centre;
      if (do_averaging) {
        cell_props.collision_freqs[reactant_idx].left =
            cell_props.rate.left / cell_props.densities[reactant_idx].left;
        cell_props.collision_freqs[reactant_idx].right =
            cell_props.rate.right / cell_props.densities[reactant_idx].right;
      }
    }
  }

  inline void add_densities_and_mass_action(Ind3D i, CellProps& cell_props,
                                            bool do_averaging) {
    // Reset mass action products
    cell_props.mass_action.centre = 1.0;
    cell_props.mass_action.left = 1.0;
    cell_props.mass_action.right = 1.0;

    // Store interpolated densities and mass action
    for (std::size_t reactant_idx = 0; reactant_idx < this->num_reactants;
         ++reactant_idx) {
      auto* dens = this->reactant_densities[reactant_idx];
      cell_props.densities[reactant_idx].centre = (*dens)[i];
      cell_props.mass_action.centre *= (*dens)[i];
      if (do_averaging) {
        auto ym = i.ym();
        auto yp = i.yp();
        cell_props.densities[reactant_idx].left =
            cellLeft<Limiter>((*dens)[i], (*dens)[ym], (*dens)[yp]);
        cell_props.densities[reactant_idx].right =
            cellRight<Limiter>((*dens)[i], (*dens)[ym], (*dens)[yp]);
        cell_props.mass_action.left *= cell_props.densities[reactant_idx].left;
        cell_props.mass_action.right *= cell_props.densities[reactant_idx].right;
      }
    }
  }

  /**
   * @brief Extract the value of a rate parameter at centre, left, and right positions.
   *
   * @param name name of the parameter (label in state["species"])
   * @param i central index
   * @param do_averaging whether to compute left and right values
   * @return CellData struct containing centre, left, and right values
   */
  inline DataTriple compute_rate_param(std::size_t idx, Ind3D i, bool do_averaging) {
    auto ym = i.ym();
    auto yp = i.yp();
    const auto& field = *this->rate_params[idx];
    DataTriple result;
    result.centre = field[i];
    if (do_averaging) {
      result.left = cellLeft<Limiter>(field[i], field[ym], field[yp]);
      result.right = cellRight<Limiter>(field[i], field[ym], field[yp]);
    }
    return result;
  }
};

} // namespace hermes

#endif
