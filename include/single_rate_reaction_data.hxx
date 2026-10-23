#pragma once
#ifndef SINGLE_RATE_REACTION_DATA_H
#define SINGLE_RATE_REACTION_DATA_H

#include <bout/bout_types.hxx>
#include <bout/options.hxx>
#include <string>
#include <vector>

#include "reaction_data.hxx"

namespace hermes {

const std::string SINGLE_RATE_REACTION_DATA_LBL = "single_rate";
/**
 * @brief Class to handle simple Reactions with a single fixed rate.
 *
 */
class SingleRateReactionData : public ReactionData {
public:
  /**
   * @brief Construct a new Single Rate Reaction Data object
   *
   * @param data_label ID/label for the specific data set
   * @param options Options object (Required for json database location)
   */
  SingleRateReactionData(
      const std::string& data_label, [[maybe_unused]] Options& options,
      [[maybe_unused]] const std::vector<std::string>& metadata_keys = {})
      : ReactionData(ReactionDataTypes::single_rate, data_label, {}), rate(1.0) {
    this->fit_type = RateParamsTypes::T;
  };

  /// Dummy implementation of get_coeffs
  const std::vector<std::vector<BoutReal>> dummy_coeffs = {{}};
  const std::vector<std::vector<BoutReal>>& get_coeffs() const final {
    return this->dummy_coeffs;
  };

protected:
  /**
   * @brief Evaluate <sigma . v . E> at a particular density and temperature
   *
   * @param T a temperature
   * @param n a density
   * @return BoutReal <sigma.v.E>(n,T)
   */
  BoutReal eval_sigma_vE_nT_impl([[maybe_unused]] BoutReal T,
                                 [[maybe_unused]] BoutReal n) final {
    return this->rate;
  };

  /**
   * @brief Evaluate <sigma.v> at a particular energy and temperature
   *
   * @param E a energy
   * @param T a temperature
   * @return BoutReal <sigma.v>(E,T)
   */
  BoutReal eval_sigma_v_ET_impl([[maybe_unused]] BoutReal E,
                                [[maybe_unused]] BoutReal T) final {
    return this->rate;
  };

  /**
   * @brief Evaluate <sigma.v> at a particular density and temperature
   *
   * @param T a temperature
   * @param n a density
   * @return BoutReal <sigma.v>(n,T)
   */
  BoutReal eval_sigma_v_nT_impl([[maybe_unused]] BoutReal T,
                                [[maybe_unused]] BoutReal n) final {
    return this->rate;
  };

  /**
   * @brief Evaluate <sigma.v> at a particular temperature
   *
   * @param T a temperature
   * @return BoutReal <sigma.v>(T)
   */
  BoutReal eval_sigma_v_T_impl([[maybe_unused]] BoutReal T) final { return this->rate; };

private:
  /// The single rate associated with this data object
  BoutReal rate;
};

/// Register with factory class
namespace {
RegisterReactionData<SingleRateReactionData>
    register_SingleRateReactionData(SINGLE_RATE_REACTION_DATA_LBL);
;
} // namespace

} // namespace hermes

#endif
