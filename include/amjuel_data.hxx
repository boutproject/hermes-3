#pragma once
#ifndef AMJUEL_DATA_H
#define AMJUEL_DATA_H

#include <bout/bout_types.hxx>
#include <bout/options.hxx>
#include <string>
#include <vector>

#include "reaction_data.hxx"

namespace hermes {

/**
 * @brief Class to handle reading from json files containing Amjuel data, storing the
 * coefficients and evaluating associated cross sections, etc.
 *
 */
class AmjuelData : public ReactionDataWithCoeffs {
public:
  /**
   * @brief Read an Amjuel data file; store the coefficients and metadata.
   *
   * @param data_label ID/label for the specific data set
   * @param options Options object (Required for json database location)
   */
  AmjuelData(const std::string& data_label, Options& options,
             const std::vector<std::string> metadata_keys = {});

protected:
  /**
   * @brief Evaluate <sigma . v . E> at a particular density and temperature
   *
   * @param T a temperature
   * @param n a density
   * @return BoutReal <sigma.v.E>(n,T)
   */
  BoutReal eval_sigma_vE_nT_impl(BoutReal T, BoutReal n) final;

  /**
   * @brief Evaluate <sigma.v> at a particular energy and temperature
   *
   * @param E a energy
   * @param T a temperature
   * @return BoutReal <sigma.v>(E,T)
   */
  BoutReal eval_sigma_v_ET_impl(BoutReal E, BoutReal T) final;

  /**
   * @brief Evaluate <sigma.v> at a particular density and temperature
   *
   * @param T a temperature
   * @param n a density
   * @return BoutReal <sigma.v>(n,T)
   */
  BoutReal eval_sigma_v_nT_impl(BoutReal T, BoutReal n) final;

  /**
   * @brief Evaluate <sigma.v> at a particular temperature
   *
   * @param T a temperature
   * @return BoutReal <sigma.v>(T)
   */
  BoutReal eval_sigma_v_T_impl(BoutReal T) final;
};

/// Register with factory class
namespace {
RegisterReactionData<AmjuelData> register_amjueldata("amjuel");
}

} // namespace hermes

#endif
