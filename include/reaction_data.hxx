#pragma once
#ifndef REACTION_DATA_H
#define REACTION_DATA_H

#include <bout/bout_types.hxx>
#include <bout/generic_factory.hxx>
#include <bout/options.hxx>
#include <filesystem>
#include <memory>
#include <vector>

#include "hermes_build_config.hxx"

BOUT_ENUM_CLASS(RateParamsTypes, T, ET, nT, undefined)
BOUT_ENUM_CLASS(ReactionCoeffTypes, sigma_v, sigma_v_E)
BOUT_ENUM_CLASS(ReactionDataTypes, ADAS, Amjuel, HydHel, undefined)

/**
 * @brief Get the json database location from options, or default to a standard location
 * in the repo, relative to this header
 *
 * @param options Options object
 * @return std::filesystem::path the path to the database directory
 */
std::filesystem::path get_json_db_dir(Options& options);

/**
 * @brief Abstract base class that defines the interface for all reaction data objects.
 *
 */
class ReactionData {
public:
  virtual ~ReactionData() = default;

  /**
   * @brief Construct a new Reaction Data object
   *
   * @param type The type of reaction data source (e.g., ADAS, Amjuel)
   * @param data_label Label/ID associated with the dataset (e.g. "H.2_3.1.8" for
   * type=Amjuel)
   * @param[optional] metadata_keys List of keys that a derived class should populate in
   * the metadata map.
   */
  explicit ReactionData(ReactionDataTypes type, std::string data_label,
                        std::vector<std::string> metadata_keys = {})
      : type(type), data_label(std::move(data_label)),
        metadata_keys(std::move(metadata_keys)){};

  /**
   * @brief Subclasses must implement retrieval/calculation of coefficients.
   *
   * @return const reference to the (1D or) 2D rate coefficients.
   */
  virtual const std::vector<std::vector<BoutReal>>& get_coeffs() const = 0;

  /**
   * @brief Get the fit type (e.g. T, nT, ET)
   *
   * @return The RateParamsTypes value
   */
  RateParamsTypes get_fit_type() const;

  /**
   * @brief Evaluate <sigma.v> at a particular energy and temperature.
   *
   * @details Subclasses override eval_sigma_v_ET_impl to provide the actual
   * implementation.
   *
   * @param E a energy
   * @param T a temperature
   * @return BoutReal <sigma.v>(E,T)
   */
  BoutReal eval_sigma_v_ET(BoutReal E, BoutReal T);

  /// @brief Subclasses provide an implementation of eval_sigma_v_ET
  virtual BoutReal eval_sigma_v_ET_impl(BoutReal E, BoutReal T) = 0;

  /**
   * @brief Evaluate <sigma.v> at a particular density and temperature.
   *
   * @details Subclasses override eval_sigma_v_T_impl to provide the actual
   * implementation.
   *
   * @param T a temperature
   * @return BoutReal <sigma.v>(T)
   */
  BoutReal eval_sigma_v_T(BoutReal T);

  /// @brief Subclasses provide an implementation of eval_sigma_v_T
  virtual BoutReal eval_sigma_v_T_impl(BoutReal T) = 0;

  /**
   * @brief Evaluate <sigma.v> at a particular density and temperature.
   *
   * @details Subclasses override eval_sigma_v_nT_impl to provide the actual
   * implementation.
   *
   * @param T a temperature
   * @param n a density
   * @return BoutReal <sigma.v>(n,T)
   */
  BoutReal eval_sigma_v_nT(BoutReal T, BoutReal n);

  /// @brief Subclasses provide an implementation of eval_sigma_v_nT
  virtual BoutReal eval_sigma_v_nT_impl(BoutReal T, BoutReal n) = 0;

  /**
   * @brief Evaluate <sigma.v.E> at a particular density and temperature
   *
   * @details Subclasses override eval_sigma_vE_nT_impl to provide the actual
   * implementation.
   *
   * @param T a temperature
   * @param n a density
   * @return BoutReal <sigma.v.E>(n,T)
   */
  BoutReal eval_sigma_vE_nT(BoutReal T, BoutReal n);

  /// @brief Subclasses provide an implementation of eval_sigma_vE_nT
  virtual BoutReal eval_sigma_vE_nT_impl(BoutReal T, BoutReal n) = 0;

  /// @brief Get a (real-valued) metadata value by key
  BoutReal get_metadata(const std::string& key) const;

  /// @brief Get the type of this reaction data (e.g. Amjuel, ADAS)
  ReactionDataTypes get_type() const { return this->type; }

  std::string src_str();

protected:
  /// The type of parameterization used for the reaction rate fit
  RateParamsTypes fit_type = RateParamsTypes::undefined;

  /// The type of reaction data source
  const ReactionDataTypes type;

  /// Label/ID associated with the dataset (e.g. "H.2_3.1.8" for type=Amjuel)
  const std::string data_label;

  const std::vector<std::string> metadata_keys;

  /**
   * @brief Allow subclasses to set metadata values.
   *
   * @param key
   * @param value
   */
  void set_metadata(const std::string& key, BoutReal value) {
    this->metadata[key] = value;
  }

private:
  /// Reaction metadata; only real values supported for now
  std::map<std::string, BoutReal> metadata;
};

/**
 * @brief Base for reaction data classes that store coefficients.
 *
 */
class ReactionDataWithCoeffs : public ReactionData {
public:
  /**
   * @brief Construct a ReactionDataWithCoeffs object
   *
   * @param type The type of reaction data source (e.g., ADAS, Amjuel)
   * @param data_label Label/ID associated with the dataset (e.g. "H.2_3.1.8" for
   * type=Amjuel)
   */
  explicit ReactionDataWithCoeffs(ReactionDataTypes type, std::string data_label,
                                  std::vector<std::string> metadata_keys = {});

  /**
   * @brief Get the coefficients
   *
   * @return const reference to the coefficient matrix
   */
  const std::vector<std::vector<BoutReal>>& get_coeffs() const final;

protected:
  /// The coefficients for the reaction rate fit
  std::vector<std::vector<BoutReal>> coeffs;
};

/// Factory for creating ReactionData objects
class ReactionDataFactory
    : public Factory<ReactionData, ReactionDataFactory, const std::string&, Options&,
                     const std::vector<std::string>> {
public:
  static constexpr auto type_name = "ReactionData";
  static constexpr auto section_name = "reactions";
  static constexpr auto option_name = "type";
  static constexpr auto default_type = "none";
};

/// More concise alias for registering a reaction data type
template <typename DerivedType>
using RegisterReactionData = ReactionDataFactory::RegisterInFactory<DerivedType>;

#endif
