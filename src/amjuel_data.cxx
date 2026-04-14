#include "amjuel_data.hxx"

#include <algorithm>
#include <cmath>
#include <fstream>

#include "../external/json.hxx"
#include <bout/boutexception.hxx>

namespace hermes {

/// @brief Helper functions for evaluating Amjuel fits.

/**
 * @brief Evaluate an Amjuel double polynomial fit in n and T, given a table of
 * coefficients (see page 20 of amjuel.pdf).
 *
 * @param T temperature in eV
 * @param n number density in m^-3
 * @param coeff_table polynomial fit coefficients in a vector-of-vectors (outer index T,
 * inner index n)
 * @return BoutReal the fit in SI, units m^3/s, or eV m^3/s for energy loss
 */
static BoutReal
eval_amjuel_nT_fit(BoutReal T, BoutReal n,
                   const std::vector<std::vector<BoutReal>>& coeff_table) {
  // Enforce range of validity
  n = std::clamp(n, 1e14, 1e22); // 1e8 - 1e16 cm^-3
  T = std::clamp(T, 0.1, 1e4);

  BoutReal logntilde = log(n / 1e14); // Note: 1e8 cm^-3
  BoutReal logT = log(T);
  BoutReal result = 0.0;

  BoutReal logT_n = 1.0; // log(T) ** n
  std::size_t num_n = coeff_table.size();
  std::size_t num_T = coeff_table[0].size();
  for (size_t n_idx = 0; n_idx < num_n; ++n_idx) {
    BoutReal logn_m = 1.0; // log(ntilde) ** m
    const auto& nrow = coeff_table[n_idx];
    for (size_t T_idx = 0; T_idx < num_T; ++T_idx) {
      result += nrow[T_idx] * logn_m * logT_n;
      logn_m *= logntilde;
    }
    logT_n *= logT;
  }
  return exp(result) * 1e-6; // Note: convert cm^3 to m^3
}

/**
 * @brief Evaluate an Amjuel single polynomial fit in T, given a table of
 * coefficients (see page 20 of amjuel.pdf).
 *
 * @param T temperature in eV
 * @param coeff_table a table of polynomial fit coefficients (index T)
 * @return BoutReal the fit in SI, units m^3/s, or eV m^3/s for energy loss
 */
static BoutReal eval_amjuel_T_fit(BoutReal T, const std::vector<BoutReal>& coeff_table) {
  const BoutReal lnT = log(T);
  BoutReal ln_sigmav = coeff_table[0];
  BoutReal lnT_n = lnT; // (lnT)^n
  const std::size_t num_T = coeff_table.size();
  for (std::size_t T_idx = 1; T_idx < num_T; T_idx++) {
    ln_sigmav += coeff_table[T_idx] * lnT_n;
    lnT_n *= lnT;
  }

  return exp(ln_sigmav);
}

///
AmjuelData::AmjuelData(const std::string& data_label, Options& options,
                       const std::vector<std::string> metadata)
    : ReactionDataWithCoeffs(ReactionDataTypes::Amjuel, data_label, metadata) {

  // Construct file path, checking that the directory exists
  auto data_dir = get_json_db_dir(options);
  if (!std::filesystem::is_directory(data_dir)) {
    throw BoutException("No json database found at '{:s}'", data_dir.string());
  }
  std::filesystem::path file_path = data_dir / ("AMJUEL_" + data_label + ".json");

  // Read the data file
  std::ifstream json_file(file_path);
  if (!json_file.good()) {
    throw BoutException("Failed to read Amjuel data file at '{:s}'",
                        std::string(file_path));
  }

  // Parse the data
  nlohmann::json data;
  json_file >> data;

  try {
    // Extract fit type (rate params)
    std::string str_fit_type = data["info"]["fit_type"];
    std::transform(str_fit_type.begin(), str_fit_type.end(), str_fit_type.begin(),
                   ::tolower);
    this->fit_type = RateParamsTypesFromString(str_fit_type);

    // Extract other metadata
    for (const auto& key : metadata) {
      if (data["info"].contains(key)) {
        set_metadata(key, data["info"][key]);
      } else {
        throw BoutException(
            fmt::format("Failed to read metadata '{}' from Amjuel data file at '{:s}'",
                        key, file_path.string()));
      }
    }

    // Extract coeff table
    this->coeffs = data["coeffs"];
  } catch (std::exception& e) {
    throw BoutException(
        "Failed to extract Amjuel reaction data from json file at '{:s}'. Error was {:s}",
        file_path.string(), e.what());
  }
}

///
BoutReal AmjuelData::eval_sigma_vE_nT_impl(BoutReal T, BoutReal n) {
  return eval_amjuel_nT_fit(T, n, this->coeffs);
}

///
BoutReal AmjuelData::eval_sigma_v_ET_impl([[maybe_unused]] BoutReal E,
                                          [[maybe_unused]] BoutReal T) {
  throw BoutException("AmjuelData::eval_sigma_v_ET not yet implemented.");
}

///
BoutReal AmjuelData::eval_sigma_v_nT_impl(BoutReal T, BoutReal n) {
  return eval_amjuel_nT_fit(T, n, this->coeffs);
}

///
BoutReal AmjuelData::eval_sigma_v_T_impl(BoutReal T) {
  return eval_amjuel_T_fit(T, this->coeffs[0]);
}

} // namespace hermes