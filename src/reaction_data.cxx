#include "reaction_data.hxx"

/// ============================= Helper functions ============================
///
std::filesystem::path get_json_db_dir(Options& options) {
  static std::filesystem::path default_json_db_dir =
      std::filesystem::path(__FILE__).parent_path().parent_path() / "json_database";

  std::string json_db_dir =
      options["json_database_dir"]
          .doc("Path to directory containing reaction data json files.")
          .withDefault(default_json_db_dir.string());

  return std::filesystem::path(json_db_dir);
}

/// ============================== ReactionData ===============================

///
BoutReal ReactionData::eval_sigma_vE_nT(BoutReal T, BoutReal n) {
  if (this->fit_type != RateParamsTypes::nT) {
    throw BoutException(
        fmt::format("Trying to call eval_sigma_vE_nT, but reaction data fit type is {}",
                    toString(this->fit_type)));
  }
  return this->eval_sigma_vE_nT_impl(T, n);
}

///
BoutReal ReactionData::eval_sigma_v_ET(BoutReal E, BoutReal T) {
  if (this->fit_type != RateParamsTypes::ET) {
    throw BoutException(
        fmt::format("Trying to call eval_sigma_v_ET, but reaction data fit type is {}",
                    toString(this->fit_type)));
  }
  return eval_sigma_v_ET_impl(E, T);
}

///
BoutReal ReactionData::eval_sigma_v_nT(BoutReal T, BoutReal n) {
  if (this->fit_type != RateParamsTypes::nT) {
    throw BoutException(
        fmt::format("Trying to call eval_sigma_v_nT, but reaction data fit type is {}",
                    toString(this->fit_type)));
  }
  return eval_sigma_v_nT_impl(T, n);
}

///
BoutReal ReactionData::eval_sigma_v_T(BoutReal T) {
  if (this->fit_type != RateParamsTypes::T) {
    throw BoutException(
        fmt::format("Trying to call eval_sigma_v_T, but reaction data fit type is {}",
                    toString(this->fit_type)));
  }
  return eval_sigma_v_T_impl(T);
}

///
RateParamsTypes ReactionData::get_fit_type() const { return this->fit_type; }

///
BoutReal ReactionData::get_metadata(const std::string& key) const {
  if (this->metadata.find(key) == this->metadata.end()) {
    throw BoutException(fmt::format(
        "Trying to access reaction metadata '{}' that hasn't been set.", key));
  }
  return this->metadata.at(key);
}

///
std::string ReactionData::src_str() {
  return fmt::format("{}_{}", toString(this->type), this->data_label);
}

/// ========================= ReactionDataWithCoeffs ==========================
///
ReactionDataWithCoeffs::ReactionDataWithCoeffs(ReactionDataTypes data_type,
                                               std::string data_label,
                                               std::vector<std::string> metadata_keys)
    : ReactionData(data_type, std::move(data_label), std::move(metadata_keys)), coeffs() {
}

///
const std::vector<std::vector<BoutReal>>& ReactionDataWithCoeffs::get_coeffs() const {
  return this->coeffs;
}

constexpr decltype(ReactionDataFactory::type_name) ReactionDataFactory::type_name;
constexpr decltype(ReactionDataFactory::section_name) ReactionDataFactory::section_name;
constexpr decltype(ReactionDataFactory::option_name) ReactionDataFactory::option_name;
constexpr decltype(ReactionDataFactory::default_type) ReactionDataFactory::default_type;