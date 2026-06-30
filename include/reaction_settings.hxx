#ifndef REACTION_SETTINGS_HXX
#define REACTION_SETTINGS_HXX

#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include <bout/boutexception.hxx>
#include <fmt/format.h>

/**
 * @brief Helper functions for extracting and processing reaction settings from
 * the options
 *
 */

/// @brief String returned if no default data id is found
static const std::string NO_DATA_ID_DEFAULT_FOUND = "__NO_DEFAULT_ID_FOUND__";

/**
 * @brief Map of default data ids for each reaction data type and reaction string.
 * - First key is the reaction data type
 * - Second key is the reaction string
 * - Third key is the reaction coefficient type
 * - Value is the default data id.
 */
static std::map<ReactionDataTypes,
                std::map<std::string, std::map<ReactionCoeffTypes, std::string>>>
    default_data_ids;

/**
 * @brief Get the default data id for a given reaction data type and reaction string,
 * if one exists. Otherwise return a fixed string indicating no default was found.
 * @param reaction_str The reaction string (e.g. "H + e -> H+ + 2e")
 * @param data_type reaction data type enum (e.g. ReactionDataTypes::Amjuel)
 * @param coeff_type reaction coefficient type enum (e.g. ReactionCoeffTypes::sigma_v)
 * @return std::string The corresponding default data id (e.g. "H.2_3.1.8"), or a
 * fixed string if no default is found.
 */
std::string get_default_data_id(const std::string& reaction_str,
                                const ReactionDataTypes& data_type,
                                const ReactionCoeffTypes& coeff_type);

static inline void add_default_id(ReactionDataTypes data_type,
                                  const std::string& reaction_str,
                                  ReactionCoeffTypes coeff_type,
                                  const std::string& data_id) {
  default_data_ids[data_type][reaction_str][coeff_type] = data_id;
}

/**
 * @brief Populate a map of default data ids for some of the most frequently used
 * reactions. Removes the need for users to specify data ids in the .inp file in most
 * cases.
 *
 */
static inline void generate_default_data_ids_map() {
  constexpr char H_isotopes[] = {'h', 'd', 't'};

  // Various atomic and molecular reactions involving H isotopes
  std::string reaction_str;
  for (const char& isotope : H_isotopes) {
    reaction_str = fmt::format("{} + e -> {}+ + 2e", isotope, isotope);
    // Ionisation of atoms
    add_default_id(ReactionDataTypes::Amjuel, reaction_str, ReactionCoeffTypes::sigma_v,
                   "H.4_2.1.5");
    add_default_id(ReactionDataTypes::Amjuel, reaction_str, ReactionCoeffTypes::sigma_v_E,
                   "H.10_2.1.5");

    // Recombination of ions
    reaction_str = fmt::format("{}+ + e -> {}", isotope, isotope);
    add_default_id(ReactionDataTypes::Amjuel, reaction_str, ReactionCoeffTypes::sigma_v,
                   "H.4_2.1.8");
    add_default_id(ReactionDataTypes::Amjuel, reaction_str, ReactionCoeffTypes::sigma_v_E,
                   "H.10_2.1.8");

    // Dissociation of neutral molecules
    reaction_str = fmt::format("{}2 + e -> 2{} + e", isotope, isotope);
    add_default_id(ReactionDataTypes::Amjuel, reaction_str, ReactionCoeffTypes::sigma_v,
                   "H.4_2.2.5g");

    // Non-dissociative ionisation of neutral molecules
    reaction_str = fmt::format("{}2 + e -> {}2+ + 2e", isotope, isotope);
    add_default_id(ReactionDataTypes::Amjuel, reaction_str, ReactionCoeffTypes::sigma_v,
                   "H.4_2.2.9");

    // Dissociative ionisation of neutral molecules
    reaction_str = fmt::format("{}2 + e -> {} + {}+ + 2e", isotope, isotope, isotope);
    add_default_id(ReactionDataTypes::Amjuel, reaction_str, ReactionCoeffTypes::sigma_v,
                   "H.4_2.2.10");

    // Dissociative ionisation of (singly-) ionised molecules
    reaction_str = fmt::format("{}2+ + e -> 2{}+ + 2e", isotope, isotope, isotope);
    add_default_id(ReactionDataTypes::Amjuel, reaction_str, ReactionCoeffTypes::sigma_v,
                   "H.4_2.2.11");

    // Dissociative excitation of (singly-) ionised molecules
    reaction_str = fmt::format("{}2+ + e -> {} + {}+ + e", isotope, isotope, isotope);
    add_default_id(ReactionDataTypes::Amjuel, reaction_str, ReactionCoeffTypes::sigma_v,
                   "H.4_2.2.12");

    // Dissociative recombination of (singly-) ionised molecules
    reaction_str = fmt::format("{}2+ + e -> 2{}", isotope, isotope);
    add_default_id(ReactionDataTypes::Amjuel, reaction_str, ReactionCoeffTypes::sigma_v,
                   "H.4_2.2.14");
  }

  // CX and elastic collisions involving H isotopes
  for (const char& isotope1 : H_isotopes) {
    for (const char& isotope2 : H_isotopes) {
      // Atom-ion CX
      std::string reaction_str =
          fmt::format("{} + {}+ -> {}+ + {}", isotope1, isotope2, isotope1, isotope2);
      add_default_id(ReactionDataTypes::Amjuel, reaction_str, ReactionCoeffTypes::sigma_v,
                     "H.2_3.1.8");
      // Molecule-ion CX
      reaction_str =
          fmt::format("{}2 + {}+ -> {}2+ + {}", isotope1, isotope2, isotope1, isotope2);
      add_default_id(ReactionDataTypes::Amjuel, reaction_str, ReactionCoeffTypes::sigma_v,
                     "H.2_3.2.3");
      // Molecule-ion elastic collisions
      reaction_str =
          fmt::format("{}2 + {}+ -> {}2 + {}+", isotope1, isotope2, isotope1, isotope2);
      add_default_id(ReactionDataTypes::Amjuel, reaction_str, ReactionCoeffTypes::sigma_v,
                     "H.2_0.3T");
    }
  }

  // He izn
  add_default_id(ReactionDataTypes::Amjuel, "he + e -> he+ + 2e",
                 ReactionCoeffTypes::sigma_v, "H.4_2.3.9a");
  add_default_id(ReactionDataTypes::Amjuel, "he + e -> he+ + 2e",
                 ReactionCoeffTypes::sigma_v_E, "H.10_2.3.9a");
  // He rec
  add_default_id(ReactionDataTypes::Amjuel, "he+ + e -> he", ReactionCoeffTypes::sigma_v,
                 "H.4_2.3.13a");
  add_default_id(ReactionDataTypes::Amjuel, "he+ + e -> he",
                 ReactionCoeffTypes::sigma_v_E, "H.10_2.3.13a");
}

/**
 * @brief Get the default data id for a given reaction data type and reaction string,
 * if it exists. Otherwise return a fixed string indicating no default was found.
 * @param reaction_str The reaction string (e.g. "H + e -> H+ + 2e")
 * @param data_type reaction data type enum (e.g. ReactionDataTypes::Amjuel)
 * @param coeff_type reaction coefficient type enum (e.g. ReactionCoeffTypes::sigma_v)
 * @return std::string The corresponding default data id (e.g. "H.2_3.1.8"), or a
 * fixed string if no default is found.
 */
inline std::string get_default_data_id(const std::string& reaction_str,
                                       const ReactionDataTypes& data_type,
                                       const ReactionCoeffTypes& coeff_type) {
  if (default_data_ids.empty()) {
    generate_default_data_ids_map();
  }
  auto data_type_it = default_data_ids.find(data_type);
  if (data_type_it != default_data_ids.end()) {
    auto reaction_it = data_type_it->second.find(reaction_str);
    if (reaction_it != data_type_it->second.end()) {
      auto coeff_it = reaction_it->second.find(coeff_type);
      if (coeff_it != reaction_it->second.end()) {
        return coeff_it->second;
      }
    }
  }
  // Return a fixed string if no default id was found
  return NO_DATA_ID_DEFAULT_FOUND;
}

static inline std::string
get_default_data_ids_str(const std::vector<std::string> reaction_strs,
                         const std::vector<ReactionDataTypes> data_types,
                         const ReactionCoeffTypes coeff_type) {
  std::stringstream default_ids_str;
  default_ids_str << get_default_data_id(reaction_strs[0], data_types[0], coeff_type);
  for (std::size_t i = 1; i < reaction_strs.size(); i++) {
    default_ids_str << ","
                    << get_default_data_id(reaction_strs[i], data_types[i], coeff_type);
  }
  return default_ids_str.str();
}

/**
 * @brief Split comma-separated strings of the form (a,b,c,d), removing parentheses and
 * whitespace.
 *
 * @param csv_str input string
 * @return std::vector<std::string>
 */
inline std::vector<std::string> split_csv_str(std::string csv_str) {
  std::vector<std::string> str_vec;
  std::stringstream ss(csv_str);
  std::string item;

  while (std::getline(ss, item, ',')) {
    // Remove parentheses
    std::regex match_parentheses("\\(|\\)");
    item = std::regex_replace(item, match_parentheses, "");
    // Remove whitespace
    item.erase(0, item.find_first_not_of(" \t\n\r"));
    item.erase(item.find_last_not_of(" \t\n\r") + 1);

    if (!item.empty()) {
      str_vec.push_back(item);
    }
  }
  return str_vec;
}

/**
 * @brief Split comma-separated strings of the form (a,b,c,d), removing parentheses and
 * whitespace. If the result has \p num_expected values, return, if it has only one value,
 * duplicate the value \p num_expected times. If it's any other size, throw an exception.
 *
 * @param csv_str the string to split
 * @param num_expected the expected number of values
 * @param lbl label used to refer to the string if there's an error
 * @return std::vector<std::string>
 */
inline std::vector<std::string> split_csv_str(std::string csv_str,
                                              const std::size_t num_expected,
                                              const std::string lbl) {
  std::vector<std::string> str_vec = split_csv_str(csv_str);
  const std::size_t num_set = str_vec.size();
  if (num_set == num_expected) {
    // Already of expected size
  } else if (num_set == 1) {
    // Duplicate single value to make result have size num_expected
    str_vec.resize(num_expected, str_vec[0]);
  } else {
    throw BoutException(fmt::format("Expected either one {:s}, or one "
                                    "per reaction string ({:d}). Got {:d}! ",
                                    lbl, num_expected, num_set));
  }
  return str_vec;
}

#endif // REACTION_SETTINGS_HXX
