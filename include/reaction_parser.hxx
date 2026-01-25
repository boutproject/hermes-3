#pragma once
#ifndef REACTION_PARSER_H
#define REACTION_PARSER_H

#include <algorithm>
#include <regex>

#include "bout/utils.hxx"

/// An enum class with which to identify various useful species subsets
enum class species_filter {
  consumed,
  electron,
  ion,
  neutral,
  reactants,
  produced,
  products,
  heavy
};

static inline std::map<std::string, int> count_species(std::string expr) {
  // std::regex_iterator instead?
  std::map<std::string, int> counts;
  const std::regex pattern("([0-9]*)([a-zA-Z]*[0-9]*\\+?\\-?[0-9]*)");
  for (auto el : strsplit(expr, ' ')) {
    if (el.compare("+") != 0) {
      std::smatch matches;
      [[maybe_unused]] const bool has_matches = std::regex_search(el, matches, pattern);
      ASSERT1(has_matches);

      int el_count = (matches[1].length() == 0) ? 1 : stringToInt(matches[1]);
      std::string name = matches[2];
      auto it = counts.find(name);
      if (it == counts.end()) {
        counts[name] = el_count;
      } else {
        counts[name] += el_count;
      }
    }
  }
  return counts;
}

/**
 * @brief Util function to get the keys of a std::string => T map
 *
 * @tparam T the type of the map values
 * @param map the map
 * @return std::vector<std::string> vector of keys
 */
template <typename T>
static inline std::vector<std::string> get_str_keys(const std::map<std::string, T>& map) {
  std::vector<std::string> keys;
  std::transform(map.begin(), map.end(), std::back_inserter(keys),
                 [](const std::pair<std::string, T>& pair) { return pair.first; });
  return keys;
}

/**
 * @brief A class to parse reaction strings and extract the stoichiometric vector (net
 * population changes for each species). Also used to retrieve the names of species in
 * various useful subsets (reactants, products, non-electron species, etc.)
 *
 */
class ReactionParser {
public:
  ReactionParser(const std::string& reaction_str);

  /// Public getter for underlying reaction string
  const std::string get_reaction_str() const { return this->reaction_str; }

  /// Get the stoichiometric vector for this reaction (as a species_name=>pop_change map)
  const std::map<std::string, int>& get_stoich() { return this->stoich; }

  /**
   * @brief Get the names of all species identified by the parser
   *
   * @return std::vector<std::string> the list of species names
   */
  std::vector<std::string> get_species() const;

  /**
   * @brief Apply a filter to the list of species identified by the parser.
   *
   * @param filter the filter to apply
   * @return std::vector<std::string> the filterered list of species names
   */
  std::vector<std::string> get_species(species_filter filter) const;

  /**
   * @brief Apply a filter to a list of species names.
   *
   * @param species_names the list of species names to filter
   * @param filter the filter to apply
   * @return std::vector<std::string> the filtered list of names
   */
  std::vector<std::string> get_species(std::vector<std::string> species_names,
                                       species_filter filter) const;

  /**
   * @brief Apply multiple filters to a list of species names.
   *
   * @details Variadic so that it can be applied recursively.
   *
   * @tparam FilterTypes
   * @param species_names the list of species names to filter
   * @param first_filter the first filter
   * @param other_filters other filters
   * @return std::vector<std::string> the filtered list of names
   */
  template <typename... FilterTypes>
  std::vector<std::string> get_species(std::vector<std::string> species_names,
                                       species_filter first_filter,
                                       FilterTypes... other_filters) const {
    std::vector<std::string> first_filter_applied =
        get_species(species_names, first_filter);
    return get_species(first_filter_applied, other_filters...);
  }

  /**
   * @brief Apply multiple filters to the list of species identified by the parser
   *
   * @details Variadic so that it can be applied recursively.
   *
   * @tparam FilterTypes
   * @param species_names the list of species names to filter
   * @param first_filter the first filter
   * @param other_filters other filters
   * @return std::vector<std::string> the filtered list of names
   */
  template <typename... FilterTypes>
  std::vector<std::string> get_species(FilterTypes... filters) const {
    return get_species(get_species(), filters...);
  }

private:
  /// The reaction string
  const std::string reaction_str;
  /// Map of species name => population change for reactants
  std::map<std::string, int> reactants;
  /// Map of species name => population change for products
  std::map<std::string, int> products;
  /// Stoichiometric 'vector' (map of species name => population change)
  std::map<std::string, int> stoich;

  /**
   * @brief Util function to compute the stoichiometric 'vector' (map) by taking the
   * difference between the reactant and product population changes.
   *
   * @param R the reactant population changes
   * @param P the product population changes
   */
  void diff_reactants_products(const std::map<std::string, int>& R,
                               const std::map<std::string, int>& P) {
    this->stoich = std::map<std::string, int>(P);
    for (const auto& [sp_name, pop_change] : R) {
      auto it = P.find(sp_name);
      if (it == P.end()) {
        stoich[sp_name] = -pop_change;
      } else {
        stoich[sp_name] -= pop_change;
      }
    }
  }
};

#endif
