#pragma once
#ifndef REACTION_PARSER_H
#define REACTION_PARSER_H

#include <algorithm>
#include <regex>

#include "hermes_utils.hxx"
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

  /**
   * @brief Apply one or more filters to a list of species names and return the one and
   * only match.
   *
   * @details Variadic so that it can be applied recursively.
   * @throw BoutException if there isn't exactly one match.
   *
   * @tparam FilterTypes
   * @param species_names the list of species names to filter
   * @param first_filter the first filter
   * @param other_filters zero or more other filters
   * @return std::string the matching species name
   */
  template <typename... FilterTypes>
  std::string get_single_species(std::vector<std::string> species_names,
                                 species_filter first_filter,
                                 FilterTypes... other_filters) const {
    std::vector<std::string> matches =
        get_species(species_names, first_filter, other_filters...);
    ASSERT0(matches.size() == 1);
    return matches[0];
  }

  /**
   * @brief Apply one or more filters to the list of species identified by the parser and
   * return the one and only match.
   * @throws BoutException if there is more than one match.
   *
   * @tparam FilterTypes
   * @param filters one or more instances of species_filter
   * @return std::string the matching species name
   */
  template <typename... FilterTypes>
  std::string get_single_species(FilterTypes... filters) const {
    std::vector<std::string> matches = get_species(filters...);
    ASSERT0(matches.size() == 1);
    return matches[0];
  }

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
   * @param filters one or more instances of species_filter
   * @return std::vector<std::string> the filtered list of names
   */
  template <typename... FilterTypes>
  std::vector<std::string> get_species(FilterTypes... filters) const {
    return get_species(get_species(), filters...);
  }

  /**
   * @brief Return a map of population changes to use when computing momentum and energy
   * transfer. For non-symmetric reactions, this is the normal stoichiometry vector; for
   * symmetric reactions, it returns the 'split' version, where species names are
   * repeated, with -ve values for reactants and +ve values for products.
   *
   * @return const std::map<std::string, int>&
   */
  const std::multimap<std::string, int>& get_mom_energy_pop_changes() const;

  /**
   * @brief Get the overall population change of a species.
   *
   * @param sp_name the species name
   * @return int the population change
   */
  int pop_change(const std::string sp_name) const;

  /**
   * @brief Get the population change of a product species.
   * If the left and right sides of the reaction string are the same, return the change on
   * the PRODUCT SIDE ONLY, otherwise return the usual net population change.
   *
   * @param sp_name the species name
   * @return int the population change
   */
  int pop_change_product(const std::string sp_name) const;

  /**
   * @brief Get the population change of a reactant species.
   * If the left and right sides of the reaction string are the same, return the change on
   * the REACTANT SIDE ONLY, otherwise return the usual net population change.
   *
   * @param sp_name the species name
   * @return int the population change
   */
  int pop_change_reactant(const std::string sp_name) const;

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
   * 'Split' stoichiometric values. useful when dealing with symmetric reactions.
   *   Species names appear twice, once with the -ve (reactant)
   *   pop. change, once with the +ve (product) pop. change.
   */
  std::multimap<std::string, int> mom_energy_stoich;

  /// Flag to identify reactions where LHS == RHS (e.g. symmetric CX)
  bool is_symmetric;

  /**
   * @brief Util function to compute the stoichiometric 'vector' (map) by taking the
   * difference between the reactant and product population changes. Also computes a
   * separate version used for momentum and energy sources. This second map differs from
   * the first if the reaction is symmetric, in which case each species appear twice, once
   * with the (-ve) reactant pop. change and once with the (+ve) product pop. change.
   *
   * @param R the reactant population changes
   * @param P the product population changes
   */
  void diff_reactants_products(const std::map<std::string, int>& R,
                               const std::map<std::string, int>& P);
};

#endif
