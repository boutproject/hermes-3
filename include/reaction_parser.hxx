#pragma once
#ifndef REACTION_PARSER_H
#define REACTION_PARSER_H

#include <algorithm>
#include <regex>

#include "bout/utils.hxx"

static inline std::map<std::string, int> count_species(std::string expr) {
  // std::regex_iterator instead?
  std::map<std::string, int> counts;
  std::regex pattern("([0-9]*)([a-zA-Z]*[0-9]*\\+?\\-?[0-9]*)");
  for (auto el : strsplit(expr, ' ')) {
    if (el.compare("+") != 0) {
      std::smatch matches;
      bool has_matches = std::regex_search(el, matches, pattern);
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

template <typename T>
static inline std::vector<std::string> get_str_keys(const std::map<std::string, T>& map) {
  std::vector<std::string> keys;
  std::transform(map.begin(), map.end(), std::back_inserter(keys),
                 [](const std::pair<std::string, T>& pair) { return pair.first; });
  return keys;
}

class ReactionParser {
public:
  ReactionParser(const std::string& reaction_str) : reaction_str(reaction_str) {
    // Assume reactants, products are separated by '->'
    const std::string rp_sep{"->"};
    const std::size_t sep_len = rp_sep.length();
    ASSERT1(reaction_str.length() >= sep_len + 2);
    auto sep_idx = reaction_str.find(rp_sep);
    ASSERT1(sep_idx > sep_len);
    std::string R = trim(reaction_str.substr(0, sep_idx));
    std::string P = trim(reaction_str.substr(sep_idx + sep_len));

    // Count species in reactants, products
    this->reactants = count_species(R);
    this->products = count_species(P);

    diff_reactants_products(this->reactants, this->products);
  }

  const std::map<std::string, int>& get_stoich() { return this->stoich; }

  std::vector<std::string> get_product_species() { return get_str_keys(this->products); }
  std::vector<std::string> get_reactant_species() {
    return get_str_keys(this->reactants);
  }

private:
  const std::string reaction_str;
  std::map<std::string, int> reactants;
  std::map<std::string, int> products;
  std::map<std::string, int> stoich;

  void diff_reactants_products(const std::map<std::string, int>& R,
                               const std::map<std::string, int>& P) {
    this->stoich = std::map<std::string, int>(P);
    for (auto R_el : R) {
      auto it = P.find(R_el.first);
      if (it == P.end()) {
        stoich[R_el.first] = -R_el.second;
      } else {
        stoich[R_el.first] -= R_el.second;
      }
    }
  }
};

#endif