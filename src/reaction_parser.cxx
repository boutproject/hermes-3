#include <algorithm>

#include "hermes_utils.hxx"

#include "reaction_parser.hxx"

ReactionParser::ReactionParser(const std::string& reaction_str)
    : reaction_str(reaction_str) {
  // Assume reactants, products are separated by '->'
  const std::string rp_sep{"->"};
  const std::size_t sep_len = rp_sep.length();
  if (reaction_str.length() < sep_len + 2) {
    throw BoutException("Reaction string not long enough to include at least on reactant "
                        "and one product!");
  }
  auto sep_idx = reaction_str.find(rp_sep);
  if (sep_idx == std::string::npos) {
    throw BoutException("Failed to find '->' separator in the reaction string");
  }
  std::string R = trim(reaction_str.substr(0, sep_idx));
  std::string P = trim(reaction_str.substr(sep_idx + sep_len));

  // Count species in reactants, products
  this->reactants = count_species(R);
  if (this->reactants.size() < 1) {
    throw BoutException("Failed to find any reactants in the reaction string");
  }
  this->products = count_species(P);
  if (this->products.size() < 1) {
    throw BoutException("Failed to find any products in the reaction string");
  }

  diff_reactants_products(this->reactants, this->products);
}

void ReactionParser::diff_reactants_products(const std::map<std::string, int>& R,
                                             const std::map<std::string, int>& P) {
  // Construct standard pop changes map
  this->stoich = std::map<std::string, int>(P);
  for (const auto& [sp_name, pop_change] : R) {
    auto it = P.find(sp_name);
    if (it == P.end()) {
      this->stoich[sp_name] = -pop_change;
    } else {
      this->stoich[sp_name] -= pop_change;
    }
  }

  // If all population changes are zero, mark as "symmetric"
  is_symmetric = std::all_of(this->stoich.begin(), this->stoich.end(),
                             [](const auto& entry) { return entry.second == 0; });

  // Construct momentum-energy pop changes differently if reaction is symmetric
  if (this->is_symmetric) {
    this->mom_energy_stoich.insert(P.begin(), P.end());
    std::transform(
        R.begin(), R.end(),
        std::inserter(this->mom_energy_stoich, this->mom_energy_stoich.end()),
        [](const auto& entry) { return make_pair(entry.first, -entry.second); });
  } else {
    this->mom_energy_stoich.insert(this->stoich.begin(), this->stoich.end());
  }
}

std::vector<std::string> ReactionParser::get_species() const {
  return str_keys(this->stoich);
}

std::vector<std::string> ReactionParser::get_species(species_filter filter) const {
  return get_species(get_species(), filter);
}

std::vector<std::string>
ReactionParser::get_species(std::vector<std::string> species_names,
                            species_filter filter) const {

  std::vector<std::string> filtered_species_names;
  switch (filter) {
  case species_filter::consumed: {
    // -ve population change
    std::copy_if(species_names.begin(), species_names.end(),
                 std::back_inserter(filtered_species_names),
                 [this](std::string sp_name) { return this->stoich.at(sp_name) < 0; });
    break;
  }
  case species_filter::electron:
    // Select only electrons ('e'). Mainly useful to detect cases where no electrons
    // are involved in a reaction.
    std::copy_if(species_names.begin(), species_names.end(),
                 std::back_inserter(filtered_species_names), [](std::string sp_name) {
                   return identifySpeciesType(sp_name) == SpeciesType::electron;
                 });
    break;
  case species_filter::heavy:
    // Filter out electrons ('e')
    std::copy_if(species_names.begin(), species_names.end(),
                 std::back_inserter(filtered_species_names), [](std::string sp_name) {
                   return identifySpeciesType(sp_name) != SpeciesType::electron;
                 });
    break;
  case species_filter::ion:
    std::copy_if(species_names.begin(), species_names.end(),
                 std::back_inserter(filtered_species_names), [](std::string sp_name) {
                   return identifySpeciesType(sp_name) == SpeciesType::ion;
                 });
    break;
  case species_filter::neutral:
    std::copy_if(species_names.begin(), species_names.end(),
                 std::back_inserter(filtered_species_names), [](std::string sp_name) {
                   return identifySpeciesType(sp_name) == SpeciesType::neutral;
                 });
    break;
  case species_filter::produced: {
    // +ve population change
    std::copy_if(species_names.begin(), species_names.end(),
                 std::back_inserter(filtered_species_names),
                 [this](std::string sp_name) { return this->stoich.at(sp_name) > 0; });
    break;
  }
  case species_filter::products: {
    std::vector<std::string> product_species = str_keys(this->products);
    std::sort(product_species.begin(), product_species.end());
    std::sort(species_names.begin(), species_names.end());
    std::set_intersection(species_names.begin(), species_names.end(),
                          product_species.begin(), product_species.end(),
                          std::back_inserter(filtered_species_names));
    break;
  }

  case species_filter::reactants: {
    std::vector<std::string> reactant_species = str_keys(this->reactants);
    std::sort(reactant_species.begin(), reactant_species.end());
    std::sort(species_names.begin(), species_names.end());
    std::set_intersection(species_names.begin(), species_names.end(),
                          reactant_species.begin(), reactant_species.end(),
                          std::back_inserter(filtered_species_names));
    break;
  }
  }
  return filtered_species_names;
}

const std::multimap<std::string, int>&
ReactionParser::get_mom_energy_pop_changes() const {
  return this->mom_energy_stoich;
}

int ReactionParser::pop_change(const std::string sp_name) const {
  return this->stoich.at(sp_name);
}

int ReactionParser::pop_change_product(const std::string sp_name) const {
  if (this->is_symmetric) {
    return this->products.at(sp_name);
  } else {
    return pop_change(sp_name);
  }
}

int ReactionParser::pop_change_reactant(const std::string sp_name) const {
  if (this->is_symmetric) {
    return -1 * this->reactants.at(sp_name);
  } else {
    return pop_change(sp_name);
  }
}
