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

std::vector<std::string> ReactionParser::get_species() const {
  return get_str_keys(this->stoich);
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
    // Filter out electrons ('e')
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
    // Filter out electrons ('e')
    std::copy_if(species_names.begin(), species_names.end(),
                 std::back_inserter(filtered_species_names), [](std::string sp_name) {
                   return identifySpeciesType(sp_name) == SpeciesType::ion;
                 });
    break;
  case species_filter::neutral:
    // Filter out electrons ('e')
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
    std::vector<std::string> product_species = get_str_keys(this->products);
    std::sort(product_species.begin(), product_species.end());
    std::sort(species_names.begin(), species_names.end());
    std::set_intersection(species_names.begin(), species_names.end(),
                          product_species.begin(), product_species.end(),
                          std::back_inserter(filtered_species_names));
    break;
  }

  case species_filter::reactants: {
    std::vector<std::string> reactant_species = get_str_keys(this->reactants);
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
