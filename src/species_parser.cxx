#include "species_parser.hxx"

#include <regex>

#include "bout/utils.hxx"
namespace {
const std::map<std::string, std::string> long_names = {
    {"e", "electron"}, {"h", "hydrogen"}, {"d", "deuterium"},
    {"t", "tritium"},  {"he", "helium"},
};

std::string construct_species_str(std::string base_species, int charge) {
  if (base_species == "e") {
    // Special case for electrons
    if (charge != -1) {
      throw BoutException(
          fmt::format("Unexpected charge for electron species! ({})", charge));
    }
    return "e";
  }

  std::string charge_str("");
  switch (charge) {
  case 0:
    charge_str = "";
    break;
  case 1:
    charge_str = "+";
    break;
  case -1:
    charge_str = "-";
    break;
  default:
    if (charge > 0) {
      charge_str = "+" + std::to_string(charge);
    } else if (charge < 0) {
      charge_str = "-" + std::to_string(-charge);
    }
  }
  return base_species + charge_str;
}

} // namespace

SpeciesTypeAlt SpeciesParser::get_type(const std::string& species_str) {
  return get_type(SpeciesParser(species_str));
}

SpeciesTypeAlt SpeciesParser::get_type(const SpeciesParser& species) {
  if (species.get_element() == "e") {
    return SpeciesTypeAlt::electron;
  } else if (species.is_molecule()) {
    if (species.get_charge() == 0) {
      return SpeciesTypeAlt::molecule;
    } else {
      return SpeciesTypeAlt::molecular_ion;
    }
  } else if (species.get_charge() == 0) {
    return SpeciesTypeAlt::neutral;
  } else {
    return SpeciesTypeAlt::ion;
  }
}

///
SpeciesParser::SpeciesParser(const std::string& species_str) {

  // Extract base species, charge and ionisation state with regex
  // Any number preceding the base species is discarded
  std::regex pattern("^([0-9]*)([a-zA-Z]{1,2}[0-9]*)([\\+|\\-]?)([0-9]*)$");
  std::smatch matches;
  bool has_matches = std::regex_search(species_str, matches, pattern);
  // String must provide at least a base species
  if (!has_matches || !matches[1].matched) {
    throw BoutException(fmt::format(
        "Unable to extract species properties from string '{}'", species_str));
  }

  // Store base species; always lower case
  this->base_species = matches[2];
  std::transform(this->base_species.begin(), this->base_species.end(),
                 this->base_species.begin(), ::tolower);

  std::regex mol_pattern("^([a-zA-Z]{1,2})([0-9]{0,2})$");
  std::smatch mol_matches;
  bool has_mol_matches = std::regex_search(this->base_species, mol_matches, mol_pattern);
  // base_species must at least contain an element
  if (!has_mol_matches || !mol_matches[1].matched) {
    throw BoutException(
        fmt::format("Unable to extract element from species name {}", species_str));
  }
  this->element = mol_matches[1];
  this->_is_molecule = !mol_matches[2].str().empty();

  // Extract charge, electron is a special case
  if (species_str == "e") {
    this->charge = -1;
  } else {
    int sign = matches[3] == "+" ? 1 : matches[3] == "-" ? -1 : 0;
    int num = (matches[4].length() == 0) ? 1 : stringToInt(matches[4]);
    this->charge = sign * num;
  }

  // Stored species string discards any leading number and is always lower case
  this->species_str = construct_species_str(this->base_species, this->charge);
}

///
SpeciesParser::SpeciesParser(const SpeciesParser& other, const int charge)
    : base_species(other.base_species), charge(charge), element(other.element),
      _is_molecule(other._is_molecule),
      species_str(construct_species_str(other.base_species, charge)) {}

///
std::string SpeciesParser::long_name() const {
  auto it = long_names.find(this->base_species);
  if (it != long_names.end()) {
    return it->second;
  } else {
    return this->base_species;
  }
}

///
SpeciesParser SpeciesParser::ionised() {
  if (this->base_species == "e") {
    throw BoutException("Cannot change electron charge!");
  }
  return SpeciesParser(*this, this->charge + 1);
}

///
SpeciesParser SpeciesParser::recombined() {
  if (this->base_species == "e") {
    throw BoutException("Cannot change electron charge!");
  }
  return SpeciesParser(*this, this->charge - 1);
}
