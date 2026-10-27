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

SpeciesType SpeciesParser::get_type(const std::string& species_str) {
  return SpeciesParser(species_str).get_type();
}

/**
 * @brief Map the parsed species to a SpeciesType.
 *
 * @details Recognise special-case strings like "i" and "e"; othewise use the
 * charge and molecular status extracted from the species string to determine type.
 *
 * @return SpeciesType
 */
SpeciesType SpeciesParser::get_type() const {
  // (Recognise 'e' as an electron)
  if (get_element() == "e") {
    return SpeciesType::electron;
  }

  // (Recognise legacy generic ion name, which has no explicit charge suffix)
  if (get_str() == "i") {
    return SpeciesType::ion;
  }
  // Molecules and molecular ions
  if (is_molecule()) {
    if (get_charge() == 0) {
      return SpeciesType::molecule;
    }
    return SpeciesType::molecular_ion;
  }
  // Monatomic neutral species => atom
  if (get_charge() == 0) {
    return SpeciesType::neutral;
  }
  // Monatomic species with non-zero charge => ion
  return SpeciesType::ion;
}

///
SpeciesParser::SpeciesParser(const std::string& species_str) {

  // Capture 1. base species (element with optional molecule count), 2. charge sign and
  // 3. charge magnitude. Any leading number is discarded in a non-capturing group.
  static const std::regex pattern(
      "(?:[1-9][0-9]*)?(([a-zA-Z]{1,2})([1-9][0-9]*)?)([+-]?)([1-9][0-9]*)?");
  std::smatch all_components;
  if (!std::regex_match(species_str, all_components, pattern)) {
    throw BoutException(fmt::format(
        "Unable to extract species properties from string '{}'", species_str));
  }

  // Store base species (atom or molecule) in lower case
  this->base_species = all_components[1];
  std::transform(this->base_species.begin(), this->base_species.end(),
                 this->base_species.begin(), ::tolower);
  // Store element, also in lower case
  this->element = this->base_species.substr(
      0, static_cast<std::string::size_type>(all_components[2].length()));
  this->_is_molecule = all_components[3].matched;

  // Store charge. Electron is a special case
  if (species_str == "e") {
    this->charge = -1;
  } else {
    int sign = all_components[4] == "+" ? 1 : all_components[4] == "-" ? -1 : 0;
    int num = (all_components[5].length() == 0) ? 1 : stringToInt(all_components[5]);
    this->charge = sign * num;
  }

  // The stored species string is a sanitized version of the input;
  // Any leading number is discarded and the string is always lower case
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
