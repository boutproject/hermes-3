#include "species_parser.hxx"

#include <regex>

#include "bout/utils.hxx"

namespace {
static std::map<std::string, std::string> long_names = {
    {"e", "electron"}, {"h", "hydrogen"}, {"d", "deuterium"},
    {"t", "tritium"},  {"he", "helium"},
};

std::string construct_species_str(std::string element, int charge) {
  if (element == "e") {
    // Special case for electrons
    if (charge != -1) {
      throw BoutException(
          fmt::format("Unexpected charge for electron species! ({})", charge));
    }
    return "e";
  } else {

    std::string charge_str;
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
    return element + charge_str;
  }
}

} // namespace

///
SpeciesParser::SpeciesParser(const std::string& species_str) {

  // Extract element name, charge and ionisation state with regex
  // Any number preceding the element is discarded
  std::regex pattern("^([0-9]*)([a-zA-Z]{1,2})([\\+|\\-]?)([0-9]*)$");
  std::smatch matches;
  bool has_matches = std::regex_search(species_str, matches, pattern);
  // String must provide at least an element name
  if (!has_matches || !matches[1].matched) {
    throw BoutException(
        fmt::format("Unable to extract charge from species name {}", species_str));
  }

  // Store element name; always lower case
  this->element = matches[2];
  if (this->element.empty()) {
    throw BoutException(
        fmt::format("Unable to extract element name from species name {}", species_str));
  }
  std::transform(this->element.begin(), this->element.end(), this->element.begin(),
                 ::tolower);

  // Extract charge, electron is a special case
  if (species_str == "e") {
    this->charge = -1;
  } else {
    int sign = matches[3] == "+" ? 1 : matches[3] == "-" ? -1 : 0;
    int num = (matches[4].length() == 0) ? 1 : stringToInt(matches[4]);
    this->charge = sign * num;
  }

  // Stored species string discards any leading number and is always lower case
  this->species_str = construct_species_str(this->element, this->charge);
}

///
SpeciesParser::SpeciesParser(const std::string element, int charge)
    : element(element), charge(charge) {
  this->species_str = construct_species_str(element, charge);
}

///
std::string SpeciesParser::long_name() const {
  auto it = long_names.find(this->element);
  if (it != long_names.end()) {
    return it->second;
  } else {
    return this->element;
  }
}

///
SpeciesParser SpeciesParser::ionised() {
  if (this->element == "e") {
    throw BoutException("Cannot change electron charge!");
  }
  int new_charge = this->charge + 1;
  return SpeciesParser(this->element, new_charge);
}

///
SpeciesParser SpeciesParser::recombined() {
  if (this->element == "e") {
    throw BoutException("Cannot change electron charge!");
  }
  int new_charge = this->charge - 1;
  return SpeciesParser(this->element, new_charge);
}
