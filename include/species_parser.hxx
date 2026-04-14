#pragma once
#ifndef SPECIES_PARSER_H
#define SPECIES_PARSER_H

#include <string>

/**
 * @brief Class to parse species strings and extract charge information.
 *
 */
class SpeciesParser {
public:
  /**
   * @brief Construct a new SpeciesParser object by extracting the element and charge from
   * a species string.
   * @details Valid string requirements:
   * - Element name must be 1 or 2 letters and is stored in lower case.
   * - Any integer can precede the element name, but is discarded
   *     e.g. "2he" is parsed as "he"
   * - A "+" or a "-", optionally followed by an integer is interpreted as the charge.
   *     e.g. "he-1" (-1), "H" (0) or "ne+8" (+8)
   * - Electrons are a special case and can be specified as "e" or "e-".
   *
   * @param species_str The species string to parse
   * @throws BoutException if the string cannot be parsed or contains invalid information
   */
  SpeciesParser(const std::string& species_str);

  int get_charge() const { return this->charge; }
  std::string get_element() const { return this->element; }
  std::string get_str() const { return this->species_str; }

  /**
   * @brief Get an appropriate long name for the element in this species (e.g. "hydrogen"
   * for "h").
   *
   * @return std::string The long name of the species, or the element name if no long name
   * is defined
   */
  std::string long_name() const;

  /**
   * @brief Construct a new object with the charge increased by 1.
   *
   * @return SpeciesParser Ionised species object
   */
  SpeciesParser ionised();

  /**
   * @brief Construct a new object with the charge reduced by 1.
   *
   * @return SpeciesParser Recombined species object
   */
  SpeciesParser recombined();

private:
  SpeciesParser(const std::string element, int charge);

  /// Species string
  std::string species_str;

  /// Atomic element of the species (or 'e' for electrons)
  std::string element;

  /// Charge of the species
  int charge;
};

#endif // SPECIES_PARSER_H
