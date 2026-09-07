#pragma once
#ifndef SPECIES_PARSER_H
#define SPECIES_PARSER_H

#include <string>

#include "bout/bout_enum_class.hxx"

BOUT_ENUM_CLASS(SpeciesTypeAlt, electron, ion, molecular_ion, neutral, molecule);

/**
 * @brief Class to parse species strings and extract charge information.
 *
 */
class SpeciesParser {

public:
  static SpeciesTypeAlt get_type(const std::string& species_str);
  static SpeciesTypeAlt get_type(const SpeciesParser& species);
  /**
   * @brief Construct a new SpeciesParser object by extracting the base species and charge
   * from a species string.
   * @details Valid string requirements:
   * - Base species (atom, molecule, or electron) consists of 1 or 2 letters optionally
   * followed by an integer, and is stored in lower case. e.g. "he" or "h2"
   * - Any integer can *precede* the base species, but is discarded
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
  std::string get_base_species() const { return this->base_species; }
  std::string get_str() const { return this->species_str; }

  /**
   * @brief Get an appropriate long name for the base species (e.g. "hydrogen" for "h").
   *
   * @return std::string The long name of the species, or the base species name if no long
   * name is defined
   */
  std::string long_name() const;

  /**
   * @brief Construct a new object with the charge increased by 1.
   *
   * @return SpeciesParser Ionised species object
   */
  SpeciesParser ionised();

  std::string get_element() const { return this->element; }

  bool is_molecule() const { return this->_is_molecule; }

  bool is_neutral() const { return this->charge == 0; }

  /**
   * @brief Construct a new object with the charge reduced by 1.
   *
   * @return SpeciesParser Recombined species object
   */
  SpeciesParser recombined();

private:
  /// Almost-copy constructor that allows charge to be changed. Used by ionised() and recombined().
  SpeciesParser(const SpeciesParser& other, const int charge);

  /// Base species. Can be an atom, a molecule (made up of a single element), or an electron
  std::string base_species;

  /// Charge of the species
  int charge;

  /// Element symbol (e.g. "H" for hydrogen)
  std::string element;

  /// Flag indicating if the species is a molecule
  bool _is_molecule;

  /// Species string
  std::string species_str;
};

#endif // SPECIES_PARSER_H
