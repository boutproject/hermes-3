#pragma once
#ifndef ADAS_LITHIUM_H
#define ADAS_LITHIUM_H

#include "adas_reaction.hxx"

#include <array>
#include <initializer_list>
#include <string_view>

/// Ionisation energies in eV
/// from https://www.webelements.com/lithium/atoms.html
/// Conversion 1 kJ mol‑1 = 1.0364e-2 eV
/// These are added (removed) from the electron energy during recombination (ionisation)
constexpr std::array<BoutReal, 3> lithium_ionisation_energy{5.39, 75.64, 122.45};

/// The name of the species. This initializer list can be passed to a string
/// constructor, or used to index into an Options tree.
///
/// li, li+, li+2, ...
///
/// Special cases for level=0, 1 and 3
///
/// @tparam level  The ionisation level: 0 is neutral, 3 is fully stripped.
template <int level>
constexpr char lithium_species_name[5] = {'l', 'i', '+', char('0' + level), '\0'};

template <>
constexpr char lithium_species_name<3>[5] = {'l', 'i', '+', '3', '\0'};

template <>
constexpr char lithium_species_name<1>[5] = {'l', 'i', '+', '\0'};

template <>
constexpr char lithium_species_name<0>[5] = {'l', 'i', '\0'};

template <int level>
struct lithium_species {
  static constexpr auto value = lithium_species_name<level>;
};

/// ADAS effective ionisation (ADF11)
///
/// @tparam level  The ionisation level of the ion on the left of the reaction
template <int level>
struct ADASLithiumIonisation : public OpenADAS {
  ADASLithiumIonisation(const ADASLithiumIonisation&) = delete;
  ADASLithiumIonisation(ADASLithiumIonisation&&) = delete;
  ADASLithiumIonisation& operator=(const ADASLithiumIonisation&) = delete;
  ADASLithiumIonisation& operator=(ADASLithiumIonisation&&) = delete;
  ADASLithiumIonisation(std::string name, Options& alloptions, Solver*)
      : OpenADAS(name, alloptions["units"], "scd96_li.json", "plt96_li.json",
                 lithium_species_name<level>, lithium_species_name<level + 1>, level,
                 -lithium_ionisation_energy[level]) {}
  static constexpr auto type = izn_component_name_v<lithium_species, level>;

  std::string typeName() const final { return std::string(type); }

private:
  void transform_impl(GuardedOptions& state) override {
    calculate_rates(
        state["species"]["e"],                            // Electrons
        state["species"][lithium_species_name<level>],    // From this ionisation state
        state["species"][lithium_species_name<level + 1>] // To this state
    );
  }
};

/////////////////////////////////////////////////

/// ADAS effective recombination coefficients (ADF11)
///
/// @tparam level  The ionisation level of the ion on the right of the reaction
template <int level>
struct ADASLithiumRecombination : public OpenADAS {
  /// @param alloptions  The top-level options. Only uses the ["units"] subsection.
  ADASLithiumRecombination(std::string name, Options& alloptions, Solver*)
      : OpenADAS(std::move(name), alloptions["units"], "acd96_li.json", "prb96_li.json",
                 lithium_species_name<level + 1>, lithium_species_name<level>, level,
                 lithium_ionisation_energy[level]) {}

  static constexpr auto type = rec_component_name_v<lithium_species, level>;

  std::string typeName() const final { return std::string(type); }

private:
  void transform_impl(GuardedOptions& state) override {
    calculate_rates(
        state["species"]["e"],                             // Electrons
        state["species"][lithium_species_name<level + 1>], // From this ionisation state
        state["species"][lithium_species_name<level>]      // To this state
    );
  }
};

/// @tparam level     The ionisation level of the ion on the right of the reaction
/// @tparam Hisotope  The hydrogen isotope ('h', 'd' or 't')
template <int level, char Hisotope>
struct ADASLithiumCX : public OpenADASChargeExchange {
  /// @param alloptions  The top-level options. Only uses the ["units"] subsection.
  ADASLithiumCX(std::string name, Options& alloptions, Solver*)
      : OpenADASChargeExchange(std::move(name), alloptions["units"], "ccd89_li.json",
                               lithium_species_name<level + 1>, {Hisotope},
                               lithium_species_name<level>, {Hisotope, '+'}, level) {}

  static constexpr auto type = cx_component_name_v<lithium_species, level, Hisotope>;

  std::string typeName() const final { return std::string(type); }

private:
  void transform_impl(GuardedOptions& state) override {
    GuardedOptions species = state["species"];
    calculate_rates(
        species["e"],                             // Electrons
        species[lithium_species_name<level + 1>], // From this ionisation state
        species[{Hisotope}],                      //    and this neutral hydrogen atom
        species[lithium_species_name<level>],     // To this state
        species[{Hisotope, '+'}]                  //    and this hydrogen ion
    );
  }
};

namespace {
// Ionisation by electron-impact
RegisterComponent<ADASLithiumIonisation<0>> register_ionisation_li0;
RegisterComponent<ADASLithiumIonisation<1>> register_ionisation_li1;
RegisterComponent<ADASLithiumIonisation<2>> register_ionisation_li2;

// Recombination
RegisterComponent<ADASLithiumRecombination<0>> register_recombination_li0;
RegisterComponent<ADASLithiumRecombination<1>> register_recombination_li1;
RegisterComponent<ADASLithiumRecombination<2>> register_recombination_li2;

// Charge exchange
RegisterComponent<ADASLithiumCX<0, 'h'>> register_cx_li0h;
RegisterComponent<ADASLithiumCX<1, 'h'>> register_cx_li1h;
RegisterComponent<ADASLithiumCX<2, 'h'>> register_cx_li2h;

RegisterComponent<ADASLithiumCX<0, 'd'>> register_cx_li0d;
RegisterComponent<ADASLithiumCX<1, 'd'>> register_cx_li1d;
RegisterComponent<ADASLithiumCX<2, 'd'>> register_cx_li2d;

RegisterComponent<ADASLithiumCX<0, 't'>> register_cx_li0t;
RegisterComponent<ADASLithiumCX<1, 't'>> register_cx_li1t;
RegisterComponent<ADASLithiumCX<2, 't'>> register_cx_li2t;
} // namespace

#endif // ADAS_LITHIUM_H
