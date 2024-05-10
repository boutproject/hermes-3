#pragma once
#ifndef ADAS_LITHIUM_H
#define ADAS_LITHIUM_H

#include "adas_reaction.hxx"

#include <array>
#include <initializer_list>

/// Ionisation energies in eV
/// from https://www.webelements.com/lithium/atoms.html
/// Conversion 1 kJ mol‑1 = 1.0364e-2 eV
/// These are added (removed) from the electron energy during recombination (ionisation)
constexpr std::array<BoutReal, 3> lithium_ionisation_energy{
    5.39, 75.64, 122.45};

/// The name of the species. This initializer list can be passed to a string
/// constructor, or used to index into an Options tree.
///
/// li, li+, li+2, ...
///
/// Special cases for level=0, 1 and 3
///
/// @tparam level  The ionisation level: 0 is neutral, 3 is fully stripped.
template <int level>
constexpr std::initializer_list<char> lithium_species_name{'l', 'i', '+', '0' + level};

template <>
constexpr std::initializer_list<char> lithium_species_name<3>{'l', 'i', '+', '3'};

template <>
constexpr std::initializer_list<char> lithium_species_name<1>{'l', 'i', '+'};

template <>
constexpr std::initializer_list<char> lithium_species_name<0>{'l', 'i'};

/// ADAS effective ionisation (ADF11)
///
/// @tparam level  The ionisation level of the ion on the left of the reaction
template <int level>
struct ADASLithiumIonisation : public OpenADAS {
  ADASLithiumIonisation(std::string, Options& alloptions, Solver*)
      : OpenADAS(alloptions["units"], "scd96_li.json", "plt96_li.json", level,
                 -lithium_ionisation_energy[level]) {}

  void transform(Options& state) override {
    calculate_rates(
        state["species"]["e"],                         // Electrons
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
  ADASLithiumRecombination(std::string, Options& alloptions, Solver*)
      : OpenADAS(alloptions["units"], "acd96_li.json", "prb96_li.json", level,
                 lithium_ionisation_energy[level]) {}

  void transform(Options& state) override {
    calculate_rates(
        state["species"]["e"],                          // Electrons
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
  ADASLithiumCX(std::string, Options& alloptions, Solver*)
      : OpenADASChargeExchange(alloptions["units"], "ccd89_li.json", level) {}

  void transform(Options& state) override {
    Options& species = state["species"];
    calculate_rates(
        species["e"],                          // Electrons
        species[lithium_species_name<level + 1>], // From this ionisation state
        species[{Hisotope}],                   //    and this neutral hydrogen atom
        species[lithium_species_name<level>],     // To this state
        species[{Hisotope, '+'}]               //    and this hydrogen ion
    );
  }
};

namespace {
// Ionisation by electron-impact
RegisterComponent<ADASLithiumIonisation<0>> register_ionisation_li0("li + e -> li+ + 2e");
RegisterComponent<ADASLithiumIonisation<1>> register_ionisation_li1("li+ + e -> li+2 + 2e");
RegisterComponent<ADASLithiumIonisation<2>> register_ionisation_li2("li+2 + e -> li+3 + 2e");

// Recombination
RegisterComponent<ADASLithiumRecombination<0>> register_recombination_li0("li+ + e -> li");
RegisterComponent<ADASLithiumRecombination<1>> register_recombination_li1("li+2 + e -> li+");
RegisterComponent<ADASLithiumRecombination<2>> register_recombination_li2("li+3 + e -> li+2");

// Charge exchange
RegisterComponent<ADASLithiumCX<0, 'h'>> register_cx_li0h("li+ + h -> li + h+");
RegisterComponent<ADASLithiumCX<1, 'h'>> register_cx_li1h("li+2 + h -> li+ + h+");
RegisterComponent<ADASLithiumCX<2, 'h'>> register_cx_li2h("li+3 + h -> li+2 + h+");

RegisterComponent<ADASLithiumCX<0, 'd'>> register_cx_li0d("li+ + d -> li + d+");
RegisterComponent<ADASLithiumCX<1, 'd'>> register_cx_li1d("li+2 + d -> li+ + d+");
RegisterComponent<ADASLithiumCX<2, 'd'>> register_cx_li2d("li+3 + d -> li+2 + d+");

RegisterComponent<ADASLithiumCX<0, 't'>> register_cx_li0t("li+ + t -> li + t+");
RegisterComponent<ADASLithiumCX<1, 't'>> register_cx_li1t("li+2 + t -> li+ + t+");
RegisterComponent<ADASLithiumCX<2, 't'>> register_cx_li2t("li+3 + t -> li+2 + t+");
} // namespace

#endif // ADAS_LITHIUM_H
