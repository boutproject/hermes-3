#pragma once
#ifndef ADAS_CARBON_H
#define ADAS_CARBON_H

#include "adas_reaction.hxx"

#include <array>
#include <initializer_list>

/// Ionisation energies in eV
/// from https://www.webelements.com/carbon/atoms.html
/// Conversion 1 kJ mol‑1 = 1.0364e-2 eV
/// These are added (removed) from the electron energy during recombination (ionisation)
constexpr std::array<BoutReal, 6> carbon_ionisation_energy{
    11.26, 24.38, 47.89, 64.49, 392.09, 489.99};

/// The name of the species. This initializer list can be passed to a string
/// constructor, or used to index into an Options tree.
///
/// c, c+, c+2, c+3, ...
///
/// Special cases for level=0, 1
///
/// @tparam level  The ionisation level: 0 is neutral, 6 is fully stripped
template <int level>
constexpr std::initializer_list<char> carbon_species_name{'c', '+', '0' + level};

template <>
constexpr std::initializer_list<char> carbon_species_name<1>{'c', '+'};

template <>
constexpr std::initializer_list<char> carbon_species_name<0>{'c'};

/// ADAS effective ionisation (ADF11)
///
/// @tparam level  The ionisation level of the ion on the left of the reaction
template <int level>
struct ADASCarbonIonisation : public OpenADAS {
  ADASCarbonIonisation(std::string, Options& alloptions, Solver*)
      : OpenADAS(alloptions["units"], "scd96_c.json", "plt96_c.json", level,
                 -carbon_ionisation_energy[level]) {}

  void transform(Options& state) override {
    calculate_rates(
        state["species"]["e"],                           // Electrons
        state["species"][carbon_species_name<level>],    // From this ionisation state
        state["species"][carbon_species_name<level + 1>] // To this state
    );
  }
};

/////////////////////////////////////////////////

/// ADAS effective recombination coefficients (ADF11)
///
/// @tparam level  The ionisation level of the ion on the right of the reaction
template <int level>
struct ADASCarbonRecombination : public OpenADAS {
  /// @param alloptions  The top-level options. Only uses the ["units"] subsection.
  ADASCarbonRecombination(std::string, Options& alloptions, Solver*)
      : OpenADAS(alloptions["units"], "acd96_c.json", "prb96_c.json", level,
                 carbon_ionisation_energy[level]) {}

  void transform(Options& state) override {
    calculate_rates(
        state["species"]["e"],                            // Electrons
        state["species"][carbon_species_name<level + 1>], // From this ionisation state
        state["species"][carbon_species_name<level>]      // To this state
    );
  }
};

/// @tparam level     The ionisation level of the ion on the right of the reaction
/// @tparam Hisotope  The hydrogen isotope ('h', 'd' or 't')
template <int level, char Hisotope>
struct ADASCarbonCX : public OpenADASChargeExchange {
  /// @param alloptions  The top-level options. Only uses the ["units"] subsection.
  ADASCarbonCX(std::string, Options& alloptions, Solver*)
      : OpenADASChargeExchange(alloptions["units"], "ccd89_c.json", level) {}

  void transform(Options& state) override {
    Options& species = state["species"];
    calculate_rates(
        species["e"],                            // Electrons
        species[carbon_species_name<level + 1>], // From this ionisation state
        species[{Hisotope}],                     //    and this neutral hydrogen atom
        species[carbon_species_name<level>],     // To this state
        species[{Hisotope, '+'}]                 //    and this hydrogen ion
    );
  }
};

namespace {
// Ionisation by electron-impact
RegisterComponent<ADASCarbonIonisation<0>> register_ionisation_0("c + e -> c+ + 2e");
RegisterComponent<ADASCarbonIonisation<1>> register_ionisation_1("c+ + e -> c+2 + 2e");
RegisterComponent<ADASCarbonIonisation<2>> register_ionisation_2("c+2 + e -> c+3 + 2e");
RegisterComponent<ADASCarbonIonisation<3>> register_ionisation_3("c+3 + e -> c+4 + 2e");
RegisterComponent<ADASCarbonIonisation<4>> register_ionisation_4("c+4 + e -> c+5 + 2e");
RegisterComponent<ADASCarbonIonisation<5>> register_ionisation_5("c+5 + e -> c+6 + 2e");

// Recombination
RegisterComponent<ADASCarbonRecombination<0>> register_recombination_0("c+ + e -> c");
RegisterComponent<ADASCarbonRecombination<1>> register_recombination_1("c+2 + e -> c+");
RegisterComponent<ADASCarbonRecombination<2>> register_recombination_2("c+3 + e -> c+2");
RegisterComponent<ADASCarbonRecombination<3>> register_recombination_3("c+4 + e -> c+3");
RegisterComponent<ADASCarbonRecombination<4>> register_recombination_4("c+5 + e -> c+4");
RegisterComponent<ADASCarbonRecombination<5>> register_recombination_5("c+6 + e -> c+5");

// Charge exchange
RegisterComponent<ADASCarbonCX<0, 'h'>> register_cx_0h("c+ + h -> c + h+");
RegisterComponent<ADASCarbonCX<1, 'h'>> register_cx_1h("c+2 + h -> c+ + h+");
RegisterComponent<ADASCarbonCX<2, 'h'>> register_cx_2h("c+3 + h -> c+2 + h+");
RegisterComponent<ADASCarbonCX<3, 'h'>> register_cx_3h("c+4 + h -> c+3 + h+");
RegisterComponent<ADASCarbonCX<4, 'h'>> register_cx_4h("c+5 + h -> c+4 + h+");
RegisterComponent<ADASCarbonCX<5, 'h'>> register_cx_5h("c+6 + h -> c+5 + h+");

RegisterComponent<ADASCarbonCX<0, 'd'>> register_cx_0d("c+ + d -> c + d+");
RegisterComponent<ADASCarbonCX<1, 'd'>> register_cx_1d("c+2 + d -> c+ + d+");
RegisterComponent<ADASCarbonCX<2, 'd'>> register_cx_2d("c+3 + d -> c+2 + d+");
RegisterComponent<ADASCarbonCX<3, 'd'>> register_cx_3d("c+4 + d -> c+3 + d+");
RegisterComponent<ADASCarbonCX<4, 'd'>> register_cx_4d("c+5 + d -> c+4 + d+");
RegisterComponent<ADASCarbonCX<5, 'd'>> register_cx_5d("c+6 + d -> c+5 + d+");

RegisterComponent<ADASCarbonCX<0, 't'>> register_cx_0t("c+ + t -> c + t+");
RegisterComponent<ADASCarbonCX<1, 't'>> register_cx_1t("c+2 + t -> c+ + t+");
RegisterComponent<ADASCarbonCX<2, 't'>> register_cx_2t("c+3 + t -> c+2 + t+");
RegisterComponent<ADASCarbonCX<3, 't'>> register_cx_3t("c+4 + t -> c+3 + t+");
RegisterComponent<ADASCarbonCX<4, 't'>> register_cx_4t("c+5 + t -> c+4 + t+");
RegisterComponent<ADASCarbonCX<5, 't'>> register_cx_5t("c+6 + t -> c+5 + t+");

} // namespace

#endif // ADAS_CARBON_H
