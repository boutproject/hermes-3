#pragma once
#ifndef ADAS_NEON_H
#define ADAS_NEON_H

#include "adas_reaction.hxx"

#include <initializer_list>

/// Ionisation energies in eV
/// from https://www.webelements.com/neon/atoms.html
/// Conversion 1 kJ mol‑1 = 1.0364e-2 eV
/// These are added (removed) from the electron energy during recombination (ionisation)
constexpr std::array<BoutReal, 10> neon_ionisation_energy{
    21.56, 40.96, 63.42, 97.19, 126.24, 157.93, 207.27, 239.09, 1195.78, 1362.16};

/// The name of the species. This initializer list can be passed to a string
/// constructor, or used to index into an Options tree.
///
/// ne, ne+, ne+2, ne+3, ...
///
/// Special cases for level=0 and level=1
///
/// @tparam level  The ionisation level: 0 is neutral, 10 is fully stripped.
template <int level>
constexpr std::initializer_list<char> neon_species_name{'n', 'e', '+', '0' + level};

template <>
constexpr std::initializer_list<char> neon_species_name<1>{'n', 'e', '+'};

template <>
constexpr std::initializer_list<char> neon_species_name<0>{'n', 'e'};

/// ADAS effective ionisation (ADF11)
///
/// @tparam level  The ionisation level of the ion on the left of the reaction
template <int level>
struct ADASNeonIonisation : public OpenADAS {
  ADASNeonIonisation(std::string, Options& alloptions, Solver*)
      : OpenADAS(alloptions["units"], "scd96_ne.json", "plt96_ne.json", level,
                 -neon_ionisation_energy[level]) {}

  void transform(Options& state) override {
    calculate_rates(
        state["species"]["e"],                         // Electrons
        state["species"][neon_species_name<level>],    // From this ionisation state
        state["species"][neon_species_name<level + 1>] // To this state
    );
  }
};

/////////////////////////////////////////////////

/// ADAS effective recombination coefficients (ADF11)
///
/// @tparam level  The ionisation level of the ion on the right of the reaction
template <int level>
struct ADASNeonRecombination : public OpenADAS {
  /// @param alloptions  The top-level options. Only uses the ["units"] subsection.
  ADASNeonRecombination(std::string, Options& alloptions, Solver*)
      : OpenADAS(alloptions["units"], "acd96_ne.json", "prb96_ne.json", level,
                 neon_ionisation_energy[level]) {}

  void transform(Options& state) override {
    calculate_rates(
        state["species"]["e"],                          // Electrons
        state["species"][neon_species_name<level + 1>], // From this ionisation state
        state["species"][neon_species_name<level>]      // To this state
    );
  }
};

/// @tparam level     The ionisation level of the ion on the right of the reaction
/// @tparam Hisotope  The hydrogen isotope ('h', 'd' or 't')
template <int level, char Hisotope>
struct ADASNeonCX : public OpenADASChargeExchange {
  /// @param alloptions  The top-level options. Only uses the ["units"] subsection.
  ADASNeonCX(std::string, Options& alloptions, Solver*)
      : OpenADASChargeExchange(alloptions["units"], "ccd89_ne.json", level) {}

  void transform(Options& state) override {
    Options& species = state["species"];
    calculate_rates(
        species["e"],                          // Electrons
        species[neon_species_name<level + 1>], // From this ionisation state
        species[{Hisotope}],                   //    and this neutral hydrogen atom
        species[neon_species_name<level>],     // To this state
        species[{Hisotope, '+'}]               //    and this hydrogen ion
    );
  }
};

namespace {
// Ionisation by electron-impact
RegisterComponent<ADASNeonIonisation<0>> register_ionisation_0("ne + e -> ne+ + 2e");
RegisterComponent<ADASNeonIonisation<1>> register_ionisation_1("ne+ + e -> ne+2 + 2e");
RegisterComponent<ADASNeonIonisation<2>> register_ionisation_2("ne+2 + e -> ne+3 + 2e");
RegisterComponent<ADASNeonIonisation<3>> register_ionisation_3("ne+3 + e -> ne+4 + 2e");
RegisterComponent<ADASNeonIonisation<4>> register_ionisation_4("ne+4 + e -> ne+5 + 2e");
RegisterComponent<ADASNeonIonisation<5>> register_ionisation_5("ne+5 + e -> ne+6 + 2e");
RegisterComponent<ADASNeonIonisation<6>> register_ionisation_6("ne+6 + e -> ne+7 + 2e");
RegisterComponent<ADASNeonIonisation<7>> register_ionisation_7("ne+7 + e -> ne+8 + 2e");
RegisterComponent<ADASNeonIonisation<8>> register_ionisation_8("ne+8 + e -> ne+9 + 2e");
RegisterComponent<ADASNeonIonisation<9>> register_ionisation_9("ne+9 + e -> ne+10 + 2e");

// Recombination
RegisterComponent<ADASNeonRecombination<0>> register_recombination_0("ne+ + e -> ne");
RegisterComponent<ADASNeonRecombination<1>> register_recombination_1("ne+2 + e -> ne+");
RegisterComponent<ADASNeonRecombination<2>> register_recombination_2("ne+3 + e -> ne+2");
RegisterComponent<ADASNeonRecombination<3>> register_recombination_3("ne+4 + e -> ne+3");
RegisterComponent<ADASNeonRecombination<4>> register_recombination_4("ne+5 + e -> ne+4");
RegisterComponent<ADASNeonRecombination<5>> register_recombination_5("ne+6 + e -> ne+5");
RegisterComponent<ADASNeonRecombination<6>> register_recombination_6("ne+7 + e -> ne+6");
RegisterComponent<ADASNeonRecombination<7>> register_recombination_7("ne+8 + e -> ne+7");
RegisterComponent<ADASNeonRecombination<8>> register_recombination_8("ne+9 + e -> ne+8");
RegisterComponent<ADASNeonRecombination<9>> register_recombination_9("ne+10 + e -> ne+9");

// Charge exchange
RegisterComponent<ADASNeonCX<0, 'h'>> register_cx_0h("ne+ + h -> ne + h+");
RegisterComponent<ADASNeonCX<1, 'h'>> register_cx_1h("ne+2 + h -> ne+ + h+");
RegisterComponent<ADASNeonCX<2, 'h'>> register_cx_2h("ne+3 + h -> ne+2 + h+");
RegisterComponent<ADASNeonCX<3, 'h'>> register_cx_3h("ne+4 + h -> ne+3 + h+");
RegisterComponent<ADASNeonCX<4, 'h'>> register_cx_4h("ne+5 + h -> ne+4 + h+");
RegisterComponent<ADASNeonCX<5, 'h'>> register_cx_5h("ne+6 + h -> ne+5 + h+");
RegisterComponent<ADASNeonCX<6, 'h'>> register_cx_6h("ne+7 + h -> ne+6 + h+");
RegisterComponent<ADASNeonCX<7, 'h'>> register_cx_7h("ne+8 + h -> ne+7 + h+");
RegisterComponent<ADASNeonCX<8, 'h'>> register_cx_8h("ne+9 + h -> ne+8 + h+");
RegisterComponent<ADASNeonCX<9, 'h'>> register_cx_9h("ne+10 + h -> ne+9 + h+");

RegisterComponent<ADASNeonCX<0, 'd'>> register_cx_0d("ne+ + d -> ne + d+");
RegisterComponent<ADASNeonCX<1, 'd'>> register_cx_1d("ne+2 + d -> ne+ + d+");
RegisterComponent<ADASNeonCX<2, 'd'>> register_cx_2d("ne+3 + d -> ne+2 + d+");
RegisterComponent<ADASNeonCX<3, 'd'>> register_cx_3d("ne+4 + d -> ne+3 + d+");
RegisterComponent<ADASNeonCX<4, 'd'>> register_cx_4d("ne+5 + d -> ne+4 + d+");
RegisterComponent<ADASNeonCX<5, 'd'>> register_cx_5d("ne+6 + d -> ne+5 + d+");
RegisterComponent<ADASNeonCX<6, 'd'>> register_cx_6d("ne+7 + d -> ne+6 + d+");
RegisterComponent<ADASNeonCX<7, 'd'>> register_cx_7d("ne+8 + d -> ne+7 + d+");
RegisterComponent<ADASNeonCX<8, 'd'>> register_cx_8d("ne+9 + d -> ne+8 + d+");
RegisterComponent<ADASNeonCX<9, 'd'>> register_cx_9d("ne+10 + d -> ne+9 + d+");

RegisterComponent<ADASNeonCX<0, 't'>> register_cx_0t("ne+ + t -> ne + t+");
RegisterComponent<ADASNeonCX<1, 't'>> register_cx_1t("ne+2 + t -> ne+ + t+");
RegisterComponent<ADASNeonCX<2, 't'>> register_cx_2t("ne+3 + t -> ne+2 + t+");
RegisterComponent<ADASNeonCX<3, 't'>> register_cx_3t("ne+4 + t -> ne+3 + t+");
RegisterComponent<ADASNeonCX<4, 't'>> register_cx_4t("ne+5 + t -> ne+4 + t+");
RegisterComponent<ADASNeonCX<5, 't'>> register_cx_5t("ne+6 + t -> ne+5 + t+");
RegisterComponent<ADASNeonCX<6, 't'>> register_cx_6t("ne+7 + t -> ne+6 + t+");
RegisterComponent<ADASNeonCX<7, 't'>> register_cx_7t("ne+8 + t -> ne+7 + t+");
RegisterComponent<ADASNeonCX<8, 't'>> register_cx_8t("ne+9 + t -> ne+8 + t+");
RegisterComponent<ADASNeonCX<9, 't'>> register_cx_9t("ne+10 + t -> ne+9 + t+");
} // namespace

#endif // ADAS_NEON_H
