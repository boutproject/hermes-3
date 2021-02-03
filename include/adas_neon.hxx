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
template <int level>
constexpr std::initializer_list<char> neon_species_name {'n', 'e', '+', '0' + level};

template <>
constexpr std::initializer_list<char> neon_species_name<1> {'n', 'e', '+'};

template <>
constexpr std::initializer_list<char> neon_species_name<0> {'n', 'e'};

/// ADAS effective ionisation (ADF11)
///
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
template <int level>
struct ADASNeonRecombination : public OpenADAS {
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
} // namespace

#endif // ADAS_NEON_H
