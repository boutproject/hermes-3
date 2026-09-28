#pragma once
#ifndef ADAS_NEON_H
#define ADAS_NEON_H

#include "adas_reaction.hxx"

#include <array>
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
/// Special cases for level=0, 1 and 10
///
/// @tparam level  The ionisation level: 0 is neutral, 10 is fully stripped.
template <int level>
constexpr char neon_species_name[6]{'n', 'e', '+', '0' + level, '\0'};

template <>
constexpr char neon_species_name<10>[6]{'n', 'e', '+', '1', '0', '\0'};

template <>
constexpr char neon_species_name<1>[6]{'n', 'e', '+', '\0'};

template <>
constexpr char neon_species_name<0>[6]{'n', 'e', '\0'};

template <int level>
struct neon_species {
  static constexpr auto value = neon_species_name<level>;
};

/// ADAS effective ionisation (ADF11)
///
/// @tparam level  The ionisation level of the ion on the left of the reaction
template <std::size_t level>
struct ADASNeonIonisation : public OpenADAS {
  ADASNeonIonisation(std::string name, Options& alloptions, Solver*)
      : OpenADAS(std::move(name), alloptions["units"], "scd96_ne.json", "plt96_ne.json",
                 neon_species_name<level>, neon_species_name<level + 1>, level,
                 -neon_ionisation_energy[level]) {}
  static constexpr auto type = izn_component_name_v<neon_species, level>;

  std::string typeName() const final { return std::string(type); }

private:
  void transform_impl(GuardedOptions& state) override {
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
template <std::size_t level>
struct ADASNeonRecombination : public OpenADAS {
  /// @param alloptions  The top-level options. Only uses the ["units"] subsection.
  ADASNeonRecombination(std::string name, Options& alloptions, Solver*)
      : OpenADAS(std::move(name), alloptions["units"], "acd96_ne.json", "prb96_ne.json",
                 neon_species_name<level + 1>, neon_species_name<level>, level,
                 neon_ionisation_energy[level]) {}
  static constexpr auto type = rec_component_name_v<neon_species, level>;

  std::string typeName() const final { return std::string(type); }

private:
  void transform_impl(GuardedOptions& state) override {
    calculate_rates(
        state["species"]["e"],                          // Electrons
        state["species"][neon_species_name<level + 1>], // From this ionisation state
        state["species"][neon_species_name<level>]      // To this state
    );
  }
};

/// @tparam level     The ionisation level of the ion on the right of the reaction
/// @tparam Hisotope  The hydrogen isotope ('h', 'd' or 't')
template <std::size_t level, char Hisotope>
struct ADASNeonCX : public OpenADASChargeExchange {
  /// @param alloptions  The top-level options. Only uses the ["units"] subsection.
  ADASNeonCX(std::string name, Options& alloptions, Solver*)
      : OpenADASChargeExchange(std::move(name), alloptions["units"], "ccd89_ne.json",
                               neon_species_name<level + 1>, {Hisotope},
                               neon_species_name<level>, {Hisotope, '+'}, level) {}
  static constexpr auto type = cx_component_name_v<neon_species, level, Hisotope>;

  std::string typeName() const final { return std::string(type); }

private:
  void transform_impl(GuardedOptions& state) override {
    GuardedOptions species = state["species"];
    calculate_rates(species["e"],                          // Electrons
                    species[neon_species_name<level + 1>], // From this ionisation state
                    species[{Hisotope}], //    and this neutral hydrogen atom
                    species[neon_species_name<level>], // To this state
                    species[{Hisotope, '+'}]           //    and this hydrogen ion
    );
  }
};

namespace {
// Ionisation by electron-impact
RegisterComponent<ADASNeonIonisation<0>> register_ionisation_ne0;
RegisterComponent<ADASNeonIonisation<1>> register_ionisation_ne1;
RegisterComponent<ADASNeonIonisation<2>> register_ionisation_ne2;
RegisterComponent<ADASNeonIonisation<3>> register_ionisation_ne3;
RegisterComponent<ADASNeonIonisation<4>> register_ionisation_ne4;
RegisterComponent<ADASNeonIonisation<5>> register_ionisation_ne5;
RegisterComponent<ADASNeonIonisation<6>> register_ionisation_ne6;
RegisterComponent<ADASNeonIonisation<7>> register_ionisation_ne7;
RegisterComponent<ADASNeonIonisation<8>> register_ionisation_ne8;
RegisterComponent<ADASNeonIonisation<9>> register_ionisation_ne9;

// Recombination
RegisterComponent<ADASNeonRecombination<0>> register_recombination_ne0;
RegisterComponent<ADASNeonRecombination<1>> register_recombination_ne1;
RegisterComponent<ADASNeonRecombination<2>> register_recombination_ne2;
RegisterComponent<ADASNeonRecombination<3>> register_recombination_ne3;
RegisterComponent<ADASNeonRecombination<4>> register_recombination_ne4;
RegisterComponent<ADASNeonRecombination<5>> register_recombination_ne5;
RegisterComponent<ADASNeonRecombination<6>> register_recombination_ne6;
RegisterComponent<ADASNeonRecombination<7>> register_recombination_ne7;
RegisterComponent<ADASNeonRecombination<8>> register_recombination_ne8;
RegisterComponent<ADASNeonRecombination<9>> register_recombination_ne9;

// Charge exchange
RegisterComponent<ADASNeonCX<0, 'h'>> register_cx_ne0h;
RegisterComponent<ADASNeonCX<1, 'h'>> register_cx_ne1h;
RegisterComponent<ADASNeonCX<2, 'h'>> register_cx_ne2h;
RegisterComponent<ADASNeonCX<3, 'h'>> register_cx_ne3h;
RegisterComponent<ADASNeonCX<4, 'h'>> register_cx_ne4h;
RegisterComponent<ADASNeonCX<5, 'h'>> register_cx_ne5h;
RegisterComponent<ADASNeonCX<6, 'h'>> register_cx_ne6h;
RegisterComponent<ADASNeonCX<7, 'h'>> register_cx_ne7h;
RegisterComponent<ADASNeonCX<8, 'h'>> register_cx_ne8h;
RegisterComponent<ADASNeonCX<9, 'h'>> register_cx_ne9h;

RegisterComponent<ADASNeonCX<0, 'd'>> register_cx_ne0d;
RegisterComponent<ADASNeonCX<1, 'd'>> register_cx_ne1d;
RegisterComponent<ADASNeonCX<2, 'd'>> register_cx_ne2d;
RegisterComponent<ADASNeonCX<3, 'd'>> register_cx_ne3d;
RegisterComponent<ADASNeonCX<4, 'd'>> register_cx_ne4d;
RegisterComponent<ADASNeonCX<5, 'd'>> register_cx_ne5d;
RegisterComponent<ADASNeonCX<6, 'd'>> register_cx_ne6d;
RegisterComponent<ADASNeonCX<7, 'd'>> register_cx_ne7d;
RegisterComponent<ADASNeonCX<8, 'd'>> register_cx_ne8d;
RegisterComponent<ADASNeonCX<9, 'd'>> register_cx_ne9d;

RegisterComponent<ADASNeonCX<0, 't'>> register_cx_ne0t;
RegisterComponent<ADASNeonCX<1, 't'>> register_cx_ne1t;
RegisterComponent<ADASNeonCX<2, 't'>> register_cx_ne2t;
RegisterComponent<ADASNeonCX<3, 't'>> register_cx_ne3t;
RegisterComponent<ADASNeonCX<4, 't'>> register_cx_ne4t;
RegisterComponent<ADASNeonCX<5, 't'>> register_cx_ne5t;
RegisterComponent<ADASNeonCX<6, 't'>> register_cx_ne6t;
RegisterComponent<ADASNeonCX<7, 't'>> register_cx_ne7t;
RegisterComponent<ADASNeonCX<8, 't'>> register_cx_ne8t;
RegisterComponent<ADASNeonCX<9, 't'>> register_cx_ne9t;
} // namespace

#endif // ADAS_NEON_H
