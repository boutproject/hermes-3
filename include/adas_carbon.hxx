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
constexpr std::array<BoutReal, 6> carbon_ionisation_energy{11.26, 24.38,  47.89,
                                                           64.49, 392.09, 489.99};

/// The name of the species. This initializer list can be passed to a string
/// constructor, or used to index into an Options tree.
///
/// c, c+, c+2, c+3, ...
///
/// Special cases for level=0, 1
///
/// @tparam level  The ionisation level: 0 is neutral, 6 is fully stripped
template <int level>
constexpr char carbon_species_name[4] = {'c', '+', '0' + level, '\0'};

template <>
constexpr char carbon_species_name<1>[4]{'c', '+', '\0'};

template <>
constexpr char carbon_species_name<0>[4]{'c', '\0'};

template <int level>
struct carbon_species {
  static constexpr auto value = carbon_species_name<level>;
};

/// ADAS effective ionisation (ADF11)
///
/// @tparam level  The ionisation level of the ion on the left of the reaction
template <int level>
struct ADASCarbonIonisation : public OpenADAS {
  ADASCarbonIonisation(std::string name, Options& alloptions, Solver*)
      : OpenADAS(std::move(name), alloptions["units"], "scd96_c.json", "plt96_c.json",
                 carbon_species_name<level>, carbon_species_name<level + 1>, level,
                 -carbon_ionisation_energy[level]) {}

  static constexpr auto type = izn_component_name_v<carbon_species, level>;

  std::string typeName() const final { return std::string(type); }

private:
  void transform_impl(GuardedOptions& state) override {
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
  ADASCarbonRecombination(std::string name, Options& alloptions, Solver*)
      : OpenADAS(std::move(name), alloptions["units"], "acd96_c.json", "prb96_c.json",
                 carbon_species_name<level + 1>, carbon_species_name<level>, level,
                 carbon_ionisation_energy[level]) {}

  static constexpr auto type = rec_component_name_v<carbon_species, level>;

  std::string typeName() const final { return std::string(type); }

private:
  void transform_impl(GuardedOptions& state) override {
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
  ADASCarbonCX(std::string name, Options& alloptions, Solver*)
      : OpenADASChargeExchange(std::move(name), alloptions["units"], "ccd96_c.json",
                               carbon_species_name<level + 1>, {Hisotope},
                               carbon_species_name<level>, {Hisotope, '+'}, level) {}

  static constexpr auto type = cx_component_name_v<carbon_species, level, Hisotope>;

  std::string typeName() const final { return std::string(type); }

private:
  void transform_impl(GuardedOptions& state) override {
    GuardedOptions species = state["species"];
    calculate_rates(species["e"],                            // Electrons
                    species[carbon_species_name<level + 1>], // From this ionisation state
                    species[{Hisotope}], //    and this neutral hydrogen atom
                    species[carbon_species_name<level>], // To this state
                    species[{Hisotope, '+'}]             //    and this hydrogen ion
    );
  }
};

namespace {
// Ionisation by electron-impact
RegisterComponent<ADASCarbonIonisation<0>> register_ionisation_c0;
RegisterComponent<ADASCarbonIonisation<1>> register_ionisation_c1;
RegisterComponent<ADASCarbonIonisation<2>> register_ionisation_c2;
RegisterComponent<ADASCarbonIonisation<3>> register_ionisation_c3;
RegisterComponent<ADASCarbonIonisation<4>> register_ionisation_c4;
RegisterComponent<ADASCarbonIonisation<5>> register_ionisation_c5;

// Recombination
RegisterComponent<ADASCarbonRecombination<0>> register_recombination_c0;
RegisterComponent<ADASCarbonRecombination<1>> register_recombination_c1;
RegisterComponent<ADASCarbonRecombination<2>> register_recombination_c2;
RegisterComponent<ADASCarbonRecombination<3>> register_recombination_c3;
RegisterComponent<ADASCarbonRecombination<4>> register_recombination_c4;
RegisterComponent<ADASCarbonRecombination<5>> register_recombination_c5;

// Charge exchange
RegisterComponent<ADASCarbonCX<0, 'h'>> register_cx_c0h;
RegisterComponent<ADASCarbonCX<1, 'h'>> register_cx_c1h;
RegisterComponent<ADASCarbonCX<2, 'h'>> register_cx_c2h;
RegisterComponent<ADASCarbonCX<3, 'h'>> register_cx_c3h;
RegisterComponent<ADASCarbonCX<4, 'h'>> register_cx_c4h;
RegisterComponent<ADASCarbonCX<5, 'h'>> register_cx_c5h;

RegisterComponent<ADASCarbonCX<0, 'd'>> register_cx_c0d;
RegisterComponent<ADASCarbonCX<1, 'd'>> register_cx_c1d;
RegisterComponent<ADASCarbonCX<2, 'd'>> register_cx_c2d;
RegisterComponent<ADASCarbonCX<3, 'd'>> register_cx_c3d;
RegisterComponent<ADASCarbonCX<4, 'd'>> register_cx_c4d;
RegisterComponent<ADASCarbonCX<5, 'd'>> register_cx_c5d;

RegisterComponent<ADASCarbonCX<0, 't'>> register_cx_c0t;
RegisterComponent<ADASCarbonCX<1, 't'>> register_cx_c1t;
RegisterComponent<ADASCarbonCX<2, 't'>> register_cx_c2t;
RegisterComponent<ADASCarbonCX<3, 't'>> register_cx_c3t;
RegisterComponent<ADASCarbonCX<4, 't'>> register_cx_c4t;
RegisterComponent<ADASCarbonCX<5, 't'>> register_cx_c5t;

} // namespace

#endif // ADAS_CARBON_H
