#pragma once
#ifndef ADAS_ARGON_H
#define ADAS_ARGON_H

#include "adas_reaction.hxx"

#include <array>
#include <initializer_list>

/// Argon rates and radiation from the OpenADAS ADF11 files scd89_ar, plt89_ar,
/// acd89_ar and prb89_ar. These files are not included in Hermes-3; create the
/// JSON files with
///
///     python3 json_database/download_openadas.py ar

/// Ionisation energies in eV
/// from the NIST Atomic Spectra Database (ASD), Ar I to Ar XVIII
/// https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
/// These are added (removed) from the electron energy during recombination (ionisation)
constexpr std::array<BoutReal, 18> argon_ionisation_energy{
    15.7596119, 27.62967, 40.735, 59.58, 74.84,  91.290, 124.41,  143.4567,   422.60,
    479.76,     540.4,    619.0,  685.5, 755.13, 855.5,  918.375, 4120.66559, 4426.22407};

/// The name of the species. This initializer list can be passed to a string
/// constructor, or used to index into an Options tree.
///
/// ar, ar+, ar+2, ar+3, ...
///
/// Special cases for level=0, 1 and 10 to 18
///
/// @tparam level  The ionisation level: 0 is neutral, 18 is fully stripped.
template <int level>
constexpr char argon_species_name[6]{'a', 'r', '+', '0' + level, '\0'};

template <>
constexpr char argon_species_name<18>[6]{'a', 'r', '+', '1', '8', '\0'};

template <>
constexpr char argon_species_name<17>[6]{'a', 'r', '+', '1', '7', '\0'};

template <>
constexpr char argon_species_name<16>[6]{'a', 'r', '+', '1', '6', '\0'};

template <>
constexpr char argon_species_name<15>[6]{'a', 'r', '+', '1', '5', '\0'};

template <>
constexpr char argon_species_name<14>[6]{'a', 'r', '+', '1', '4', '\0'};

template <>
constexpr char argon_species_name<13>[6]{'a', 'r', '+', '1', '3', '\0'};

template <>
constexpr char argon_species_name<12>[6]{'a', 'r', '+', '1', '2', '\0'};

template <>
constexpr char argon_species_name<11>[6]{'a', 'r', '+', '1', '1', '\0'};

template <>
constexpr char argon_species_name<10>[6]{'a', 'r', '+', '1', '0', '\0'};

template <>
constexpr char argon_species_name<1>[6]{'a', 'r', '+', '\0'};

template <>
constexpr char argon_species_name<0>[6]{'a', 'r', '\0'};

template <int level>
struct argon_species {
  static constexpr auto value = argon_species_name<level>;
};

/// ADAS effective ionisation (ADF11)
///
/// @tparam level  The ionisation level of the ion on the left of the reaction
template <std::size_t level>
struct ADASArgonIonisation : public OpenADAS {
  ADASArgonIonisation(std::string name, Options& alloptions, Solver*)
      : OpenADAS(std::move(name), alloptions, "scd89_ar.json", "plt89_ar.json",
                 argon_species_name<level>, argon_species_name<level + 1>, level,
                 -argon_ionisation_energy[level]) {}
  static constexpr auto type = izn_component_name_v<argon_species, level>;

  std::string typeName() const final { return std::string(type); }

private:
  void transform_impl(GuardedOptions& state) override {
    calculate_rates(
        state["species"]["e"],                          // Electrons
        state["species"][argon_species_name<level>],    // From this ionisation state
        state["species"][argon_species_name<level + 1>] // To this state
    );
  }
};

/////////////////////////////////////////////////

/// ADAS effective recombination coefficients (ADF11)
///
/// @tparam level  The ionisation level of the ion on the right of the reaction
template <std::size_t level>
struct ADASArgonRecombination : public OpenADAS {
  /// @param alloptions  The top-level options.
  ADASArgonRecombination(std::string name, Options& alloptions, Solver*)
      : OpenADAS(std::move(name), alloptions, "acd89_ar.json", "prb89_ar.json",
                 argon_species_name<level + 1>, argon_species_name<level>, level,
                 argon_ionisation_energy[level]) {}
  static constexpr auto type = rec_component_name_v<argon_species, level>;

  std::string typeName() const final { return std::string(type); }

private:
  void transform_impl(GuardedOptions& state) override {
    calculate_rates(
        state["species"]["e"],                           // Electrons
        state["species"][argon_species_name<level + 1>], // From this ionisation state
        state["species"][argon_species_name<level>]      // To this state
    );
  }
};

namespace {
// Ionisation by electron-impact
RegisterComponent<ADASArgonIonisation<0>> register_ionisation_ar0;
RegisterComponent<ADASArgonIonisation<1>> register_ionisation_ar1;
RegisterComponent<ADASArgonIonisation<2>> register_ionisation_ar2;
RegisterComponent<ADASArgonIonisation<3>> register_ionisation_ar3;
RegisterComponent<ADASArgonIonisation<4>> register_ionisation_ar4;
RegisterComponent<ADASArgonIonisation<5>> register_ionisation_ar5;
RegisterComponent<ADASArgonIonisation<6>> register_ionisation_ar6;
RegisterComponent<ADASArgonIonisation<7>> register_ionisation_ar7;
RegisterComponent<ADASArgonIonisation<8>> register_ionisation_ar8;
RegisterComponent<ADASArgonIonisation<9>> register_ionisation_ar9;
RegisterComponent<ADASArgonIonisation<10>> register_ionisation_ar10;
RegisterComponent<ADASArgonIonisation<11>> register_ionisation_ar11;
RegisterComponent<ADASArgonIonisation<12>> register_ionisation_ar12;
RegisterComponent<ADASArgonIonisation<13>> register_ionisation_ar13;
RegisterComponent<ADASArgonIonisation<14>> register_ionisation_ar14;
RegisterComponent<ADASArgonIonisation<15>> register_ionisation_ar15;
RegisterComponent<ADASArgonIonisation<16>> register_ionisation_ar16;
RegisterComponent<ADASArgonIonisation<17>> register_ionisation_ar17;

// Recombination
RegisterComponent<ADASArgonRecombination<0>> register_recombination_ar0;
RegisterComponent<ADASArgonRecombination<1>> register_recombination_ar1;
RegisterComponent<ADASArgonRecombination<2>> register_recombination_ar2;
RegisterComponent<ADASArgonRecombination<3>> register_recombination_ar3;
RegisterComponent<ADASArgonRecombination<4>> register_recombination_ar4;
RegisterComponent<ADASArgonRecombination<5>> register_recombination_ar5;
RegisterComponent<ADASArgonRecombination<6>> register_recombination_ar6;
RegisterComponent<ADASArgonRecombination<7>> register_recombination_ar7;
RegisterComponent<ADASArgonRecombination<8>> register_recombination_ar8;
RegisterComponent<ADASArgonRecombination<9>> register_recombination_ar9;
RegisterComponent<ADASArgonRecombination<10>> register_recombination_ar10;
RegisterComponent<ADASArgonRecombination<11>> register_recombination_ar11;
RegisterComponent<ADASArgonRecombination<12>> register_recombination_ar12;
RegisterComponent<ADASArgonRecombination<13>> register_recombination_ar13;
RegisterComponent<ADASArgonRecombination<14>> register_recombination_ar14;
RegisterComponent<ADASArgonRecombination<15>> register_recombination_ar15;
RegisterComponent<ADASArgonRecombination<16>> register_recombination_ar16;
RegisterComponent<ADASArgonRecombination<17>> register_recombination_ar17;
} // namespace

#endif // ADAS_ARGON_H
