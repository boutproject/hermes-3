#pragma once
#ifndef ADAS_REACTION_H
#define ADAS_REACTION_H

#include "component.hxx"
#include "reaction.hxx"

#include <cstddef>
#include <vector>

/// Represent a 2D rate coefficient table (T,n)
/// Reads data from a file, then interpolates at required values.
struct OpenADASRateCoefficient {
  /// Read the file, extracting data for the given ionisation level
  /// @param filename   The file to read. Path relative to run working directory
  /// @param level      The first index in the log coefficient array
  ///                   (ionisation level)
  OpenADASRateCoefficient(const std::string& filename, std::size_t level);

  std::vector<std::vector<BoutReal>> log_coeff;
  std::vector<BoutReal> log_temperature;
  std::vector<BoutReal> log_density;

  BoutReal Tmin, Tmax; ///< Range of T  [eV]
  BoutReal nmin, nmax; ///< Range of density [m^-3]

  /// Inputs:
  /// @param  n  Electron density in m^-3
  /// @param  T  Electron temperature in eV
  ///
  /// @returns rate in units of m^3/s or eV m^3/s
  BoutReal evaluate(BoutReal T, BoutReal n);
};

/// Read in and perform calculations with OpenADAS data
/// https://open.adas.ac.uk/
///
/// Uses the JSON files produced by:
///   https://github.com/TBody/OpenADAS_to_JSON
struct OpenADAS : public ReactionBase {
  ///
  /// Inputs
  /// ------
  /// @param units       Options tree containing normalisation constants
  /// @param rate_file   A JSON file containing reaction rate <σv> rates (e.g. SCD, ACD)
  /// @param radiation_file   A JSON file containing radiation loss rates (e.g. PLT, PRB)
  /// @param level       The lower ionisation state in the transition
  ///               e.g. 0 for neutral -> 1st ionisation
  ///               and 1st -> neutral recombination
  /// @param electron_heating   The heating of the electrons per reaction [eV]
  ///               This is the ionisation energy, positive for recombination
  ///               and negative for ionisation
  ///
  /// Notes
  ///  - The rate and radiation file names have "json_database/" prepended
  ///
  OpenADAS(const Options& units, const std::string& rate_file,
           const std::string& radiation_file, std::string from_ion, std::string to_ion,
           std::size_t level, BoutReal electron_heating)
      : ReactionBase({readIfSet("species:{sp}:charge"), readOnly("species:{sp}:AA"),
                      readOnly("species:{from_ion}:{val}"), readOnly("species:e:{e_val}"),
                      readWrite("species:{sp}:{w_val}"),
                      readWrite("species:e:{ew_val}")}),
        rate_coef(std::string("json_database/") + rate_file, level),
        radiation_coef(std::string("json_database/") + radiation_file, level),
        electron_heating(electron_heating) {
    // Get the units
    Tnorm = get<BoutReal>(units["eV"]);
    Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
    FreqNorm = 1. / get<BoutReal>(units["seconds"]);
    substitutePermissions("val", {"density", "temperature", "velocity"});
    substitutePermissions("e_val", {"density", "temperature"});
    substitutePermissions("w_val",
                          {"density_source", "momentum_source", "energy_source"});
    // FIXME: electron density_source only written if from_ion charge != to_ion charge.
    substitutePermissions("ew_val",
                          {"density_source", "momentum_source", "energy_source"});
    substitutePermissions("sp", {from_ion, to_ion});
    substitutePermissions("from_ion", {from_ion});
  }

  /// Perform the calculation of rates, and transfer of particles/momentum/energy
  ///
  /// @param electron  The electron species e.g. state["species"]["e"]
  /// @param from_ion  The ion on the left of the reaction
  /// @param to_ion    The ion on the right of the reaction
  void calculate_rates(GuardedOptions && electron, GuardedOptions && from_ion, GuardedOptions && to_ion);
private:
  OpenADASRateCoefficient rate_coef;      ///< Reaction rate coefficient
  OpenADASRateCoefficient radiation_coef; ///< Energy loss (radiation) coefficient

  BoutReal electron_heating; ///< Heating per reaction [eV]

  BoutReal Tnorm, Nnorm, FreqNorm; ///< Normalisations
};

struct OpenADASChargeExchange : public ReactionBase {
  OpenADASChargeExchange(const Options& units, const std::string& rate_file,
                         std::string from_A, std::string from_B, std::string to_A,
                         std::string to_B, std::size_t level)
      : ReactionBase({readIfSet("species:{sp}:charge"), readOnly("species:{sp}:AA"),
                      readOnly("species:{from_ion}:{val}"), readOnly("species:e:{e_val}"),
                      readWrite("species:{sp}:{w_val}")}),
        rate_coef(std::string("json_database/") + rate_file, level) {
    // Get the units
    Tnorm = get<BoutReal>(units["eV"]);
    Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
    FreqNorm = 1. / get<BoutReal>(units["seconds"]);
    substitutePermissions("val", {"density", "temperature", "velocity"});
    substitutePermissions("e_val", {"density", "temperature"});
    substitutePermissions("w_val",
                          {"density_source", "momentum_source", "energy_source"});
    substitutePermissions("sp", {from_A, from_B, to_A, to_B});
    substitutePermissions("from_ion", {from_A, from_B});
  }
  /// Perform charge exchange
  ///
  /// from_A + from_B -> to_A + to_B
  ///
  /// from_A and to_A must have the same atomic mass
  /// from_B and to_B must have the same atomic mass
  /// The charge of from_A + from_B must equal the charge of to_A + to_B
  void calculate_rates(GuardedOptions && electron, GuardedOptions && from_A, GuardedOptions && from_B, GuardedOptions && to_A,
                       GuardedOptions && to_B);

private:
  OpenADASRateCoefficient rate_coef;      ///< Reaction rate coefficient
  BoutReal Tnorm, Nnorm, FreqNorm; ///< Normalisations
};

#endif // ADAS_REACTION_H
