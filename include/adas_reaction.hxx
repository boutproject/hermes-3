#pragma once
#ifndef ADAS_REACTION_H
#define ADAS_REACTION_H

#include "component.hxx"
#include <vector>

struct OpenADASRateCoefficient {
  /// Read the file, extracting data for the given ionisation level
  OpenADASRateCoefficient(const std::string& filename, int level);

  std::vector<std::vector<BoutReal>> log_coeff;
  std::vector<BoutReal> log_temperature;
  std::vector<BoutReal> log_density;

  BoutReal Tmin, Tmax; // Range of T
  BoutReal nmin, nmax; // Range of density

  /// Input in units:
  ///     n in m^-3
  ///     T in eV
  ///
  /// Output in units m^3/s or eV m^3/s
  BoutReal evaluate(BoutReal T, BoutReal n);
};

/// Read in and perform calculations with OpenADAS data
/// https://open.adas.ac.uk/
///
/// Uses the JSON files produced by:
///   https://github.com/TBody/OpenADAS_to_JSON
struct OpenADAS : public Component {
  ///
  /// Inputs
  /// ------
  ///   units       Options tree containing normalisation constants
  ///   rate_file   A JSON file containing reaction rate <σv> rates (e.g. SCD, ACD)
  ///   radiation_file   A JSON file containing radiation loss rates (e.g. PLT, PRB)
  ///   level       The lower ionisation state in the transition
  ///               e.g. 0 for neutral -> 1st ionisation
  ///               and 1st -> neutral recombination
  ///   electron_heating   The heating of the electrons per reaction [eV]
  ///               This is the ionisation energy, positive for recombination
  ///               and negative for ionisation
  ///
  /// Notes
  ///  - The rate and radiation file names have "json_database/" prepended
  /// 
  OpenADAS(const Options& units, const std::string& rate_file,
           const std::string& radiation_file, int level, BoutReal electron_heating)
      : rate_coef(std::string("json_database/") + rate_file, level),
        radiation_coef(std::string("json_database/") + radiation_file, level),
        electron_heating(electron_heating) {
    // Get the units
    Tnorm = get<BoutReal>(units["eV"]);
    Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
    FreqNorm = 1. / get<BoutReal>(units["seconds"]);
  }

  /// Perform the calculation of rates, and transfer of particles/momentum/energy
  void calculate_rates(Options& electron, Options& from_ion, Options& to_ion);

private:
  OpenADASRateCoefficient rate_coef;      ///< Reaction rate coefficient
  OpenADASRateCoefficient radiation_coef; ///< Energy loss (radiation) coefficient

  BoutReal electron_heating; // Heating per reaction [eV]

  BoutReal Tnorm, Nnorm, FreqNorm; // Normalisations
};

#endif // ADAS_REACTION_H
