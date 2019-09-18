#pragma once

#ifndef IMPURITYSPECIES_H // Preprocessor directives to prevent multiple
                          // definitions
#define IMPURITYSPECIES_H

#include <map>
#include <string>
#include "RateCoefficient.hxx"

#include <memory>

class ImpuritySpecies {
  // # For storing OpenADAS data related to a particular impurity species
  // # Loosely based on cfe316/atomic/atomic_data.py/AtomicData class (although
  // with much less code since
  // # all of the F77 importing is done in the seperate <<make json_update>>
  // code since BOUT++ protocol
  // # requires fortran code be isolated from main operation)
public:
  /**
   * @brief ImpuritySpecies Constructor
   * @details Uses the "impurity_user_input" file to set hard-coded impurity
   * parameters corresponding to the impurity represented by the symbol
   * Determines the OpenADAS file-dependencies and adds them to .adas_data_dict
   * Makes RateCoefficient objects corresponding to the files and adds smart
   * pointers for them to .rate_coefficients
   *
   * @param impurity_symbol_supplied typically should be of length 1 (i.e. "c"
   * for carbon, "n" for nitrogen)
   */
  ImpuritySpecies(std::string &impurity_symbol_supplied);
  /**
   * @brief Determines the OpenADAS json files which the impurity data is given
   * in and adds them to .adas_data_dict
   * Will throw a runtime error if the file isn't found in the database
   *
   * @param physics_process a std::string corresponding to a physics process
   * @param filetype_code the code used by OpenADAS to represent this process
   * @param json_database_path relative or absolute path to where the json data
   * files from OpenADAS are located
   */
  void addJSONFiles(const std::string &physics_process, const std::string &filetype_code,
                    const std::string &json_database_path);
  /**
   * @brief Uses the OpenADAS files determined in addJSONFiles to construct a
   * map between a physics_process string and a smart-pointer to a corresponding
   * RateCoefficient object
   *
   */
  void makeRateCoefficients();
  std::string get_symbol();
  std::string get_name();
  int get_year();
  bool get_has_charge_exchange();
  int get_atomic_number();
  std::map<std::string, std::string> get_adas_files_dict();
  std::map<std::string, shared_ptr<RateCoefficient>> get_rate_coefficients();
  /**
   * @brief Accesses the value of the rate_coefficient map corresponding
   * to the supplied string key
   *
   * @param key string corresponding to a physics_process
   * @return shared (smart) pointer to a RateCoefficient object
   */
  shared_ptr<RateCoefficient> get_rate_coefficient(const std::string &key);

private:
  // Data fields
  std::string symbol;
  std::string name;
  int year;
  bool has_charge_exchange;
  int atomic_number;
  std::map<std::string, std::string> adas_files_dict;
  std::map<std::string, shared_ptr<RateCoefficient>> rate_coefficients;
};
std::string get_json_database_path();
std::string get_impurity_user_input();

#endif
