#pragma once
#ifndef AMJUEL_REACTION_H
#define AMJUEL_REACTION_H

#include <algorithm>
#include <filesystem>
#include <string>

#include "amjueldata.hxx"
#include "component.hxx"
#include "integrate.hxx"
#include "reaction.hxx"

/**
 * @brief Extract the json database location from options, or fall back to a standard
 * location in the repo (relative to this header).
 *
 * @param options Options object
 * @return std::filesystem::path the path to be passed to AmjuelData
 */
static inline std::filesystem::path get_json_db_dir(Options& options) {
  static std::filesystem::path default_json_db_dir =
      std::filesystem::path(__FILE__).parent_path().parent_path() / "json_database";

  std::string json_db_dir =
      options["json_database_dir"]
          .doc("Path to directory containing reaction data json files.")
          .withDefault(default_json_db_dir.string());

  return std::filesystem::path(json_db_dir);
}

/**
 * @brief Base class for reactions that read data using AmjuelData.
 *
 */
struct AmjuelReaction : public Reaction {
  AmjuelReaction(std::string name, std::string short_reaction_type,
                 std::string amjuel_lbl, std::string from_species, std::string to_species,
                 Options& alloptions)
      : Reaction(name, alloptions),
        amjuel_src(std::string("Amjuel ") + amjuel_lbl),
        short_reaction_type(short_reaction_type), from_species(from_species),
        to_species(to_species), amjuel_data(get_json_db_dir(alloptions), short_reaction_type, amjuel_lbl) {

    this->includes_sigma_v_e = amjuel_data.includes_sigma_v_e;
  }

protected:
  // Store some strings for use in attribute docstrings
  const std::string amjuel_src;
  const std::string short_reaction_type;
  const std::string from_species;
  const std::string to_species;

  /// Functions to calculate Amjuel rates from underlying tables
  BoutReal eval_amjuel_fit(BoutReal T, BoutReal n,
                           const std::vector<std::vector<BoutReal>>& coeff_table);
  virtual BoutReal eval_sigma_v_E(BoutReal T, BoutReal n) override final;
  virtual BoutReal eval_sigma_v(BoutReal T, BoutReal n) override final;

  virtual void transform_additional(Options& state,
                                    Field3D& reaction_rate) override final;

private:
  /// Data object
  const AmjuelData amjuel_data;

  /// Directory containing json data for this reaction
  std::filesystem::path json_db_dir;
};

#endif // AMJUEL_REACTION_H
