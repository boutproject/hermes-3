#pragma once
#ifndef AMJUEL_REACTION_H
#define AMJUEL_REACTION_H

#include <filesystem>
#include <string>

#include "amjueldata.hxx"
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
                 std::string amjuel_lbl, Options& alloptions)
      : Reaction(name, alloptions), amjuel_src(std::string("Amjuel ") + amjuel_lbl),
        short_reaction_type(short_reaction_type),
        amjuel_data(get_json_db_dir(alloptions), short_reaction_type, amjuel_lbl) {

    this->includes_sigma_v_e = amjuel_data.includes_sigma_v_e;
    // Most of the access information we need is inherited from the parent Reaction class.
    // The electron velocity will be read if it is set
    setPermissions(readIfSet("species:e:velocity"));
    // The energy source is set for electrons
    setPermissions(readWrite("species:e:energy_source"));
    std::string heavy_reactant = this->parser->get_species(species_filter::reactants,
                                                           species_filter::heavy)[0],
                heavy_product = this->parser->get_species(species_filter::products,
                                                          species_filter::heavy)[0],
                neutral = this->parser->get_species(species_filter::neutral)[0];
    setPermissions(
        readWrite(fmt::format("species:{}:collision_frequencies:{}_{}_{}", neutral,
                              heavy_reactant, heavy_product, short_reaction_type)));
  }

protected:
  // Store some strings for use in attribute docstrings
  const std::string amjuel_src;
  const std::string short_reaction_type;

  /// Functions to calculate Amjuel rates from underlying tables
  BoutReal eval_amjuel_nT_fit(BoutReal T, BoutReal n,
                              const std::vector<std::vector<BoutReal>>& coeff_table);
  BoutReal eval_amjuel_T_fit(BoutReal T, const std::vector<BoutReal>& coeff_table);

  /// Override the various rate-evaluator functions
  BoutReal eval_sigma_vE_nT(BoutReal T, BoutReal n) final;
  BoutReal eval_sigma_v_nT(BoutReal T, BoutReal n) final;
  BoutReal eval_sigma_v_T(BoutReal T) final;

  RateParamsTypes get_rate_params_type() const final;

  void transform_additional(GuardedOptions& state,
                            const RateData& rate_calc_results) override;

private:
  /// Data object
  const AmjuelData amjuel_data;

  /// Directory containing json data for this reaction
  std::filesystem::path json_db_dir;
};

#endif // AMJUEL_REACTION_H
