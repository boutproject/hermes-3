#pragma once
#ifndef AMJUEL_REACTION_H
#define AMJUEL_REACTION_H

#include <filesystem>
#include <string>

#include "../external/json.hxx"
#include "component.hxx"
#include "integrate.hxx"
#include "reaction.hxx"

/**
 * @brief Handle reading and storage of Amjuel reaction data.
 *
 */
struct AmjuelData {
  friend class AmjuelReaction;

private:
  AmjuelData(const std::string& short_reaction_type, const std::string& data_label) {
    AUTO_TRACE();

    std::filesystem::path file_path =
        std::filesystem::path(__FILE__).parent_path().parent_path() / "json_database"
        / (short_reaction_type + "_AMJUEL_" + data_label + ".json");

    // Read the data file
    std::ifstream json_file(file_path);

    if (!json_file.good()) {
      throw BoutException("Could not read Amjuel data file '{}'", std::string(file_path));
    }

    // Parse the data
    nlohmann::json data;
    json_file >> data;

    // Extract coeff tables into member vars
    std::vector<std::vector<double>> rate_coeffs_tmp = data["rate_coeffs"];
    this->rate_coeffs = rate_coeffs_tmp;
    std::vector<std::vector<double>> rad_coeffs_tmp = data["radiation_coeffs"];
    this->rad_coeffs = rad_coeffs_tmp;

    // Extract electron heating value into member var
    double electron_heating_tmp = data["electron_heating"];
    this->electron_heating = electron_heating_tmp;
  }

  // N.B. E-index varies fastest, so coefficient indices are [T][n]
  std::vector<std::vector<BoutReal>> rate_coeffs;
  std::vector<std::vector<BoutReal>> rad_coeffs;

  BoutReal electron_heating;
};

struct AmjuelReaction : public Reaction {
  AmjuelReaction(std::string name, std::string short_reaction_type,
                 std::string amjuel_lbl, std::string from_species, std::string to_species,
                 Options& alloptions)
      : Reaction(name, alloptions), amjuel_data(short_reaction_type, amjuel_lbl),
        amjuel_src(std::string("Amjuel ") + amjuel_lbl), from_species(from_species),
        short_reaction_type(short_reaction_type), to_species(to_species) {}

protected:
  // Store some strings for use in attribute docstrings
  const std::string amjuel_src;
  const std::string short_reaction_type;
  const std::string from_species;
  const std::string to_species;

  // For diagnostics
  Field3D S; ///< Particle exchange
  Field3D F; ///< Momentum exchange
  Field3D E; ///< Energy exchange
  Field3D R; ///< Radiation loss

  /// Functions to calculate Amjuel rates from underlying tables
  virtual BoutReal eval_electron_energy_loss_rate(BoutReal T, BoutReal n) override final;
  BoutReal eval_rate(BoutReal T, BoutReal n,
                     const std::vector<std::vector<BoutReal>>& coeff_table);
  virtual BoutReal eval_reaction_rate(BoutReal T, BoutReal n) override final;

  // Getter functions so that data object can be private
  const BoutReal get_electron_heating() const;
  const std::vector<std::vector<BoutReal>>& get_rad_coeffs() const;
  const std::vector<std::vector<BoutReal>>& get_rate_coeffs() const;

  virtual void transform_additional(Options& state, Field3D& reaction_rate,
                                    Field3D& momentum_exchange, Field3D& energy_exchange,
                                    Field3D& energy_loss) override final;

private:
  const AmjuelData amjuel_data;
};

#endif // AMJUEL_REACTION_H
