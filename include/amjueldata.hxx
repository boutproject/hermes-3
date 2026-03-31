#pragma once
#ifndef AMJUELDATA_H
#define AMJUELDATA_H

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "../external/json.hxx"
#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/msg_stack.hxx>

/**
 * @brief Handle reading and storage of Amjuel reaction data.
 *
 */
struct AmjuelData {
  friend struct AmjuelReaction;

private:
  AmjuelData(const std::filesystem::path& data_dir,
             const std::string& short_reaction_type, const std::string& data_label) {

    if (!std::filesystem::is_directory(data_dir)) {
      throw BoutException(fmt::format("No json database found at ", data_dir.string()));
    }

    std::filesystem::path file_path =
        data_dir / (short_reaction_type + "_AMJUEL_" + data_label + ".json");

    // Read the data file
    std::ifstream json_file(file_path);

    if (!json_file.good()) {
      throw BoutException("Could not read Amjuel data file '{}'", std::string(file_path));
    }

    // Parse the data
    nlohmann::json data;
    json_file >> data;

    // Extract fit type (rate params)
    this->fit_type = data["info"]["fit_type"];

    try {
      includes_sigma_v_e = data["info"]["includes_sigma_v_e"];
    } catch (nlohmann::json::type_error& e) {
      throw BoutException(fmt::format("Amjuel json data at '{:s}' doesn't contain key "
                                      "'info/includes_sigma_v_e'. Error was {:s}",
                                      file_path.string(), e.what()));
    }

    try {
      // Extract <sigma v> coeff table
      this->sigma_v_coeffs = data["sigma_v_coeffs"];
      if (this->includes_sigma_v_e) {
        // Extract <sigma v E> coeff table
        this->sigma_v_E_coeffs = data["sigma_v_E_coeffs"];
        // Extract electron heating value
        this->electron_heating = data["electron_heating"];
      }
    } catch (nlohmann::json::type_error& e) {
      throw BoutException(fmt::format("json file at '{:s}' doesn't contain valid Amjuel "
                                      "coefficient data. Error was {:s}",
                                      file_path.string(), e.what()));
    }
  }

  std::string fit_type;

  // N.B. E-index varies fastest, so coefficient indices are [T][n]
  std::vector<std::vector<BoutReal>> sigma_v_coeffs;
  std::vector<std::vector<BoutReal>> sigma_v_E_coeffs;

  BoutReal electron_heating;

  bool includes_sigma_v_e;
};

#endif
