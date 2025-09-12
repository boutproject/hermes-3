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
  friend class AmjuelReaction;

private:
  AmjuelData(const std::filesystem::path& data_dir,
             const std::string& short_reaction_type, const std::string& data_label) {
    AUTO_TRACE();

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

    try {
      includes_sigma_v_e = data["info"]["includes_sigma_v_e"];
    } catch (nlohmann::json::type_error e) {
      throw BoutException(fmt::format("Amjuel json data at '{:s}' doesn't contain key "
                                      "'info/includes_sigma_v_e'. Error was {:s}",
                                      file_path.string(), e.what()));
    }

    try {
      // Extract <sigma v> coeff table
      std::vector<std::vector<double>> sigma_v_coeffs_tmp = data["sigma_v_coeffs"];
      this->sigma_v_coeffs = sigma_v_coeffs_tmp;
      if (this->includes_sigma_v_e) {
        // Extract <sigma v E> coeff table
        std::vector<std::vector<double>> sigma_v_E_coeffs_tmp = data["sigma_v_E_coeffs"];
        this->sigma_v_E_coeffs = sigma_v_E_coeffs_tmp;

        // Extract electron heating value
        double electron_heating_tmp = data["electron_heating"];
        this->electron_heating = electron_heating_tmp;
      }
    } catch (nlohmann::json::type_error e) {
      throw BoutException(fmt::format("json file at '{:s}' doesn't contain valid Amjuel "
                                      "coefficient data. Error was {:s}",
                                      file_path.string(), e.what()));
    }
  }

  // N.B. E-index varies fastest, so coefficient indices are [T][n]
  std::vector<std::vector<BoutReal>> sigma_v_coeffs;
  std::vector<std::vector<BoutReal>> sigma_v_E_coeffs;

  BoutReal electron_heating;

  bool includes_sigma_v_e;
};

#endif