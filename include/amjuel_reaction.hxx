#pragma once
#ifndef AMJUEL_REACTION_H
#define AMJUEL_REACTION_H

#include "../external/json.hxx"
#include "component.hxx"
#include "integrate.hxx"
#include "reaction.hxx"
#include <filesystem>
#include <string>

struct AmjuelData {

  AmjuelData(const std::string& amjuel_label) {
    AUTO_TRACE();

    std::filesystem::path filename =
        std::filesystem::path(__FILE__).parent_path().parent_path() / "json_database"
        / (std::string(amjuel_label) + ".json");

    // Read the data file
    std::ifstream json_file(filename);

    if (!json_file.good()) {
      throw BoutException("Could not read Amjuel data file '{}'", std::string(filename));
    }

    nlohmann::json data;
    json_file >> data;

    // Extract coeffs and copy them to the member vars
    std::vector<std::vector<double>> rate_coeffs_tmp = data["rate_coeffs"];
    this->rate_coeffs = rate_coeffs_tmp;
    std::vector<std::vector<double>> rad_coeffs_tmp = data["radiation_coeffs"];
    this->rad_coeffs = rad_coeffs_tmp;

    // Extract electron heating value
    double electron_heating_tmp = data["electron_heating"];
    this->electron_heating = electron_heating_tmp;
  }

  // N.B. E-index varies fastest, so coefficient is [T][n]
  std::vector<std::vector<BoutReal>> rate_coeffs;
  std::vector<std::vector<BoutReal>> rad_coeffs;

  BoutReal electron_heating;

  // BoutReal Tmin, Tmax; ///< Range of T  [eV]
  // BoutReal nmin, nmax; ///< Range of density [m^-3]
};

struct AmjuelReaction : public Reaction {
  AmjuelReaction(std::string name, std::string amjuel_label, std::string from_species,
                 std::string to_species, Options& alloptions)
      : Reaction(name, alloptions), amjuel_data(amjuel_label),
        amjuel_src(std::string("amjuel_") + amjuel_label), from_species(from_species),
        to_species(to_species) {}

protected:
  const AmjuelData amjuel_data;

  const std::string amjuel_src;
  const std::string from_species;
  const std::string to_species;

  // For diagnostics
  Field3D S; ///< Particle exchange
  Field3D F; ///< Momentum exchange
  Field3D E; ///< Energy exchange
  Field3D R; ///< Radiation loss

  const BoutReal get_electron_heating() const { return amjuel_data.electron_heating; }

  const std::vector<std::vector<BoutReal>>& get_rad_coeffs() const {
    return amjuel_data.rad_coeffs;
  }
  const std::vector<std::vector<BoutReal>>& get_rate_coeffs() const {
    return amjuel_data.rate_coeffs;
  }

  BoutReal clip(BoutReal value, BoutReal min, BoutReal max) {
    if (value < min)
      return min;
    if (value > max)
      return max;
    return value;
  }

  /**
   * @brief Evaluate an Amjuel rate, given a table of polynomial fit coefficients
   *
   *
   * @param T temperature
   * @param n number density
   * @param coeff_table a table of polynomial fit coefficients
   * @return BoutReal the rate
   */
  BoutReal eval_rate(BoutReal T, BoutReal n,
                     const std::vector<std::vector<BoutReal>>& coeff_table) {
    // Enforce range of validity
    n = clip(n, 1e14, 1e22); // 1e8 - 1e16 cm^-3
    T = clip(T, 0.1, 1e4);

    BoutReal logntilde = log(n / 1e14); // Note: 1e8 cm^-3
    BoutReal logT = log(T);
    BoutReal result = 0.0;

    BoutReal logT_n = 1.0; // log(T) ** n
    for (size_t n = 0; n < coeff_table.size(); ++n) {
      BoutReal logn_m = 1.0; // log(ntilde) ** m
      for (size_t m = 0; m < coeff_table[n].size(); ++m) {
        result += coeff_table[n][m] * logn_m * logT_n;
        logn_m *= logntilde;
      }
      logT_n *= logT;
    }
    return exp(result) * 1e-6; // Note: convert cm^3 to m^3
  }

  /**
   * @brief Evaluate electron energy loss rate as a function of temperature and density.
   *
   * @param T temperature
   * @param n number density
   * @return BoutReal the electron energy loss rate
   */
  virtual BoutReal eval_electron_energy_loss_rate(BoutReal T, BoutReal n) override final {
    return eval_rate(T, n, get_rad_coeffs());
  }

  /**
   * @brief Evaluate reaction rate as a function of temperature and density.
   *
   * @param T temperature
   * @param n number density
   * @return BoutReal the reaction rate
   */
  virtual BoutReal eval_reaction_rate(BoutReal T, BoutReal n) override final {
    return eval_rate(T, n, get_rate_coeffs());
  }

  virtual void transform_additional(Options& state, Field3D& reaction_rate,
                                    Field3D& momentum_exchange, Field3D& energy_exchange,
                                    Field3D& energy_loss) override final {

    Options& electron = state["species"]["e"];
    Field3D T_e = get<Field3D>(electron["temperature"]);

    std::vector<std::string> reactant_species =
        parser->get_species(species_filter::reactants);

    // Restrict to 2 reactants for now;
    ASSERT1(reactant_species.size() == 2);
    Options& r1 = state["species"][reactant_species[0]];
    Options& r2 = state["species"][reactant_species[1]];
    Field3D n_r1 = get<Field3D>(r1["density"]);
    Field3D n_r2 = get<Field3D>(r2["density"]);

    // Get heavy reactant species
    std::vector<std::string> heavy_reactant_species;
    std::copy_if(reactant_species.begin(), reactant_species.end(),
                 std::back_inserter(heavy_reactant_species),
                 [](std::string sp_name) { return sp_name != "e"; });
    Options& rh = state["species"][heavy_reactant_species[0]];
    Field3D n_rh = get<Field3D>(rh["density"]);

    // Get heavy product species
    std::vector<std::string> heavy_product_species =
        parser->get_species(species_filter::heavy, species_filter::products);

    // Get the velocity of the heavy product
    Options& ph = state["species"][heavy_product_species[0]];

    // Energy source for electrons (Does this need to be separate?)
    const int e_pop_change = this->parser->get_stoich().at("e");
    if (e_pop_change != 0) {
      ASSERT1(electron.isSet("velocity"));
      auto v_e = get<Field3D>(electron["velocity"]);
      auto m_e = get<BoutReal>(electron["AA"]);
      add(electron["energy_source"], 0.5 * m_e * e_pop_change * reaction_rate * SQ(v_e));
    }

    // Electron energy loss (radiation, ionisation potential)
    Field3D n_e = get<Field3D>(electron["density"]);
    energy_loss = cellAverage(
        [&](BoutReal nr1, BoutReal nr2, BoutReal ne, BoutReal te) {
          return nr1 * nr2 * eval_electron_energy_loss_rate(te * Tnorm, ne * Nnorm)
                 * Nnorm / (Tnorm * FreqNorm) * radiation_multiplier;
        },
        n_e.getRegion("RGN_NOBNDRY"))(n_r1, n_r2, n_e, T_e);

    // Loss is reduced by heating
    energy_loss -=
        (get_electron_heating() / Tnorm) * reaction_rate * radiation_multiplier;

    subtract(electron["energy_source"], energy_loss);

    // Collision frequencies

    // Same as reaction_rate but without the n1 factor: returns atom/ion collisionality in
    // [s^-1]
    Field3D heavy_particle_frequency = cellAverage(
        [&](BoutReal ne, BoutReal te) {
          return ne * eval_reaction_rate(te * Tnorm, ne * Nnorm) * Nnorm / FreqNorm
                 * rate_multiplier;
        },
        n_e.getRegion("RGN_NOBNDRY"))(n_e, T_e);

    // Same as reaction_rate but without the ne factor: returns electron collisionality in
    // [s^-1]
    Field3D electron_frequency = cellAverage(
        [&](BoutReal ne, BoutReal n1, BoutReal te) {
          return n1 * eval_reaction_rate(te * Tnorm, ne * Nnorm) * Nnorm / FreqNorm
                 * rate_multiplier;
        },
        n_e.getRegion("RGN_NOBNDRY"))(n_e, n_rh, T_e);

    // Add individual reaction collision frequency to each species
    std::string reaction_type;
    if (rh.name().find("+") != std::string::npos) {
      reaction_type = "rec";
    } else {
      reaction_type = "iz";
    }
    if (reaction_type == "iz") {
      set(rh["collision_frequencies"]
            [rh.name() + std::string("_") + ph.name() + std::string("_iz")],
          heavy_particle_frequency);
    } else if (reaction_type == "rec") {
      set(ph["collision_frequencies"]
            [rh.name() + std::string("_") + ph.name() + std::string("_rec")],
          heavy_particle_frequency);
    }
  }
};

#endif // AMJUEL_REACTION_H
