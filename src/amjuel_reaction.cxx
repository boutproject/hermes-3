#include "amjuel_reaction.hxx"
#include "hermes_utils.hxx"

/**
 * @brief Evaluate electron energy loss rate as a function of temperature and density.
 *
 * @param T temperature
 * @param n number density
 * @return BoutReal the electron energy loss rate
 */
BoutReal AmjuelReaction::eval_electron_energy_loss_rate(BoutReal T, BoutReal n) {
  return eval_rate(T, n, get_rad_coeffs());
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
BoutReal
AmjuelReaction::eval_rate(BoutReal T, BoutReal n,
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
 * @brief Evaluate reaction rate as a function of temperature and density.
 *
 * @param T temperature
 * @param n number density
 * @return BoutReal the reaction rate
 */
BoutReal AmjuelReaction::eval_reaction_rate(BoutReal T, BoutReal n) {
  return eval_rate(T, n, get_rate_coeffs());
}

const BoutReal AmjuelReaction::get_electron_heating() const {
  return amjuel_data.electron_heating;
}

const std::vector<std::vector<BoutReal>>& AmjuelReaction::get_rad_coeffs() const {
  return amjuel_data.rad_coeffs;
}

const std::vector<std::vector<BoutReal>>& AmjuelReaction::get_rate_coeffs() const {
  return amjuel_data.rate_coeffs;
}

void AmjuelReaction::transform_additional(Options& state, Field3D& reaction_rate,
                                          Field3D& momentum_exchange,
                                          Field3D& energy_exchange,
                                          Field3D& energy_loss) {

  // Amjuel-based reactions are assumed to have 2 reactants, for now;
  std::vector<std::string> reactant_species =
      parser->get_species(species_filter::reactants);
  ASSERT1(reactant_species.size() == 2);

  // Extract heavy reactant properties
  std::vector<std::string> heavy_reactant_species;
  std::copy_if(reactant_species.begin(), reactant_species.end(),
               std::back_inserter(heavy_reactant_species),
               [](std::string sp_name) { return sp_name != "e"; });
  Options& rh = state["species"][heavy_reactant_species[0]];
  BoutReal AA_rh = get<BoutReal>(rh["AA"]);
  Field3D n_rh = get<Field3D>(rh["density"]);
  Field3D v_rh = get<Field3D>(rh["velocity"]);

  // Extract heavy product properties
  std::vector<std::string> heavy_product_species =
      parser->get_species(species_filter::heavy, species_filter::products);
  Options& ph = state["species"][heavy_product_species[0]];
  Field3D v_ph = get<Field3D>(ph["velocity"]);

  // Energy source for heavy product
  add(ph["energy_source"], 0.5 * AA_rh * reaction_rate * SQ(v_rh - v_ph));

  // Energy source for electrons due to pop change
  Options& electron = state["species"]["e"];
  Field3D T_e = get<Field3D>(electron["temperature"]);
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
      [&](BoutReal nrh, BoutReal ne, BoutReal te) {
        return nrh * ne * eval_electron_energy_loss_rate(te * Tnorm, ne * Nnorm) * Nnorm
               / (Tnorm * FreqNorm) * radiation_multiplier;
      },
      n_e.getRegion("RGN_NOBNDRY"))(n_rh, n_e, T_e);

  // Loss is reduced by heating
  energy_loss -= (get_electron_heating() / Tnorm) * reaction_rate * radiation_multiplier;

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
  set(rh["collision_frequencies"]
        [rh.name() + "_" + ph.name() + "_" + this->short_reaction_type],
      heavy_particle_frequency);
}