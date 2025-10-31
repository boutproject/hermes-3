#include "amjuel_reaction.hxx"
#include "hermes_utils.hxx"

/**
 * @brief Evaluate an Amjuel double polynomial fit in n and T, given a table of
 * coefficients (see page 20 of amjuel.pdf).
 *
 * @param T temperature in eV
 * @param n number density in m^-3
 * @param coeff_table a table of polynomial fit coefficients (coefs[T][n])
 * @return BoutReal the fit in SI, units m^3/s, or eV m^3/s for energy loss
 */
BoutReal
AmjuelReaction::eval_amjuel_fit(BoutReal T, BoutReal n,
                                const std::vector<std::vector<BoutReal>>& coeff_table) {
  // Enforce range of validity
  n = std::clamp(n, 1e14, 1e22); // 1e8 - 1e16 cm^-3
  T = std::clamp(T, 0.1, 1e4);

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
 * @brief Evaluate <sigma . v . E> at a particular density and temperature by evaluating
 * an Amjuel fit.
 *
 * @param T temperature
 * @param n number density
 * @return BoutReal <sigma . v . E>(n, T)
 */
BoutReal AmjuelReaction::eval_sigma_v_E(BoutReal T, BoutReal n) {
  return eval_amjuel_fit(T, n, amjuel_data.sigma_v_E_coeffs);
}

/**
 * @brief Evaluate <sigma . v . E> at a particular density and temperature
 * (Subclasses MAY define)
 *
 * @param T a temperature
 * @param n a density
 * @return BoutReal <sigma . v . E>(n, T)
 */
BoutReal AmjuelReaction::eval_sigma_v(BoutReal T, BoutReal n) {
  return eval_amjuel_fit(T, n, amjuel_data.sigma_v_coeffs);
}

void AmjuelReaction::transform_additional(GuardedOptions& state, Field3D& reaction_rate) {

  // Amjuel-based reactions are assumed to have exactly 2 reactants, for now.
  std::vector<std::string> reactant_species =
      parser->get_species(species_filter::reactants);
  ASSERT1(reactant_species.size() == 2);

  // Extract heavy reactant properties
  std::vector<std::string> heavy_reactant_species =
      parser->get_species(reactant_species, species_filter::heavy);
  // Amjuel-based reactions are assumed to have exactly 1 heavy reactant, for now.
  ASSERT1(heavy_reactant_species.size() == 1);
  GuardedOptions rh = state["species"][heavy_reactant_species[0]];
  BoutReal AA_rh = get<BoutReal>(rh["AA"]);
  Field3D n_rh = get<Field3D>(rh["density"]);
  Field3D v_rh = get<Field3D>(rh["velocity"]);

  // Extract heavy product properties
  std::vector<std::string> heavy_product_species =
      parser->get_species(species_filter::heavy, species_filter::products);
  // Amjuel-based reactions are assumed to have exactly 1 heavy product, for now.
  ASSERT1(heavy_product_species.size() == 1);
  GuardedOptions ph = state["species"][heavy_product_species[0]];
  Field3D v_ph = get<Field3D>(ph["velocity"]);

  // Kinetic energy transfer to thermal energy
  //
  // This is a combination of three terms in the pressure
  // (thermal energy) equation:
  //
  // d/dt(3/2 p_1) = ... - F_12 v_1 + W_12 - (1/2) m R v_1^2 // From ion ('heavy
  // reactant')
  //
  // d/dt(3/2 p_2) = ... - F_21 v_2 + W_21 + (1/2) m R v_2^2 // to_ion ('heavy product')
  //
  // where forces are F_21 = -F_12, energy transfer W_21 = -W_12,
  // and masses are equal so m_1 = m_2 = m
  //
  // Here the force F_21 = m R v_1   i.e. momentum_exchange
  //
  // As p_1 -> 0 the sum of the p_1 terms must go to zero.
  // Hence:
  //
  // W_12 = F_12 v_1 + (1/2) m R v_1^2
  //      = -R m v_1^2 + (1/2) m R v_1^2
  //      = - (1/2) m R v_1^2
  //
  // d/dt(3/2 p_2) = - m R v_1 v_2 + (1/2) m R v_1^2 + (1/2) m R v_2^2
  //               = (1/2) m R (v_1 - v_2)^2
  //
  //
  // This term accounts for the broadening of the species velocity
  // distribution when combining the two underlying distributions.
  // The greater the difference in velocities of the underlying species,
  // the wider the resultant distribution, which corresponds to an
  // increase in temperature and therefore internal energy.
  add(ph["energy_source"], 0.5 * AA_rh * reaction_rate * SQ(v_rh - v_ph));

  // Energy source for electrons due to pop change
  GuardedOptions electron = state["species"]["e"];
  Field3D T_e = get<Field3D>(electron["temperature"]);
  const int e_pop_change = this->parser->get_stoich().at("e");
  if (e_pop_change != 0) {
    if (IS_SET(electron["velocity"])) {
      // Transfer of electron kinetic to thermal energy due to density source
      // For ionisation:
      // Electrons with zero average velocity are created, diluting the kinetic energy.
      // Total energy conservation requires a corresponding internal energy source.
      //
      // For recombination:
      // Electrons with some velocity are incorporated into a neutral species with their
      // kinetic energy converted to an internal energy source of that species.
      auto v_e = get<Field3D>(electron["velocity"]);
      auto m_e = get<BoutReal>(electron["AA"]);
      add(electron["energy_source"], 0.5 * m_e * e_pop_change * reaction_rate * SQ(v_e));
    }
  }

  // Electron energy loss (radiation, ionisation potential)
  Field3D n_e = get<Field3D>(electron["density"]);
  Field3D energy_loss = cellAverage(
      [&](BoutReal nrh, BoutReal ne, BoutReal te) {
        return nrh * ne * eval_sigma_v_E(te * Tnorm, ne * Nnorm) * Nnorm
               / (Tnorm * FreqNorm) * radiation_multiplier;
      },
      n_e.getRegion("RGN_NOBNDRY"))(n_rh, n_e, T_e);

  // Loss is reduced by heating
  energy_loss -=
      (amjuel_data.electron_heating / Tnorm) * reaction_rate * radiation_multiplier;

  update_source<subtract<Field3D>>(state, "e", ReactionDiagnosticType::energy_loss,
                                   energy_loss);

  // Collision frequencies

  // Same as reaction_rate but without the n1 factor: returns atom/ion collisionality in
  // [s^-1]
  Field3D heavy_particle_frequency = cellAverage(
      [&](BoutReal ne, BoutReal te) {
        return ne * eval_sigma_v(te * Tnorm, ne * Nnorm) * Nnorm / FreqNorm
               * rate_multiplier;
      },
      n_e.getRegion("RGN_NOBNDRY"))(n_e, T_e);

  // Same as reaction_rate but without the ne factor: returns electron collisionality in
  // [s^-1]
  Field3D electron_frequency = cellAverage(
      [&](BoutReal ne, BoutReal n1, BoutReal te) {
        return n1 * eval_sigma_v(te * Tnorm, ne * Nnorm) * Nnorm / FreqNorm
               * rate_multiplier;
      },
      n_e.getRegion("RGN_NOBNDRY"))(n_e, n_rh, T_e);

  // Set collision frequency on the neutral species (should be exactly 1 of them for
  // Amjuel reactions)
  std::vector<std::string> neutral_species = parser->get_species(species_filter::neutral);
  ASSERT1(neutral_species.size() == 1);
  set(state["species"][neutral_species[0]]["collision_frequencies"]
           [heavy_reactant_species[0] + "_" + heavy_product_species[0] + "_" + this->short_reaction_type],
      heavy_particle_frequency);
}
