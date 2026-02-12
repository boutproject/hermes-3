#include "amjuel_reaction.hxx"
#include "hermes_utils.hxx"

/**
 * @brief Evaluate an Amjuel double polynomial fit in n and T, given a table of
 * coefficients (see page 20 of amjuel.pdf).
 *
 * @param T temperature in eV
 * @param n number density in m^-3
 * @param coeff_table polynomial fit coefficients in a vector-of-vectors (outer index T,
 * inner index n)
 * @return BoutReal the fit in SI, units m^3/s, or eV m^3/s for energy loss
 */
BoutReal AmjuelReaction::eval_amjuel_nT_fit(
    BoutReal T, BoutReal n, const std::vector<std::vector<BoutReal>>& coeff_table) {
  // Enforce range of validity
  n = std::clamp(n, 1e14, 1e22); // 1e8 - 1e16 cm^-3
  T = std::clamp(T, 0.1, 1e4);

  BoutReal logntilde = log(n / 1e14); // Note: 1e8 cm^-3
  BoutReal logT = log(T);
  BoutReal result = 0.0;

  BoutReal logT_n = 1.0; // log(T) ** n
  std::size_t num_n = coeff_table.size();
  std::size_t num_T = coeff_table[0].size();
  for (size_t n_idx = 0; n_idx < num_n; ++n_idx) {
    BoutReal logn_m = 1.0; // log(ntilde) ** m
    const auto& nrow = coeff_table[n_idx];
    for (size_t T_idx = 0; T_idx < num_T; ++T_idx) {
      result += nrow[T_idx] * logn_m * logT_n;
      logn_m *= logntilde;
    }
    logT_n *= logT;
  }
  return exp(result) * 1e-6; // Note: convert cm^3 to m^3
}

/**
 * @brief Evaluate an Amjuel single polynomial fit in T, given a table of
 * coefficients (see page 20 of amjuel.pdf).
 *
 * @param T temperature in eV
 * @param coeff_table a table of polynomial fit coefficients (index T)
 * @return BoutReal the fit in SI, units m^3/s, or eV m^3/s for energy loss
 */
BoutReal AmjuelReaction::eval_amjuel_T_fit(BoutReal T,
                                           const std::vector<BoutReal>& coeff_table) {
  const BoutReal lnT = log(T);
  BoutReal ln_sigmav = coeff_table[0];
  BoutReal lnT_n = lnT; // (lnT)^n
  const std::size_t num_T = coeff_table.size();
  for (std::size_t T_idx = 1; T_idx < num_T; T_idx++) {
    ln_sigmav += coeff_table[T_idx] * lnT_n;
    lnT_n *= lnT;
  }

  return exp(ln_sigmav);
}

/**
 * @brief Use Amjuel coeffs to evaluate <sigma.v.E> at a particular density and
 * temperature.
 *
 * @param T temperature
 * @param n number density
 * @return BoutReal <sigma.v.E>(n, T)
 */
BoutReal AmjuelReaction::eval_sigma_vE_nT(BoutReal T, BoutReal n) {
  return eval_amjuel_nT_fit(T, n, amjuel_data.sigma_v_E_coeffs);
}

/**
 * @brief Use Amjuel coeffs to evaluate <sigma.v> at a particular density and
 * temperature.
 *
 * @param T a temperature
 * @param n a density
 * @return BoutReal <sigma.v>(n, T)
 */
BoutReal AmjuelReaction::eval_sigma_v_nT(BoutReal T, BoutReal n) {
  return eval_amjuel_nT_fit(T, n, amjuel_data.sigma_v_coeffs);
}

/**
 * @brief Use Amjuel coeffs to evaluate <sigma.v> at a particular (effective)
 * temperature
 *
 * @param T a temperature
 * @return BoutReal <sigma.v>(T_eff)
 */
BoutReal AmjuelReaction::eval_sigma_v_T(BoutReal T) {
  return eval_amjuel_T_fit(T, amjuel_data.sigma_v_coeffs[0]);
}

/**
 * @brief Extract rate parameters type from json data
 *
 * @return RateParamsTypes the rate parameters type
 */
RateParamsTypes AmjuelReaction::get_rate_params_type() const {
  std::string fit_type_lcase = this->amjuel_data.fit_type;
  std::transform(fit_type_lcase.begin(), fit_type_lcase.end(), fit_type_lcase.begin(),
                 ::tolower);
  return RateParamsTypesFromString(fit_type_lcase);
}

void AmjuelReaction::transform_additional(GuardedOptions& state,
                                          const RateData& rate_data) {

  // Extract the rate
  Field3D rate = rate_data.rate;

  // Amjuel-based reactions are assumed to have exactly 2 reactants, for now.
  std::vector<std::string> reactant_species =
      parser->get_species(species_filter::reactants);
  ASSERT1(reactant_species.size() == 2);

  // Amjuel-based reactions are assumed to have exactly 1 heavy reactant unless this
  // function has been overridden
  std::string heavy_reactant_species =
      parser->get_single_species(reactant_species, species_filter::heavy);
  GuardedOptions rh = state["species"][heavy_reactant_species];
  const BoutReal AA_rh = get<BoutReal>(rh["AA"]);
  const Field3D& n_rh = rh["density"].GetRef<Field3D>();
  const Field3D& v_rh = rh["velocity"].GetRef<Field3D>();

  // Amjuel-based reactions are assumed to have exactly 1 heavy product unless this
  // function has been overridden
  std::string heavy_product_species =
      parser->get_single_species(species_filter::heavy, species_filter::products);
  GuardedOptions ph = state["species"][heavy_product_species];
  const Field3D& v_ph = ph["velocity"].GetRef<Field3D>();

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
  add(ph["energy_source"], 0.5 * AA_rh * rate * SQ(v_rh - v_ph));

  // Energy source for electrons due to pop change
  GuardedOptions electron = state["species"]["e"];
  const Field3D& T_e = electron["temperature"].GetRef<Field3D>();
  const int e_pop_change = this->parser->pop_change("e");
  if (e_pop_change != 0) {
    if (electron.isSet("velocity")) {
      // Transfer of electron kinetic to thermal energy due to density source
      // For ionisation:
      // Electrons with zero average velocity are created, diluting the kinetic energy.
      // Total energy conservation requires a corresponding internal energy source.
      //
      // For recombination:
      // Electrons with some velocity are incorporated into a neutral species with their
      // kinetic energy converted to an internal energy source of that species.
      const Field3D& v_e = electron["velocity"].GetRef<Field3D>();
      const BoutReal m_e = get<BoutReal>(electron["AA"]);
      add(electron["energy_source"], 0.5 * m_e * e_pop_change * rate * SQ(v_e));
    }
  }

  // Electron energy loss (radiation, ionisation potential)
  const Field3D& n_e = electron["density"].GetRef<Field3D>();
  auto region_no_bndry = n_e.getRegion("RGN_NOBNDRY");

  Field3D energy_loss = cellAverage(
      [&](BoutReal nrh, BoutReal ne, BoutReal te) {
        return nrh * ne * eval_sigma_vE_nT(te * Tnorm, ne * Nnorm) * Nnorm
               / (Tnorm * FreqNorm) * radiation_multiplier;
      },
      region_no_bndry)(n_rh, n_e, T_e);

  // Loss is reduced by heating
  energy_loss -= (amjuel_data.electron_heating / Tnorm) * rate * radiation_multiplier;

  update_source<subtract<Field3D>>(state, "e", ReactionDiagnosticType::energy_loss,
                                   energy_loss);

  // Collision frequencies [s^-1]

  // Set on the neutral species for both ionisation and recombination
  std::string neutral_species = parser->get_single_species(species_filter::neutral);
  std::string heavy_cf_lbl =
      fmt::format("{:s}_{:s}_{:s}", heavy_reactant_species, heavy_product_species,
                  this->short_reaction_type);
  set(state["species"][neutral_species]["collision_frequencies"][heavy_cf_lbl],
      rate_data.coll_freq(rh.name()));

  // N.B. nothing is done with the electron collision frequency; presumably we want to
  // update it in the state?!
  Field3D electron_frequency = rate_data.coll_freq("e");
}
