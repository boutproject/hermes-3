#include "izn_rec_reaction.hxx"
#include "hermes_utils.hxx"
#include "reaction.hxx"
#include "reaction_settings.hxx"

namespace hermes {

// ============================== IznRecReaction ==============================
///
IznRecReaction::IznRecReaction(std::string short_reaction_type, std::string name,
                               Options& options)
    : Reaction(name, options), short_reaction_type(short_reaction_type) {

  // Store heavy species names for convenience. Parser calls throw if there isn't exactly
  // one heavy reactant and one heavy product.
  this->heavy_reactant =
      this->parser->get_single_species(species_filter::heavy, species_filter::reactants);
  this->heavy_product =
      this->parser->get_single_species(species_filter::heavy, species_filter::products);

  /*
  Data for calculating the electron energy loss rate.
  Data type is assumed to be the same as the izn/rec rate data for now and the data ID is
  the default for that type.
  */
  ReactionDataTypes e_energy_loss_data_type = this->rate_data->get_type();

  const std::string e_energy_loss_data_src_id =
      get_default_data_id(this->parser->get_reaction_str(), e_energy_loss_data_type,
                          ReactionCoeffTypes::sigma_v_E);

  std::vector<std::string> metadata_keys = {"electron_heating"};
  this->e_energy_loss_data = ReactionDataFactory::getInstance().create(
      toString(e_energy_loss_data_type), e_energy_loss_data_src_id, options,
      metadata_keys);

  // Most of the access information we need is inherited from the parent Reaction
  // class. The electron velocity will be read if it is set
  setPermissions(readIfSet("species:e:velocity"));
  // The energy source is set for electrons
  setPermissions(readWrite("species:e:energy_source"));

  // Collision frequencies are keyed by the lower charge state heavy species (reactants
  // for ionisation, products for recombination)
  this->heavy_collfreq_species =
      (this->short_reaction_type == "iz") ? this->heavy_reactant : this->heavy_product;
  setPermissions(readWrite(fmt::format("species:{}:collision_frequencies:{}_{}_{}",
                                       this->heavy_collfreq_species, this->heavy_reactant,
                                       this->heavy_product, this->short_reaction_type)));

  if (this->diagnose) {
    // Set up diagnostics. Names and signs differ between ionisation and recombination.
    DiagnosticTransformerType default_transformer;
    std::string default_diag_suffix, rad_diag_suffix;
    // For izn, tweak R diagnostic name and reverse the signs of the S, F, E diagnostics
    if (this->short_reaction_type == "iz") {
      default_diag_suffix =
          fmt::format("{:s}_{:s}", this->heavy_product, short_reaction_type);
      default_transformer = identity;
      rad_diag_suffix = fmt::format("{:s}_ex", this->heavy_product);
    } else {
      default_diag_suffix =
          fmt::format("{:s}_{:s}", this->heavy_reactant, short_reaction_type);
      default_transformer = negate;
      rad_diag_suffix = default_diag_suffix;
    }
    DiagnosticTransformerType rad_transformer = identity;

    std::string long_reaction_type =
        this->short_reaction_type == "iz" ? "ionisation" : "recombination";
    std::string rate_data_src = this->rate_data->src_str();
    add_diagnostic(
        this->heavy_product, fmt::format("S{:s}", default_diag_suffix),
        fmt::format("Particle source due to {:s} of {:s} to {:s}", long_reaction_type,
                    this->heavy_reactant, this->heavy_product),
        ReactionDiagnosticType::density_src, rate_data_src, default_transformer);

    add_diagnostic(
        this->heavy_product, fmt::format("F{:s}", default_diag_suffix),
        fmt::format("Momentum transfer due to {:s} of {:s} to {:s}", long_reaction_type,
                    this->heavy_reactant, this->heavy_product),
        ReactionDiagnosticType::momentum_src, rate_data_src, default_transformer);

    add_diagnostic(
        this->heavy_product, fmt::format("E{:s}", default_diag_suffix),
        fmt::format("Energy transfer due to {:s} of {:s} to {:s}", long_reaction_type,
                    this->heavy_reactant, this->heavy_product),
        ReactionDiagnosticType::energy_src, rate_data_src, default_transformer);

    add_diagnostic("e", fmt::format("R{:s}", rad_diag_suffix),
                   fmt::format("Radiation loss due to {:s} of {:s} to {:s}",
                               long_reaction_type, this->heavy_reactant,
                               this->heavy_product),
                   ReactionDiagnosticType::energy_loss, rate_data_src, rad_transformer);
  }
}

///
void IznRecReaction::transform_additional(GuardedOptions& state,
                                          const RateData& rate_data) {

  // Extract the rate
  Field3D rate = rate_data.rate;

  // Extract heavy reactant properties
  GuardedOptions hr = state["species"][this->heavy_reactant];
  const BoutReal AA_rh = get<BoutReal>(hr["AA"]);
  const Field3D& n_rh = hr["density"].GetRef<Field3D>();
  const Field3D& v_rh = hr["velocity"].GetRef<Field3D>();

  // Extract heavy product properties
  GuardedOptions hp = state["species"][this->heavy_product];
  const Field3D& v_hp = hp["velocity"].GetRef<Field3D>();

  // Kinetic energy transfer to thermal energy
  //
  // This is a combination of three terms in the pressure
  // (thermal energy) equation:
  //
  // d/dt(3/2 p_1) = ... - F_12 v_1 + W_12 - (1/2) m R v_1^2 // From heavy reactant)
  //
  // d/dt(3/2 p_2) = ... - F_21 v_2 + W_21 + (1/2) m R v_2^2 // To heavy product)
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
  add(hp["energy_source"], 0.5 * AA_rh * rate * SQ(v_rh - v_hp));

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
      // Electrons with some velocity are incorporated into a species with their
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
      [&](BoutReal n_hr, BoutReal n_e, BoutReal T_e) {
        return n_hr * n_e
               * this->e_energy_loss_data->eval_sigma_vE_nT(T_e * Tnorm, n_e * Nnorm)
               * Nnorm / (Tnorm * FreqNorm) * radiation_multiplier;
      },
      region_no_bndry)(n_rh, n_e, T_e);

  // Loss is reduced by heating
  BoutReal electron_heating = this->e_energy_loss_data->get_metadata("electron_heating");
  energy_loss -= (electron_heating / Tnorm) * rate * radiation_multiplier;

  update_source<subtract<Field3D>>(state, "e", ReactionDiagnosticType::energy_loss,
                                   energy_loss);

  // Set collision frequency [s^-1] for heavy species
  std::string heavy_collfreq_lbl =
      fmt::format("{:s}_{:s}_{:s}", this->heavy_reactant, this->heavy_product,
                  this->short_reaction_type);
  set(state["species"][this->heavy_collfreq_species]["collision_frequencies"]
           [heavy_collfreq_lbl],
      rate_data.coll_freq(hr.name()));

  // N.B. nothing is done with the electron collision frequency at the moment.
  // Add it to the state too?
  Field3D electron_frequency = rate_data.coll_freq("e");
}

// ================================ IznReaction ===============================
///
IznReaction::IznReaction(std::string name, Options& options)
    : IznRecReaction("iz", name, options) {
  // Multiplier options are set on the heavy reactant species for ionisation
  this->rate_multiplier = options[this->heavy_reactant]["K_iz_multiplier"]
                              .doc("Scale the ionisation rate by this factor")
                              .withDefault<BoutReal>(1.0);

  this->radiation_multiplier = options[this->heavy_reactant]["R_ex_multiplier"]
                                   .doc("Scale the ionisation excitation/de-excitation "
                                        "radiation rate by this factor")
                                   .withDefault<BoutReal>(1.0);
}

// ================================ RecReaction ===============================
///
RecReaction::RecReaction(std::string name, Options& options)
    : IznRecReaction("rec", name, options) {
  // Multiplier options are set on the heavy product species for recombination
  this->rate_multiplier = options[this->heavy_product]["K_rec_multiplier"]
                              .doc("Scale the recombination rate by this factor")
                              .withDefault<BoutReal>(1.0);

  this->radiation_multiplier =
      options[this->heavy_product]["R_rec_multiplier"]
          .doc("Scale the recombination radiation (incl. 3 body) rate by this factor")
          .withDefault<BoutReal>(1.0);
}

} // namespace hermes