#include "reaction.hxx"

#include <iomanip>
#include <memory>
#include <numeric>
#include <regex>
#include <utility>

#include <bout/boutexception.hxx>

#include "integrate.hxx"

Reaction::Reaction(std::string name, Options& options)
    : ReactionBase({readOnly("species:{sp}:{r_val}"), readOnly("species:e:{e_val}"),
                    readWrite("species:{sp}:{w_val}")}),
      units(options["units"]), name(name) {

  // Extract some relevant options, units to member vars for readability
  Tnorm = get<BoutReal>(this->units["eV"]);
  Nnorm = get<BoutReal>(this->units["inv_meters_cubed"]);
  FreqNorm = 1. / get<BoutReal>(this->units["seconds"]);

  this->diagnose = options[name]["diagnose"]
                       .doc("Output additional diagnostics?")
                       .withDefault<bool>(false);

  /*
   * Awful hack to extract the correct reaction expression from the params; depends on
   * instantiation order matching the order reactions are listed in the input file. There
   * must be a better way...
   */
  std::string reaction_grp_str = options[name]["type"];
  std::regex match_parentheses("\\(|\\)");
  reaction_grp_str = std::regex_replace(reaction_grp_str, match_parentheses, "");
  std::string reaction_str;
  std::stringstream ss(reaction_grp_str);
  for (auto ii = 0; ii < this->inst_num; ii++) {
    std::getline(ss, reaction_str, ',');
  }

  // Parse the reaction string
  this->parser = std::make_unique<ReactionParser>(reaction_str);

  std::vector<std::string> species = this->parser->get_species();
  // Participation factors. All set to unity for now; could make configurable in future.
  for (const std::string& sp : species) {
    this->pfactors[sp] = 1;
  }

  // Initialise momentum/energy channel maps
  for (const std::string& reactant :
       this->parser->get_species(species_filter::heavy, species_filter::reactants)) {
    if (this->energy_channels.count(reactant) == 0) {
      energy_channels[reactant] = std::map<std::string, BoutReal>();
    }
    if (this->momentum_channels.count(reactant) == 0) {
      momentum_channels[reactant] = std::map<std::string, BoutReal>();
    }
  }

  substitutePermissions("sp", species);
  substitutePermissions("r_val", {"AA", "density", "velocity", "temperature"});
  substitutePermissions("e_val", {"density", "temperature"});
  substitutePermissions("w_val", {"momentum_source", "energy_source", "density_source"});
  setPermissions(readWrite("species:{reactant}:collision_frequency"));
  substitutePermissions("reactant", this->parser->get_species(species_filter::reactants));
}

/**
 * @brief Add a new diagnostic.
 *
 * @param sp_name Species with which the diagnostic will be associated
 * @param diag_name Label used in the output (and to store it temporarily in the state)
 * @param diag_desc Description to use as the 'long_name' output attribute
 * @param diag_type enum identifying the diagnostic type, also used to determine source
 * name
 * @param data_source Name to use as the 'source' output attribute
 * @param standard_name Optional string to use as the 'standard_name' output attribute.
 * Defaults to diag_name.
 * @param transformer Optional transformer function to use when modifying the diagnostic
 * (default is 'negate', i.e. the diagnostic has the opposite sign to the source)
 */
void Reaction::add_diagnostic(const std::string& sp_name, const std::string& diag_name,
                              const std::string& diag_desc,
                              ReactionDiagnosticType diag_type,
                              const std::string& data_source,
                              DiagnosticTransformerType transformer,
                              const std::string& standard_name) {
  std::pair<std::string, ReactionDiagnosticType> diag_key =
      std::make_pair(sp_name, diag_type);
  if (standard_name.empty()) {
    this->diagnostics.insert(
        std::make_pair(diag_key, ReactionDiagnostic(diag_name, diag_desc, diag_type,
                                                    data_source, transformer)));
  } else {
    this->diagnostics.insert(std::make_pair(
        diag_key, ReactionDiagnostic(diag_name, diag_desc, diag_type, data_source,
                                     standard_name, transformer)));
  }
  setPermissions(readWrite(diag_name));
}

/**
 * @brief Set weights for any reactant => product momentum / energy channel that hasn't
 * already been specified via set_energy_channel_weight and set_momentum_channel_weight.
 *
 * @param state current simulation state
 */
void Reaction::init_channel_weights(GuardedOptions& state) {
  std::vector<std::string> heavy_reactants =
      this->parser->get_species(species_filter::heavy, species_filter::reactants);
  std::vector<std::string> heavy_products =
      this->parser->get_species(species_filter::heavy, species_filter::products);

  // If all channels already have values, bail out
  std::size_t num_energy_channels_set = 0, num_momentum_channels_set = 0;
  const std::size_t num_channels_expected =
      heavy_reactants.size() * heavy_products.size();
  for (const auto& reactant : heavy_reactants) {
    num_energy_channels_set += this->energy_channels[reactant].size();
    num_momentum_channels_set += this->momentum_channels[reactant].size();
  }
  if (num_energy_channels_set == num_channels_expected
      && num_momentum_channels_set == num_channels_expected) {
    return;
  }

  // Compute total weights:
  BoutReal momentum_weightsum = 0, energy_weightsum = 0;
  for (const std::string& sp :
       this->parser->get_species(species_filter::heavy, species_filter::produced)) {
    int num_produced = this->parser->pop_change(sp);
    BoutReal pfac = pfactors.at(sp);
    momentum_weightsum += num_produced * pfac * get<BoutReal>(state["species"][sp]["AA"]);
    energy_weightsum += num_produced * pfac;
  }
  // Set default values for any unset channels
  for (const std::string& reactant : heavy_reactants) {
    for (const std::string& product : heavy_products) {
      if (this->energy_channels[reactant].count(product) == 0) {
        this->energy_channels[reactant][product] = this->parser->pop_change(product)
                                                   * this->pfactors.at(product)
                                                   / energy_weightsum;
      }
      if (this->momentum_channels[reactant].count(product) == 0) {
        this->momentum_channels[reactant][product] =
            this->parser->pop_change(product) * this->pfactors.at(product)
            * get<BoutReal>(state["species"][product]["AA"]) / momentum_weightsum;
      }
    }
  }

  /* Make sure we're not trying to distribute < 0 or > 1 times total momentum/energy of
   * each reactant. The total weights could be restricted to be exactly 1, but we want to
   * allow momentum/energy contributions from certain species to be turned off.
   */
  for (const std::string& reactant : heavy_reactants) {
    double total_energy_weight = std::accumulate(
        this->energy_channels[reactant].begin(), this->energy_channels[reactant].end(),
        0.0, [](double sum, const auto& pair) { return sum + pair.second; });
    ASSERT0(total_energy_weight >= 0 && total_energy_weight <= 1);
    double total_momentum_weight =
        std::accumulate(this->momentum_channels[reactant].begin(),
                        this->momentum_channels[reactant].end(), 0.0,
                        [](double sum, const auto& pair) { return sum + pair.second; });
    ASSERT0(total_momentum_weight >= 0 && total_momentum_weight <= 1);
  }
}

/**
 * @brief Copy all diagnostics into the output, setting the appropriate metadata at the
 * same time
 *
 * @param state
 */
void Reaction::outputVars(Options& state) {
  if (this->diagnose) {
    for (auto& [key, diag] : this->diagnostics) {
      diag.add_to_state(state);
    }
  }
}

void Reaction::set_energy_channel_weight(const std::string& reactant_name,
                                         const std::string& product_name,
                                         BoutReal weight) {
  this->energy_channels[reactant_name][product_name] = weight;
}

void Reaction::set_momentum_channel_weight(const std::string& reactant_name,
                                           const std::string& product_name,
                                           BoutReal weight) {
  this->momentum_channels[reactant_name][product_name] = weight;
}

/**
 * @brief Add density, momentum and energy sources that apply to all reactions (e.g.
 * those driven by species population changes), then call transform_additional() to allow
 * subclasses to add other terms.
 *
 * @param state
 */
void Reaction::transform_impl(GuardedOptions& state) {

  Field3D momentum_exchange, energy_exchange, energy_loss;
  zero_diagnostics(state);

  // Get the species name(s) of reactants, heavy reactants and products
  std::vector<std::string> reactant_names =
      parser->get_species(species_filter::reactants);
  std::vector<std::string> heavy_reactant_species =
      parser->get_species(reactant_names, species_filter::heavy);
  std::vector<std::string> heavy_product_species =
      parser->get_species(species_filter::heavy, species_filter::products);

  // First reactant; just used to get region
  Field3D first_reactant = get<Field3D>(state["species"][reactant_names[0]]["density"]);

  // Create rate helper and compute reaction rate, collision frequencies
  RateData rate_calc_results;
  RateParamsTypes rate_params_type = get_rate_params_type();
  if (rate_params_type == RateParamsTypes::ET) {
    throw BoutException("RateParamsTypes::ET not implemented");
  } else if (rate_params_type == RateParamsTypes::nT) {
    TwoDRateFunc calc_rate = [&](BoutReal mass_action, BoutReal ne, BoutReal te) {
      BoutReal result = mass_action * eval_sigma_v_nT(te * Tnorm, ne * Nnorm) * Nnorm
                        / FreqNorm * rate_multiplier;
      return result;
    };
    auto rate_helper = RateHelper<RateParamsTypes::nT>(
        state, units, reactant_names, first_reactant.getRegion("RGN_NOBNDRY"));
    rate_calc_results = rate_helper.calc_rates(calc_rate, this->do_parallel_averaging);
  } else if (rate_params_type == RateParamsTypes::T) {
    OneDRateFunc calc_rate = [&](BoutReal mass_action, BoutReal Teff) {
      BoutReal result =
          mass_action * 1e-6 * eval_sigma_v_T(Teff) * Nnorm / FreqNorm * rate_multiplier;
      return result;
    };

    auto rate_helper = RateHelper<RateParamsTypes::T>(
        state, units, reactant_names, first_reactant.getRegion("RGN_NOBNDRY"));

    rate_calc_results = rate_helper.calc_rates(calc_rate, this->do_parallel_averaging);
  } else {
    throw BoutException("Unhandled RateParamsTypes in Reaction::transform()");
  }

  // Set collision frequencies
  for (const auto& reactant_name : reactant_names) {
    update_source<set<Field3D>>(state, reactant_name,
                                ReactionDiagnosticType::collision_freq,
                                rate_calc_results.coll_freq(reactant_name));
  }

  // Subclasses perform any additional transform tasks
  transform_additional(state, rate_calc_results);

  // Use the stoichiometric values to set density sources for all species
  Field3D density_source;
  for (const auto& sp_name : this->parser->get_species()) {
    int pop_change = this->parser->pop_change(sp_name);
    if (pop_change != 0) {
      // Density sources
      density_source = pfactors.at(sp_name) * pop_change * rate_calc_results.rate;
      update_source<add<Field3D>>(state, sp_name, ReactionDiagnosticType::density_src,
                                  density_source);
    }
  }

  // Population change-driven sources for all species other than electrons
  init_channel_weights(state);
  momentum_exchange = 0.0;
  energy_exchange = 0.0;
  Field3D momentum_source, energy_source;
  for (const auto& [sp_name, pop_change_s] : this->parser->get_mom_energy_pop_changes()) {
    // No momentum, energy source for electrons due to pop change
    if (sp_name.compare("e") == 0) {
      continue;
    }
    momentum_source = 0.0;
    energy_source = 0.0;
    if (pop_change_s < 0) {
      // For species with net loss, sources follows directly from pop change
      momentum_source = pop_change_s * rate_calc_results.rate
                        * get<BoutReal>(state["species"][sp_name]["AA"])
                        * get<Field3D>(state["species"][sp_name]["velocity"]);
      energy_source = pop_change_s * rate_calc_results.rate * (3. / 2)
                      * get<Field3D>(state["species"][sp_name]["temperature"]);
    } else if (pop_change_s > 0) {
      // Species with net gain receive a proportion of the momentum and energy lost by
      // consumed reactants. See init_channel_weights() for default splitting factors.
      momentum_exchange = energy_exchange = 0;
      for (auto& rsp_name : heavy_reactant_species) {
        // All consumed (net loss) reactants can contribute
        int pop_change_r = this->parser->pop_change_reactant(rsp_name);
        if (pop_change_r < 0) {
          momentum_source += -pop_change_r * pfactors.at(rsp_name)
                             * this->momentum_channels[rsp_name][sp_name]
                             * rate_calc_results.rate
                             * get<BoutReal>(state["species"][rsp_name]["AA"])
                             * get<Field3D>(state["species"][rsp_name]["velocity"]);
          energy_source += -pop_change_r * pfactors.at(rsp_name)
                           * this->energy_channels[rsp_name][sp_name]
                           * rate_calc_results.rate * (3. / 2)
                           * get<Field3D>(state["species"][rsp_name]["temperature"]);
        }
      }
      momentum_exchange += momentum_source;
      energy_exchange += energy_source;
    } else {
      // No pop change
      continue;
    }

    // Update sources
    update_source<add<Field3D>>(state, sp_name, ReactionDiagnosticType::momentum_src,
                                momentum_source);
    update_source<add<Field3D>>(state, sp_name, ReactionDiagnosticType::energy_src,
                                energy_source);
  }
}

/**
 * @brief Reset the temporary values of the diagnostics stored in the state.
 *
 * @param state
 */
void Reaction::zero_diagnostics(GuardedOptions& state) {
  if (this->diagnose) {
    for (auto& [key, diag] : diagnostics) {
      set<Field3D>(state[diag.get_name()], 0.0);
    }
  }
}
