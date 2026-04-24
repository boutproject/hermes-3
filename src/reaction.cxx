#include "reaction.hxx"

#include <iomanip>
#include <memory>
#include <numeric>
#include <regex>
#include <utility>

#include <bout/boutexception.hxx>

#include "integrate.hxx"
#include "reaction_settings.hxx"

namespace hermes {

///
Reaction::Reaction(std::string name, Options& options)
    : ReactionBase({readOnly("species:{sp}:{r_val}"), readOnly("species:e:{e_val}"),
                    readWrite("species:{sp}:{w_val}")}),
      units(options["units"]), name(name) {

  // Extract some relevant options, units to member vars for readability
  this->Tnorm = get<BoutReal>(this->units["eV"]);
  this->Nnorm = get<BoutReal>(this->units["inv_meters_cubed"]);
  this->FreqNorm = 1. / get<BoutReal>(this->units["seconds"]);

  this->diagnose = options[name]["diagnose"]
                       .doc("Output additional diagnostics?")
                       .withDefault<bool>(false);

  std::string reaction_str, data_src_id;
  ReactionDataTypes data_src_type;

  // Extract reaction string, data source type for this reaction
  get_reaction_settings(options[name], reaction_str, data_src_type, data_src_id);

  // Parse the reaction string
  this->parser = std::make_unique<ReactionParser>(reaction_str);

  std::vector<std::string> metadata_keys = {};
  this->rate_data = ReactionDataFactory::getInstance().create(
      toString(data_src_type), data_src_id, options, metadata_keys);

  std::vector<std::string> species = this->parser->get_species();
  // Participation factors. All set to unity for now; could make
  // configurable in future.
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

///
void Reaction::add_diagnostic(const std::string& sp_name, const std::string& diag_name,
                              const std::string& diag_desc,
                              ReactionDiagnosticType diag_type,
                              const std::string& data_source,
                              DiagnosticTransformerType transformer,
                              const std::string& standard_name) {
  std::pair<std::string, ReactionDiagnosticType> diag_key =
      std::make_pair(sp_name, diag_type);
  if (standard_name.empty()) {
    this->diagnostics.insert(std::make_pair(
        diag_key, ReactionDiagnostic(diag_name, diag_desc, diag_type, data_source,
                                     this->units, transformer)));
  } else {
    this->diagnostics.insert(std::make_pair(
        diag_key, ReactionDiagnostic(diag_name, diag_desc, diag_type, data_source,
                                     standard_name, this->units, transformer)));
  }
  setPermissions(readWrite(diag_name));
}

///
void Reaction::get_reaction_settings(Options& options, std::string& reaction_str,
                                     ReactionDataTypes& data_type, std::string& data_id) {

  // Extract reaction string(s). At least one must have been specified.
  std::vector<std::string> reaction_strs = split_csv_str(
      options["type"].doc("Comma-separated list of reaction strings").as<std::string>());
  const std::size_t num_reactions = reaction_strs.size();
  ASSERT1(this->inst_num <= num_reactions);
  const std::size_t inst_idx = inst_num - 1;
  reaction_str = reaction_strs[inst_idx];

  // Extract data source type(s). Default to Amjuel for all if none was provided.
  std::vector<std::string> data_type_strs = split_csv_str(
      options["data_srcs"]
          .doc("Reaction data source type ('ADAS', 'Amjuel', etc.), either a "
               "single value for all reaction strings, "
               "or a comma-separated list with one entry per reaction string")
          .withDefault(toString(ReactionDataTypes::Amjuel)),
      num_reactions, "reaction data source type");
  std::vector<ReactionDataTypes> data_types;
  std::transform(data_type_strs.begin(), data_type_strs.end(), data_type_strs.begin(),
                 [](std::string s) {
                   std::transform(s.begin(), s.end(), s.begin(), ::tolower);
                   return s;
                 });
  std::transform(data_type_strs.begin(), data_type_strs.end(),
                 std::back_inserter(data_types), ReactionDataTypesFromString);
  data_type = data_types[inst_idx];

  // Extract data id(s). Set a sensible default if none was provided.
  std::string data_ids_str =
      options["data_ids"]
          .doc("Comma-separated list of data identifiers to use for each reaction "
               "string.")
          .withDefault(get_default_data_ids_str(reaction_strs, data_types,
                                                ReactionCoeffTypes::sigma_v));
  std::vector<std::string> data_ids =
      split_csv_str(data_ids_str, num_reactions, "reaction data source type");
  data_id = data_ids[inst_idx];
  if (data_id == NO_DATA_ID_DEFAULT_FOUND) {
    throw BoutException(fmt::format(
        "No reaction data id specified and no suitable default found for reaction "
        "string '{}'. Please provide a data id using the 'data_ids' option.",
        reaction_str));
  }
}

///
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

///
void Reaction::outputVars(Options& state) {
  if (this->diagnose) {
    for (auto& [key, diag] : this->diagnostics) {
      diag.add_to_state(state);
    }
  }
}

///
void Reaction::set_energy_channel_weight(const std::string& reactant_name,
                                         const std::string& product_name,
                                         BoutReal weight) {
  this->energy_channels[reactant_name][product_name] = weight;
}

///
void Reaction::set_momentum_channel_weight(const std::string& reactant_name,
                                           const std::string& product_name,
                                           BoutReal weight) {
  this->momentum_channels[reactant_name][product_name] = weight;
}

///
void Reaction::transform_impl(GuardedOptions& state) {
  zero_diagnostics(state);

  // Get the species name(s) of reactants, heavy reactants and products
  std::vector<std::string> reactant_names =
      parser->get_species(species_filter::reactants);
  std::vector<std::string> heavy_reactant_species =
      parser->get_species(reactant_names, species_filter::heavy);
  std::vector<std::string> heavy_product_species =
      parser->get_species(species_filter::heavy, species_filter::products);

  // All reaction sources are computed in interior region only
  Region<Ind3D> rng_nobndry = get<Field3D>(state["species"][reactant_names[0]]["density"])
                                  .getRegion("RGN_NOBNDRY");

  // Create rate helper and compute reaction rate, collision frequencies
  RateData rate_calc_results;
  RateParamsTypes rate_params_type = this->rate_data->get_fit_type();
  if (rate_params_type == RateParamsTypes::ET) {
    throw BoutException("RateParamsTypes::ET not implemented");
  } else if (rate_params_type == RateParamsTypes::nT) {
    TwoDRateFunc calc_rate = [&](BoutReal mass_action, BoutReal ne, BoutReal te) {
      BoutReal result = mass_action
                        * this->rate_data->eval_sigma_v_nT(te * Tnorm, ne * Nnorm) * Nnorm
                        / FreqNorm * rate_multiplier;
      return result;
    };
    auto rate_helper =
        RateHelper<RateParamsTypes::nT>(state, units, reactant_names, rng_nobndry);
    rate_calc_results = rate_helper.calc_rates(calc_rate, this->do_parallel_averaging);
  } else if (rate_params_type == RateParamsTypes::T) {
    OneDRateFunc calc_rate = [&](BoutReal mass_action, BoutReal Teff) {
      BoutReal result = mass_action * 1e-6 * this->rate_data->eval_sigma_v_T(Teff) * Nnorm
                        / FreqNorm * rate_multiplier;
      return result;
    };

    auto rate_helper =
        RateHelper<RateParamsTypes::T>(state, units, reactant_names, rng_nobndry);

    rate_calc_results = rate_helper.calc_rates(calc_rate, this->do_parallel_averaging);
  } else {
    throw BoutException("Unhandled RateParamsTypes in Reaction::transform_impl()");
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
  Field3D density_source(0.0);
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

///
void Reaction::zero_diagnostics(GuardedOptions& state) {
  if (this->diagnose) {
    for (auto& [key, diag] : diagnostics) {
      set<Field3D>(state[diag.get_name()], 0.0);
    }
  }
}

} // namespace hermes
