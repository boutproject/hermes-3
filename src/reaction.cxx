#include "reaction.hxx"

#include <iomanip>
#include <memory>
#include <regex>
#include <utility>

#include <bout/boutexception.hxx>

#include "integrate.hxx"

Reaction::Reaction(std::string name, Options& options) : name(name) {

  // Extract some relevant options, units to member vars for readability
  const auto& units = options["units"];
  Tnorm = get<BoutReal>(units["eV"]);
  Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
  FreqNorm = 1. / get<BoutReal>(units["seconds"]);

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
  int inst_num = get_instance_num() + 1;
  for (auto ii = 0; ii < inst_num; ii++) {
    std::getline(ss, reaction_str, ',');
  }

  // Parse the reaction string
  this->parser = std::make_unique<ReactionParser>(reaction_str);

  // Participation factors. All set to unity for now; could make configurable in future.
  for (const std::string& sp : this->parser->get_species()) {
    this->pfactors[sp] = 1;
  }

  // Initialise weight sums with dummy values. Real values are set on first call to
  // transform().
  this->momentum_weightsum = -1;
  this->energy_weightsum = -1;
}

/**
 * @brief
 *
 * @param sp_name
 * @param diag_name
 * @param long_diag_name
 * @param type
 * @param data_source
 *
 * @todo just use ReactionDiagnosticType in key directly?
 */
void Reaction::add_diagnostic(const std::string& sp_name, const std::string& diag_name,
                              const std::string& long_diag_name,
                              ReactionDiagnosticType type, const std::string& data_source,
                              DiagnosticTransformerType transformer,
                              const std::string& standard_name) {
  std::pair<std::string, ReactionDiagnosticType> diag_key = std::make_pair(sp_name, type);
  if (standard_name.empty()) {
    this->diagnostics.insert(
        std::make_pair(diag_key, ReactionDiagnostic(diag_name, long_diag_name, type,
                                                    data_source, transformer)));
  } else {
    this->diagnostics.insert(std::make_pair(
        diag_key, ReactionDiagnostic(diag_name, long_diag_name, type, data_source,
                                     standard_name, transformer)));
  }
}

/**
 * @brief Compute weight sums, if it hasn't been done already.
 *          Energy   : sum of (+ve pop change) participation factors
 *          Momentum : sum of (+ve pop change) participation factors, weighted by mass
 *
 * @param state current simulation state
 */
void Reaction::calc_weightsums(Options& state) {
  if (this->energy_weightsum < 0 || this->momentum_weightsum < 0) {
    this->momentum_weightsum = 0;
    this->energy_weightsum = 0;
    std::map<std::string, int> pop_changes = this->parser->get_stoich();
    for (const std::string& sp :
         this->parser->get_species(species_filter::heavy, species_filter::produced)) {
      int num_produced = pop_changes.at(sp);
      BoutReal pfac = pfactors.at(sp);
      this->momentum_weightsum +=
          num_produced * pfac * get<BoutReal>(state["species"][sp]["AA"]);
      this->energy_weightsum += num_produced * pfac;
    }
  }
}

/**
 * @brief Copy all diagnostics into the output
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

/**
 * @brief Use the stoichiometry matrix to compute density, momentum and energy sources.
 *
 * @param state
 */
void Reaction::transform(Options& state) {

  Field3D momentum_exchange, energy_exchange, energy_loss;
  zero_diagnostics(state);

  std::vector<std::string> reactant_names =
      parser->get_species(species_filter::reactants);

  // Extract electron properties
  Options& electron = state["species"]["e"];
  Field3D n_e = get<Field3D>(electron["density"]);
  Field3D T_e = get<Field3D>(electron["temperature"]);

  // Function passed to RateHelper to calculate reaction rate. Optionally scales by
  // multiplier.
  RateFunctionType calc_rate = [&](BoutReal mass_action, BoutReal ne, BoutReal te) {
    BoutReal result = mass_action * eval_reaction_rate(te * Tnorm, ne * Nnorm) * Nnorm
                      / FreqNorm * rate_multiplier;
    return result;
  };

  RateHelper rate_helper =
      RateHelper(state, reactant_names, calc_rate, n_e.getRegion("RGN_NOBNDRY"));
  Field3D reaction_rate = rate_helper.calc_rate();

  // Subclasses perform any additional transform tasks
  transform_additional(state, reaction_rate);

  // Use the stoichiometric values to set density sources for all species
  auto pop_changes = parser->get_stoich();
  for (const auto& [sp_name, pop_change] : pop_changes) {
    if (pop_change != 0) {
      // Density sources
      Field3D density_source = pfactors.at(sp_name) * pop_change * reaction_rate;
      update_source<add<Field3D>>(state, sp_name, ReactionDiagnosticType::density_src,
                                  density_source);
    }
  }

  // Get the species name(s) of heavy reactant, products
  std::vector<std::string> heavy_reactant_species =
      parser->get_species(reactant_names, species_filter::heavy);
  std::vector<std::string> heavy_product_species =
      parser->get_species(species_filter::heavy, species_filter::products);

  // Momentum and energy sources
  calc_weightsums(state);
  momentum_exchange = 0.0;
  energy_exchange = 0.0;
  for (const auto& [sp_name, pop_change_s] : pop_changes) {
    // No momentum, energy source for electrons due to pop change
    if (sp_name.compare("e") == 0) {
      continue;
    }
    Field3D momentum_source = 0.0;
    Field3D energy_source = 0.0;
    if (pop_change_s < 0) {
      // For species with net loss, sources follows directly from pop change
      momentum_source = pop_change_s * reaction_rate
                        * get<BoutReal>(state["species"][sp_name]["AA"])
                        * get<Field3D>(state["species"][sp_name]["velocity"]);
      energy_source = pop_change_s * reaction_rate * (3. / 2)
                      * get<Field3D>(state["species"][sp_name]["temperature"]);
    } else if (pop_change_s > 0) {
      // Species with net gain receive a proportion of the momentum and energy lost by
      // consumed reactants
      momentum_exchange = energy_exchange = 0;
      // Splitting factors - fraction of the total momentum/energy lost by consumed
      // species that will go to this product
      BoutReal momentum_split = pop_change_s * this->pfactors.at(sp_name)
                                * get<BoutReal>(state["species"][sp_name]["AA"])
                                / this->momentum_weightsum;
      BoutReal energy_split =
          pop_change_s * this->pfactors.at(sp_name) / this->energy_weightsum;
      for (auto& rsp_name : heavy_reactant_species) {
        // All consumed (net loss) reactants contribute
        int pop_change_r = pop_changes.at(rsp_name);
        if (pop_change_r < 0) {
          BoutReal momentum_split = this->pfactors.at(sp_name)
                                    * get<BoutReal>(state["species"][sp_name]["AA"])
                                    / this->momentum_weightsum;
          momentum_source += -pop_change_r * pfactors.at(rsp_name) * momentum_split
                             * reaction_rate
                             * get<BoutReal>(state["species"][rsp_name]["AA"])
                             * get<Field3D>(state["species"][rsp_name]["velocity"]);
          BoutReal energy_split = this->pfactors.at(sp_name) / this->energy_weightsum;
          energy_source += -pop_change_r * pfactors.at(rsp_name) * energy_split
                           * reaction_rate * (3. / 2)
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
 * @brief Reset the values of diagnostics stored in the state.
 *
 * @param state
 */
void Reaction::zero_diagnostics(Options& state) {
  if (this->diagnose) {
    for (auto& [key, diag] : diagnostics) {
      set<Field3D>(state[diag.name], 0.0);
    }
  }
}