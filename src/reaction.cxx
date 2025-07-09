#include <iomanip>
#include <memory>
#include <regex>

#include "integrate.hxx"

#include "reaction.hxx"

Reaction::Reaction(std::string name, Options& options) : name(name) {

  // Extract some relevant options, units to member vars for readability
  const auto& units = options["units"];
  Tnorm = get<BoutReal>(units["eV"]);
  Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
  FreqNorm = 1. / get<BoutReal>(units["seconds"]);

  diagnose = options[name]["diagnose"]
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
 * @brief Use the stoichiometry matrix to compute density, momentum and energy sources.
 * ASSUMES EXACTLY 2 REACTANTS (electron(s) + 1 ion/neutral) AND 2 PRODUCTS (electron(s) +
 * 1 ion/neutral) for now.
 *
 * @param state
 */
void Reaction::transform(Options& state) {

  Field3D momentum_exchange, energy_exchange, energy_loss;

  std::vector<std::string> reactant_names =
      parser->get_species(species_filter::reactants);

  // Restrict to 2 reactants for now;
  ASSERT1(reactant_names.size() == 2);

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

  // Use the stoichiometric values to set density sources for all species
  auto pop_changes = parser->get_stoich();
  for (auto el : pop_changes) {
    std::string sp_name = el.first;
    int pop_change = el.second;
    if (pop_change != 0) {
      // Density sources
      add(state["species"][sp_name]["density_source"], pop_change * reaction_rate);
    }
  }

  // Get the species name(s) of heavy reactant, products
  std::vector<std::string> heavy_reactant_species =
      parser->get_species(reactant_names, species_filter::heavy);
  std::vector<std::string> heavy_product_species =
      parser->get_species(species_filter::heavy, species_filter::products);

  // Momentum and energy sources
  calc_weightsums(state);
  for (auto el : pop_changes) {
    std::string sp_name = el.first;
    // No momentum, energy source for electrons due to pop change
    if (sp_name.compare("e") == 0) {
      continue;
    }
    int pop_change = el.second;
    // Species momentum
    auto Gs = get<BoutReal>(state["species"][sp_name]["AA"])
              * get<Field3D>(state["species"][sp_name]["velocity"]);
    // Species energy
    auto Ws = (3. / 2) * get<Field3D>(state["species"][sp_name]["temperature"]);

    if (pop_change < 0) {
      // For species with net loss, sources follows directly from pop change
      momentum_exchange = pop_change * reaction_rate * Gs;
      energy_exchange = pop_change * reaction_rate * Ws;
    } else if (pop_change > 0) {
      // Species with net gain receive a proportion of the momentum and energy lost by
      // consumed reactants
      momentum_exchange = energy_exchange = 0;
      for (auto& rsp_name : heavy_reactant_species) {
        // All consumed (net loss) reactants contribute
        if (pop_changes[rsp_name] < 0) {
          auto Gr = get<BoutReal>(state["species"][rsp_name]["AA"])
                    * get<Field3D>(state["species"][rsp_name]["velocity"]);
          BoutReal momentum_split = this->pfactors.at(sp_name)
                                    * get<BoutReal>(state["species"][sp_name]["AA"])
                                    / this->momentum_weightsum;
          momentum_exchange +=
              pfactors.at(rsp_name) * momentum_split * reaction_rate * Gr;

          BoutReal energy_split = this->pfactors.at(sp_name) / this->energy_weightsum;
          auto Wr = (3. / 2) * get<Field3D>(state["species"][rsp_name]["temperature"]);
          energy_exchange += pfactors.at(rsp_name) * energy_split * reaction_rate * Wr;
        }
      }
    } else {
      // No pop change
      continue;
    }

    // Update sources
    add(state["species"][sp_name]["momentum_source"], momentum_exchange);
    add(state["species"][sp_name]["energy_source"], energy_exchange);
  }

  // Subclasses perform any additional transform tasks
  transform_additional(state, reaction_rate, momentum_exchange, energy_exchange,
                       energy_loss);

  set_diagnostic_fields(reaction_rate, momentum_exchange, energy_exchange, energy_loss);
}

/**
 * @brief Construct a new RateHelper.
 *
 * @tparam LimiterType
 * @tparam RegionType
 * @param state
 * @param reactant_names
 * @param rate_calc_func
 * @param region
 */
template <typename LimiterType, typename IdxType>
RateHelper<LimiterType, IdxType>::RateHelper(
    const Options& state, const std::vector<std::string>& reactant_names,
    RateFunctionType rate_calc_func, const Region<IdxType> region)
    : rate_calc_func(rate_calc_func), region(region) {

  // Extract electron properties from state
  const Options& electron = state["species"]["e"];
  n_e = get<Field3D>(electron["density"]);
  T_e = get<Field3D>(electron["temperature"]);

  // Extract and store reactant densities
  std::transform(reactant_names.begin(), reactant_names.end(),
                 std::back_inserter(n_reactants), [&](const std::string& reactant_name) {
                   return get<Field3D>(state["species"][reactant_name]["density"]);
                 });
}

/**
 * @brief Compute the cell-averaged reaction rate, accounting for the mass action factor
 * (product of reactant densities)
 *
 * @tparam LimiterType
 * @tparam IdxType
 * @return Field3D the cell-averaged reaction rate
 */
template <typename LimiterType, typename IdxType>
Field3D RateHelper<LimiterType, IdxType>::calc_rate() {
  Field3D reaction_rate{emptyFrom(n_e)};
  auto J = reaction_rate.getCoordinates()->J;
  BOUT_FOR(i, region) {

    auto yp = i.yp();
    auto ym = i.ym();
    auto Ji = J[i];

    reaction_rate[i] =
        4. / 6 * rate_calc_func(mass_action(i), n_e[i], T_e[i])
        + (Ji + J[ym]) / (12. * Ji)
              * rate_calc_func(mass_action_left(i, ym, yp),
                               cellLeft<LimiterType>(n_e[i], n_e[ym], n_e[yp]),
                               cellLeft<LimiterType>(T_e[i], T_e[ym], T_e[yp]))
        + (Ji + J[yp]) / (12. * Ji)
              * rate_calc_func(mass_action_right(i, ym, yp),
                               cellRight<LimiterType>(n_e[i], n_e[ym], n_e[yp]),
                               cellRight<LimiterType>(T_e[i], T_e[ym], T_e[yp]));
  }
  return reaction_rate;
}

template <typename LimiterType, typename IdxType>
BoutReal RateHelper<LimiterType, IdxType>::mass_action(IdxType i) {
  BoutReal result = 1;
  for (const auto& n : n_reactants) {
    result *= n[i];
  }
  return result;
}

template <typename LimiterType, typename IdxType>
BoutReal RateHelper<LimiterType, IdxType>::mass_action_left(IdxType i, IdxType ym,
                                                            IdxType yp) {
  BoutReal result = 1;
  for (const auto& n : n_reactants) {
    result *= cellLeft<LimiterType>(n[i], n[ym], n[yp]);
  }
  return result;
}

template <typename LimiterType, typename IdxType>
BoutReal RateHelper<LimiterType, IdxType>::mass_action_right(IdxType i, IdxType ym,
                                                             IdxType yp) {
  BoutReal result = 1;
  for (const auto& n : n_reactants) {
    result *= cellRight<LimiterType>(n[i], n[ym], n[yp]);
  }
  return result;
}