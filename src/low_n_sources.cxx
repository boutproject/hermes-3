#include "low_n_sources.hxx"

#include <algorithm>
#include <string>

#include <bout/boutexception.hxx>
#include <bout/options.hxx>
#include <bout/output.hxx>
#include <bout/utils.hxx>

#include "guarded_options.hxx"
#include "species_parser.hxx"

// ======================= Helper struct and functions ========================
namespace {

/// The information needed to define a pseudo-reaction
using PseudoReactionSpec = hermes::PseudoReactionSpec;

/**
 * @brief Smoothly saturate `x` between `lo` and `hi`, avoiding the discontinuous
 * derivative of std::clamp() at the two bounds (helps Newton-based solvers).
 *
 * @param x The value to clamp
 * @param lo Lower bound, approached smoothly as x -> -inf
 * @param hi Upper bound, approached smoothly as x -> +inf
 * @param width Transition width between lo and hi; smaller is closer to a hard clamp
 */
BoutReal smooth_clamp(BoutReal x, BoutReal lo, BoutReal hi, BoutReal width) {
  return lo + (hi - lo) / (1.0 + exp(-(x - 0.5 * (lo + hi)) / width));
}

/**
 * @brief Smoothly transition from 0 (ratio << 1) to 1 (ratio >> 1), replacing a hard
 * `ratio < 1 ? 0 : 1` step (helps Newton-based solvers) with a differentiable one.
 *
 * @param ratio A value expected to vary around 1 (e.g. density / density_threshold)
 * @param width Transition width; smaller is closer to a hard step
 */
BoutReal smooth_gate(BoutReal ratio, BoutReal width) {
  return 1.0 / (1.0 + exp(-(ratio - 1.0) / width));
}

/**
 * @brief Append a (pseudo-)reaction definition if both species are present
 *
 * @param first_species The first species
 * @param second_species The second species
 * @param type The type of pseudo-reaction
 * @param reaction_defs The vector of reaction definitions to append to
 */
void append_reaction_def_if_present(
    const std::string& first_species, const std::string& second_species,
    PseudoReactionType type,
    std::vector<hermes::LowNSources::ReactionDef>& reaction_defs) {
  if (first_species.empty() or second_species.empty()) {
    return;
  }
  reaction_defs.emplace_back(first_species, second_species, type);
}

/**
 * @brief Assemble the information needed to define a pseudo-reaction: its reaction
 * string and target species (i.e. the species whose density the reaction restores).
 *
 * @param type The type of pseudo-reaction
 * @param species1 The first species involved in the reaction
 * @param species2 A possible second species involved in the reaction
 * @return PseudoReactionSpec The reaction specification
 */
PseudoReactionSpec make_reaction_spec(PseudoReactionType type,
                                      const std::string& species1,
                                      const std::string& species2) {
  PseudoReactionSpec spec;
  spec.type = type;

  switch (type) {
  case PseudoReactionType::association: {
    // Allow for molecules made up of any number of atoms (of the same element)
    const SpeciesParser molecule(species2);
    const int atoms_per_molecule =
        stringToInt(molecule.get_base_species().substr(molecule.get_element().size()));
    // N.B. SpeciesParser constructor reject molecule strings with <1 atoms, so we can assume
    // atoms_per_molecule is +ve
    spec.reaction_str = fmt::format("{}{} -> {}", atoms_per_molecule, species1, species2);
    spec.target_species = species2;
    return spec;
  }
  case PseudoReactionType::recombination:
    spec.reaction_str = fmt::format("{} + e -> {}", species1, species2);
    spec.target_species = species2;
    return spec;
  case PseudoReactionType::ionization:
    spec.reaction_str = fmt::format("{} -> {} + e", species1, species2);
    spec.target_species = species2;
    return spec;
  case PseudoReactionType::external_src:
    // N.B. Since no species in this reaction has a -ve pop. change, Reaction::transform doesn't add momentum and energy sources.
    // These are added in PseudoReaction::transform_impl instead, unless turned off with external_src_mom_energy_enabled=false.
    spec.reaction_str = fmt::format("{} -> 2{}", species1, species1);
    spec.target_species = species1;
    return spec;
  }
  throw BoutException("Unhandled PseudoReactionType in make_reaction_spec()");
}

/// @brief Convenience function to set a unique species or throw an exception if one is already set
void set_unique_species_or_throw(std::string& existing_species,
                                 const std::string& new_species,
                                 const std::string& base_species,
                                 const std::string& species_category) {
  if (!existing_species.empty()) {
    throw BoutException(
        fmt::format("LowNSources supports at most one {} species per base species. "
                    "Base species '{}' has '{}' and '{}'.",
                    species_category, base_species, existing_species, new_species));
  }
  existing_species = new_species;
}

} // namespace

namespace hermes {

// ======================= LowNSources member functions =======================
LowNSources::LowNSources(const std::string& name, Options& options,
                         [[maybe_unused]] Solver* solver)
    : LowNSources(name, options) {}

LowNSources::LowNSources(const std::string& name, Options& options)
    : NamedComponent(name, {readOnly("species:{all_species}:{read_vals}"),
                            readWrite("species:{all_species}:{write_vals}")}),
      options(options) {
  substitutePermissions("read_vals", {"AA", "density", "velocity", "temperature"});
  substitutePermissions("write_vals", {"density_source", "momentum_source",
                                       "energy_source", "collision_frequency"});

  // Read the list of components
  const std::string component_names = options["hermes"]["components"]
                                          .doc("Components in order of execution")
                                          .withDefault<std::string>("");

  // Read control parameters
  /// Floor factor
  this->floor_factor = options[name]["floor_factor"]
                           .doc("Multiplier applied to the low-density floor")
                           .withDefault<BoutReal>(10.0);
  /// Timescale over which low density sources act (input value is in normalised units;
  /// default is 1e-5 s converted to normalised units)
  const BoutReal seconds_norm = get<BoutReal>(options["units"]["seconds"]);
  this->timescale = options[name]["timescale"]
                        .doc("Timescale over which low density sources act "
                             "[normalised]")
                        .withDefault<BoutReal>(1e-5 / seconds_norm);
  /// Rate calculation strategy
  this->strategy =
      options[name]["strategy"]
          .doc("Strategy to use when generating sources to avoid low species density.")
          .withDefault(toString(LowNSrcsStrategy::standard));
  // Diagnostics output switch
  this->diagnose = options[name]["diagnose"]
                       .doc("Whether or not to output diagnostics for the low-n sources")
                       .withDefault<bool>(false);

  /// Whether to enable momentum and energy terms for external ion sources
  this->external_src_mom_energy_enabled =
      options[name]["external_src_mom_energy_enabled"]
          .doc("Whether to enable momentum and energy sources for external ion sources")
          .withDefault<bool>(true);

  // Classify species by element and by type
  std::vector<std::string> species_with_mass;
  for (const auto& component : strsplit(component_names, ',')) {
    const std::string component_name = trim(component, " \t\r()");
    // Species components must have 'AA' set
    if (component_name.empty() or !options[component_name].isSet("AA")) {
      continue;
    }
    species_with_mass.push_back(component_name);
    SpeciesParser species(component_name);
    append_species_by_type(species, species.get_type());
    this->all_species.push_back(component_name);
  }

  // LowNSources can be instantiated directly in unit tests, bypassing scheduler-level
  // declareAllSpecies() substitutions, so resolve placeholders here.
  substitutePermissions("all_species", species_with_mass);

  // Set activation thresholds for each species now that species lists are populated
  for (const auto& species_name : this->all_species) {
    const BoutReal density_floor = options[species_name]["density_floor"]
                                       .doc("Minimum density floor")
                                       .withDefault<BoutReal>(1e-7);
    this->density_thresholds[species_name] = this->floor_factor * density_floor;
  }

  instantiate_reactions();
}

void LowNSources::append_species_by_type(const SpeciesParser& species, SpeciesType type) {
  if (type != SpeciesType::electron
      and (species.get_charge() < 0 or species.get_charge() > 1)) {
    output_warn.write(
        "LowNSources: species '{}' has unsupported ion charge {}; no low-n sources "
        "will be added for this species.\n",
        species.get_str(), species.get_charge());
    return;
  }

  // The full species string
  const std::string& species_str = species.get_str();
  const std::string element_str = species.get_element();

  // First, register the species as an ion/neutral/molecule of a particular element
  auto [element_group_it, element_group_inserted] =
      this->species_by_base.try_emplace(element_str);
  if (element_group_inserted) {
    this->base_species_order.push_back(element_str);
  }
  auto& element_group = element_group_it->second;

  if (type == SpeciesType::ion) {
    set_unique_species_or_throw(element_group.ion, species_str, element_str, "ion");
  } else if (type == SpeciesType::atom) {
    set_unique_species_or_throw(element_group.neutral, species_str, element_str,
                                "neutral");
  } else if (type == SpeciesType::molecule) {
    set_unique_species_or_throw(element_group.molecule, species_str, element_str,
                                "molecule");
  }

  // If the species is molecular add a separate group
  //  (Since molecules and molecular ions may need to exchange particles)
  if (species.is_molecule()) {
    // Get the neutral molecule base (which is the same as species_str for neutral molecules, but not for molecular ions)
    const std::string& mol_str = species.get_base_species();
    auto [mol_group_it, mol_group_inserted] = this->species_by_base.try_emplace(mol_str);
    if (mol_group_inserted) {
      this->base_species_order.push_back(mol_str);
    }
    auto& mol_group = this->species_by_base.at(mol_str);
    if (type == SpeciesType::molecular_ion) {
      // Treat it as an 'ion' of the molecule base
      set_unique_species_or_throw(mol_group.ion, species_str, mol_str, "ion");
    } else if (type == SpeciesType::molecule) {
      // Treat it as a 'neutral' of the molecule base
      set_unique_species_or_throw(mol_group.neutral, species_str, mol_str, "neutral");
    }
  }
  this->reaction_defs_populated = false;
}

/**
 * @brief Getter for the reaction definitions; used for testing.
 *
 * @return std::vector<LowNSources::ReactionDef>
 */
std::vector<LowNSources::ReactionDef> LowNSources::get_reaction_defs() const {
  populate_reaction_defs();
  return this->reaction_defs;
}

void LowNSources::instantiate_reactions() {
  this->reaction_options_store.clear();
  this->reactions.clear();

  std::size_t num_reactions = 0;
  for (const auto& element : this->base_species_order) {
    const auto& groups = this->species_by_base.at(element);
    if (!groups.ion.empty() and !groups.neutral.empty()) {
      // recombination, ionization, external_src of ions
      num_reactions += 3;
    } else if (!groups.ion.empty()) {
      // no atoms; create external_src of ions only
      num_reactions += 1;
    }
    // association
    if (!groups.neutral.empty() and !groups.molecule.empty()) {
      num_reactions += 1;
    }
  }

  this->reaction_options_store.reserve(num_reactions);
  this->reactions.reserve(num_reactions);

  for (const auto& element : this->base_species_order) {
    const auto& groups = this->species_by_base.at(element);

    // Add a pseudo-reaction if both `species1` and `species2` are present.
    // `gate_species`, if given, disables the reaction entirely while it is in deficit;
    // `boost_species`, if given, increases the rate more the further it is in deficit.
    //  If `boost_requires_deficit` is also set, the reaction is disabled entirely (rather than just left un-boosted) if `boost_species` is not
    // deficient.
    auto add_reaction_if_present =
        [&](const std::string& species1, const std::string& species2,
            PseudoReactionType type, const std::string& gate_species = "",
            const std::string& boost_species = "", bool boost_requires_deficit = false) {
          // Both species must be specified
          if (species1.empty() or species2.empty()) {
            return;
          }

          PseudoReactionSpec spec = make_reaction_spec(type, species1, species2);

          // Settings required by the Reaction base class need to be in an Options object
          this->reaction_options_store.emplace_back(this->options.copy());
          Options& reaction_options = this->reaction_options_store.back();
          reaction_options[PSEUDO_REACTION_COMPONENT_NAME]["data_srcs"] =
              toString(ReactionDataTypes::None);
          reaction_options[PSEUDO_REACTION_COMPONENT_NAME]["data_ids"] =
              toString(ReactionDataTypes::None);
          reaction_options[PSEUDO_REACTION_COMPONENT_NAME]["type"] = spec.reaction_str;
          reaction_options[PSEUDO_REACTION_COMPONENT_NAME]["diagnose"] = this->diagnose;
          /// This ensures that PseudoReaction instances don't interfere with the parsing of reaction strings from the config
          //  (See Reaction, ReactionBase constructors for details)
          reaction_options[PSEUDO_REACTION_COMPONENT_NAME]["is_internal"] = "true";

          // Configure all other settings in a PseudoReactionSpec struct
          const BoutReal gate_species_threshold =
              gate_species.empty() ? 0.0 : this->density_thresholds.at(gate_species);
          const BoutReal boost_species_threshold =
              boost_species.empty() ? 0.0 : this->density_thresholds.at(boost_species);
          spec.density_threshold = this->density_thresholds.at(spec.target_species);
          spec.timescale = this->timescale;
          spec.gate_species = gate_species;
          spec.gate_species_threshold = gate_species_threshold;
          spec.boost_species = boost_species;
          spec.boost_species_threshold = boost_species_threshold;
          spec.boost_requires_deficit = boost_requires_deficit;
          spec.external_src_mom_energy_enabled = this->external_src_mom_energy_enabled;
          spec.strategy = this->strategy;

          // Instantiate and add to the (pseudo-)reaction list
          this->reactions.emplace_back(
              std::make_unique<PseudoReaction>(reaction_options, spec));

          // transform_guarded() reuses LowNSources' own permission-checked state, so any
          // permissions a PseudoReaction sets on itself (e.g. when adding diagnostics) must be copied across
          for (const auto& [varname, region] :
               this->reactions.back()->getPermissions().getVariablesWithMinimumPermission(
                   PermissionTypes::Write)) {
            setPermissions(readWrite(varname, region));
          }
        };

    // Recombination restores the atom population and is boosted further when molecules are also depleted.
    add_reaction_if_present(groups.ion, groups.neutral, PseudoReactionType::recombination,
                            "", groups.molecule);
    // Ionization restores the ion population, but only when there are enough atoms to do so.
    add_reaction_if_present(groups.neutral, groups.ion, PseudoReactionType::ionization,
                            groups.neutral);
    // Association restores the molecule population.
    add_reaction_if_present(groups.neutral, groups.molecule,
                            PseudoReactionType::association);
    if (!groups.neutral.empty()) {
      // When atoms and ions are both present, the external ion source is enabled only when BOTH are depleted.
      // It is boosted by the level of atom depletion.
      add_reaction_if_present(groups.ion, groups.neutral,
                              PseudoReactionType::external_src, "", groups.neutral, true);
    } else {
      // With no atom species, the external ion source only responds to the ion deficit.
      add_reaction_if_present(groups.ion, groups.ion, PseudoReactionType::external_src);
    }
  }
}

void LowNSources::populate_reaction_defs() const {
  if (this->reaction_defs_populated) {
    return;
  }

  this->reaction_defs.clear();

  for (const auto& base_species : this->base_species_order) {
    const auto& groups = this->species_by_base.at(base_species);

    // Add reaction definitions.  If a required species isn't present, no definition is generated
    append_reaction_def_if_present(groups.ion, groups.neutral,
                                   PseudoReactionType::recombination,
                                   this->reaction_defs);
    append_reaction_def_if_present(groups.neutral, groups.ion,
                                   PseudoReactionType::ionization, this->reaction_defs);
    append_reaction_def_if_present(groups.neutral, groups.molecule,
                                   PseudoReactionType::association, this->reaction_defs);
    if (!groups.neutral.empty()) {
      append_reaction_def_if_present(groups.ion, groups.neutral,
                                     PseudoReactionType::external_src,
                                     this->reaction_defs);
    } else {
      append_reaction_def_if_present(
          groups.ion, groups.ion, PseudoReactionType::external_src, this->reaction_defs);
    }
  }

  this->reaction_defs_populated = true;
}

void LowNSources::transform_impl(GuardedOptions& state) {
  for (auto& reaction : this->reactions) {
    reaction->transform_guarded(state);
  }
}

void LowNSources::outputVars(Options& state) {
  for (auto& reaction : this->reactions) {
    reaction->outputVars(state);
  }
}

PseudoReaction::PseudoReaction(Options& options, const PseudoReactionSpec& spec)
    : Reaction(PSEUDO_REACTION_COMPONENT_NAME, options), type(spec.type),
      external_src_mom_energy_enabled(spec.external_src_mom_energy_enabled),
      density_threshold(spec.density_threshold), timescale(spec.timescale),
      target_species(spec.target_species), gate_species(spec.gate_species),
      gate_species_threshold(spec.gate_species_threshold),
      boost_species(spec.boost_species),
      boost_species_threshold(spec.boost_species_threshold),
      boost_requires_deficit(spec.boost_requires_deficit),
      strategy(LowNSrcsStrategyFromString(spec.strategy)) {

  if (this->diagnose) {
    // Set up a diagnostic for the density source of the target (produced) species
    std::string diag_name =
        fmt::format("S{:s}_low-n_{:s}", this->target_species, toString(this->type));
    std::string description = fmt::format("Particle source of {:s} due to low-n {:s}",
                                          this->target_species, toString(this->type));
    std::string src_name = fmt::format("low_n_{:s}", toString(this->type));
    add_diagnostic(this->target_species, diag_name, description,
                   ReactionDiagnosticType::density_src, src_name, identity);
  }
}

// ===================== PseudoReaction member functions ======================

RateData PseudoReaction::get_rate_standard(GuardedOptions& state) {

  // Get target species density
  // Local variable to workaround gcc dangling ref false positive
  const GuardedOptions n_target_opt = state["species"][this->target_species]["density"];
  const Field3D& n_target = n_target_opt.GetRef<Field3D>();
  auto region = n_target.getRegion("RGN_NOBNDRY");

  // Initialise return object. Pass an empty vector of reactants to avoid constructing unnecessary
  // collisions frequency fields, which aren't set by pseudo-reactions.
  RateData result({});
  result.rate = zeroFrom(n_target);

  // Set rate such that the source term will restore the target species to its threshold density over delta_t = [timescale]
  //   (deficit = source = 0  when n_target > n_thresh)
  BOUT_FOR(i, region) {
    const BoutReal n_deficit =
        std::max<BoutReal>(0.0, this->density_threshold - n_target[i]);
    result.rate[i] = n_deficit / this->timescale;
  }

  // If this reaction has a gate species, smoothly reduce the rate towards zero as it
  // becomes depleted (a hard on/off step here made the Newton solver's job harder).
  if (!this->gate_species.empty()) {
    // Local variable to workaround gcc dangling ref false positive
    const GuardedOptions n_gate_opt = state["species"][this->gate_species]["density"];
    const Field3D& n_gate = n_gate_opt.GetRef<Field3D>();
    BOUT_FOR(i, region) {
      result.rate[i] *= smooth_gate(n_gate[i] / this->gate_species_threshold, 0.1);
    }
  }

  // If this reaction has a boost species, modify the rate based on the level of its depletion.
  // Disable reactions entirely if a boost species deficit is required (i.e. for the external ion source)
  if (!this->boost_species.empty()) {
    // Local variable to workaround gcc dangling ref false positive
    const GuardedOptions n_boost_opt = state["species"][this->boost_species]["density"];
    const Field3D& n_boost = n_boost_opt.GetRef<Field3D>();
    BOUT_FOR(i, region) {
      const BoutReal thresh_ratio = n_boost[i] / this->boost_species_threshold;
      if (thresh_ratio < 1.0) {
        // The more deficient the boost species is, the further the rate is boosted
        // (Smoothly floor the inverse boost factor at 0.01 to avoid a divide-by-zero
        // and the derivative discontinuity a hard std::clamp would introduce)
        const BoutReal inv_boost_fac = smooth_clamp(thresh_ratio, 0.01, 1.0, 0.01);
        result.rate[i] /= inv_boost_fac;
      } else if (this->boost_requires_deficit) {
        // Some reactions (e.g. an external ion source) should only run at all while
        // the boost species is deficient, rather than merely running un-boosted.
        result.rate[i] = 0.0;
      }
    }
  }

  // Set flag if all (rank-local) rates are zero so that Reaction transform can skip source term construction
  result.is_zero = std::all_of(region.begin(), region.end(),
                               [&](const auto& i) { return result.rate[i] == 0.0; });
  return result;
}

RateData PseudoReaction::get_rate(GuardedOptions& state) {
  switch (this->strategy) {
  case LowNSrcsStrategy::standard:
    return get_rate_standard(state);
  default:
    throw std::runtime_error("Unknown low-n sources strategy [" + toString(this->strategy)
                             + "]");
  }
}

/**
 * @brief Add momentum and energy for external ion sources.
 *
 * @details The external ion source abuses the reaction framework using a string like "ion -> 2ion".
 * Since the only reactant has a +ve population change, Reaction::transform won't generate momentum
 * or energy sources. Add them here instead, unless they've been explicitly disabled.
 *
 */
void PseudoReaction::transform_additional(GuardedOptions& state,
                                          const RateData& rate_calc_results) {
  // Momentum and energy sources (for the external ion source only)
  if (this->type != PseudoReactionType::external_src
      || !this->external_src_mom_energy_enabled) {
    return;
  }
  const auto species = state["species"][this->target_species];
  const Field3D momentum_source = rate_calc_results.rate * get<BoutReal>(species["AA"])
                                  * get<Field3D>(species["velocity"]);
  const Field3D energy_source =
      rate_calc_results.rate * (3. / 2) * get<Field3D>(species["temperature"]);
  update_source<add<Field3D>>(state, this->target_species,
                              ReactionDiagnosticType::momentum_src, momentum_source);
  update_source<add<Field3D>>(state, this->target_species,
                              ReactionDiagnosticType::energy_src, energy_source);
}

} // namespace hermes
