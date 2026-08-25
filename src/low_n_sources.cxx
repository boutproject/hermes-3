#include "low_n_sources.hxx"

#include <string>

#include <bout/boutexception.hxx>
#include <bout/options.hxx>
#include <bout/utils.hxx>

#include "guarded_options.hxx"
#include "single_rate_reaction_data.hxx"
#include "species_parser.hxx"

// ======================= Helper struct and functions ========================
namespace {

/// The information needed to define a pseudo-reaction
using PseudoReactionSpec = hermes::PseudoReactionSpec;

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
  case PseudoReactionType::association:
    spec.reaction_str = fmt::format("2{} -> {}", species1, species2);
    spec.target_species = species2;
    return spec;
  case PseudoReactionType::recombination:
    spec.reaction_str = fmt::format("{} + e -> {}", species1, species2);
    spec.target_species = species2;
    return spec;
  case PseudoReactionType::ionization:
    spec.reaction_str = fmt::format("{} -> {} + e", species1, species2);
    spec.target_species = species2;
    return spec;
  case PseudoReactionType::external_src:
    spec.reaction_str = fmt::format("{} -> 2{}", species1, species1);
    spec.target_species = species1;
    return spec;
  }
  throw BoutException("Unhandled PseudoReactionType in make_reaction_spec()");
}

/// @brief Convenience function to set a unique species or throw an exception if one is already set
void set_unique_species_or_throw(std::string& existing_species,
                                 const std::string& new_species,
                                 const std::string& element,
                                 const std::string& species_category) {
  if (!existing_species.empty()) {
    throw BoutException(
        fmt::format("LowNSources supports at most one {} species per element. "
                    "Element '{}' has '{}' and '{}'.",
                    species_category, element, existing_species, new_species));
  }
  existing_species = new_species;
}

} // namespace

namespace hermes {

// ======================= LowNSources member functions =======================
LowNSources::LowNSources(std::string name, Options& options,
                         [[maybe_unused]] Solver* solver)
    : LowNSources(std::move(name), options) {}

LowNSources::LowNSources(std::string name, Options& options)
    : NamedComponent(name, {readOnly("species:{all_species}:{read_vals}"),
                            readWrite("species:{all_species}:{write_vals}")}),
      options(options), atom_species{}, ion_species{}, molecule_species{} {
  substitutePermissions("read_vals", {"AA", "density", "velocity", "temperature"});
  substitutePermissions("write_vals", {"density_source", "momentum_source",
                                       "energy_source", "collision_frequency"});

  // Read the list of components
  const std::string component_names = options["hermes"]["components"]
                                          .doc("Components in order of execution")
                                          .withDefault<std::string>("");

  // Read control parameters
  /// Floor factor
  this->low_n_floor_fac = options[name]["low_n_floor_fac"]
                              .doc("Multiplier applied to the low-density floor")
                              .withDefault<BoutReal>(10.0);
  Options& root_options = Options::root();
  /// Timescale over which low density sources act
  BoutReal default_timescale = 10 * root_options["timestep"].withDefault<BoutReal>(1.0);
  this->low_n_timescale = options[name]["low_n_timescale"]
                              .doc("Timescale over which low density sources act")
                              .withDefault<BoutReal>(default_timescale);
  /// Rate calculation strategy
  this->strategy =
      options[name]["strategy"]
          .doc("Strategy to use when generating sources to avoid low species density.")
          .withDefault(toString(LowNSrcsStrategy::standard));
  // Diagnostics output switch
  this->diagnose = options[name]["diagnose"]
                       .doc("Whether or not to output diagnostics for the low-n sources")
                       .withDefault<bool>(false);

  // Classify species by element and by type
  std::vector<std::string> species_with_mass;
  for (const auto& component : strsplit(component_names, ',')) {
    const std::string component_name = trim(component, " \t\r()");
    if (component_name.empty() or !options[component_name].isSet("AA")) {
      continue;
    }
    species_with_mass.push_back(component_name);
    SpeciesParser species(component_name);
    append_species_by_type(species, SpeciesParser::get_type(species));
  }

  // LowNSources can be instantiated directly in unit tests, bypassing scheduler-level
  // declareAllSpecies() substitutions, so resolve placeholders here.
  substitutePermissions("all_species", species_with_mass);

  // Set activation thresholds for each species now that species lists are populated
  for (const auto& species_group :
       {&this->ion_species, &this->atom_species, &this->molecule_species}) {
    for (const auto& species_name : *species_group) {
      const BoutReal density_floor = options[species_name]["density_floor"]
                                         .doc("Minimum density floor")
                                         .withDefault<BoutReal>(1e-7);
      this->density_thresholds[species_name] = this->low_n_floor_fac * density_floor;
    }
  }

  instantiate_reactions();
}

void LowNSources::append_species_by_type(const SpeciesParser& species,
                                         SpeciesTypeAlt type) {
  const std::string& species_name = species.get_str();
  const std::string element = species.get_element();

  auto [it, inserted] = this->species_by_element.try_emplace(element);
  if (inserted) {
    this->element_order.push_back(element);
  }

  auto& groups = it->second;
  if (type == SpeciesTypeAlt::ion or type == SpeciesTypeAlt::molecular_ion) {
    set_unique_species_or_throw(groups.ion, species_name, element, "ion");
    this->ion_species.push_back(species_name);
  } else if (type == SpeciesTypeAlt::neutral) {
    set_unique_species_or_throw(groups.atom, species_name, element, "atom");
    this->atom_species.push_back(species_name);
  } else if (type == SpeciesTypeAlt::molecule) {
    set_unique_species_or_throw(groups.molecule, species_name, element, "molecule");
    this->molecule_species.push_back(species_name);
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
  for (const auto& element : this->element_order) {
    const auto& groups = this->species_by_element.at(element);
    if (!groups.ion.empty() and !groups.atom.empty()) {
      // recombination, ionization, external_src of ions
      num_reactions += 3;
    } else if (!groups.ion.empty()) {
      // no atoms; create external_src of ions only
      num_reactions += 1;
    }
    // association
    if (!groups.atom.empty() and !groups.molecule.empty()) {
      num_reactions += 1;
    }
  }

  this->reaction_options_store.reserve(num_reactions);
  this->reactions.reserve(num_reactions);

  for (const auto& element : this->element_order) {
    const auto& groups = this->species_by_element.at(element);

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
              SINGLE_RATE_REACTION_DATA_LBL;
          reaction_options[PSEUDO_REACTION_COMPONENT_NAME]["data_ids"] = "pseudo_rate";
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
          spec.low_n_timescale = this->low_n_timescale;
          spec.gate_species = gate_species;
          spec.gate_species_threshold = gate_species_threshold;
          spec.boost_species = boost_species;
          spec.boost_species_threshold = boost_species_threshold;
          spec.boost_requires_deficit = boost_requires_deficit;
          spec.strategy = this->strategy;

          // Instantiate and add to the (pseudo-)reaction list
          this->reactions.emplace_back(
              std::make_unique<PseudoReaction>(reaction_options, spec));

          // transform_guarded() reuses LowNSources' own permission-checked state, so any
          // permissions a pseudoreaction sets on itself (e.g. when adding diagnostics) must be copied across
          for (const auto& [varname, region] :
               this->reactions.back()->getPermissions().getVariablesWithPermission(
                   PermissionTypes::Write, false)) {
            setPermissions(readWrite(varname, region));
          }
        };

    // Recombination restores the atom population and is boosted further when molecules are also depleted.
    add_reaction_if_present(groups.ion, groups.atom, PseudoReactionType::recombination,
                            "", groups.molecule);
    // Ionization restores the ion population, but only when there are enough atoms to do so.
    add_reaction_if_present(groups.atom, groups.ion, PseudoReactionType::ionization,
                            groups.atom);
    // Association restores the molecule population.
    add_reaction_if_present(groups.atom, groups.molecule,
                            PseudoReactionType::association);
    if (!groups.atom.empty()) {
      // When atoms and ions are both present, the external ion source is enabled only when BOTH are depleted.
      // It is boosted by the level of atom depletion.
      add_reaction_if_present(groups.ion, groups.atom, PseudoReactionType::external_src,
                              "", groups.atom, true);
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

  for (const auto& element : this->element_order) {
    const auto& groups = this->species_by_element.at(element);

    // Add reaction definitions.  If a required species isn't present, no definition is generated
    append_reaction_def_if_present(
        groups.ion, groups.atom, PseudoReactionType::recombination, this->reaction_defs);
    append_reaction_def_if_present(groups.atom, groups.ion,
                                   PseudoReactionType::ionization, this->reaction_defs);
    append_reaction_def_if_present(groups.atom, groups.molecule,
                                   PseudoReactionType::association, this->reaction_defs);
    if (!groups.atom.empty()) {
      append_reaction_def_if_present(
          groups.ion, groups.atom, PseudoReactionType::external_src, this->reaction_defs);
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
    : Reaction(PSEUDO_REACTION_COMPONENT_NAME, options) {
  this->type = spec.type;
  this->density_threshold = spec.density_threshold;
  this->low_n_timescale = spec.low_n_timescale;
  this->target_species = spec.target_species;
  this->gate_species = spec.gate_species;
  this->gate_species_threshold = spec.gate_species_threshold;
  this->boost_species = spec.boost_species;
  this->boost_species_threshold = spec.boost_species_threshold;
  this->boost_requires_deficit = spec.boost_requires_deficit;
  this->strategy = LowNSrcsStrategyFromString(spec.strategy);

  if (this->diagnose) {
    // Set up a diagnostic for the density source of the target (produced) species
    std::string label =
        fmt::format("S{:s}_low-n_{:s}", this->target_species, toString(this->type));
    std::string description = fmt::format("Particle source of {:s} due to low-n {:s}",
                                          toString(this->type), this->target_species);
    std::string src_name = fmt::format("low_n_{:s}", toString(this->type));
    add_diagnostic(this->target_species, label, description,
                   ReactionDiagnosticType::density_src, src_name, identity);
  }
}

// ===================== PseudoReaction member functions ======================
RateData PseudoReaction::get_rate_standard(GuardedOptions& state) {

  // Get target species density
  const Field3D& n_target =
      state["species"][this->target_species]["density"].GetRef<Field3D>();
  auto region = n_target.getRegion("RGN_NOBNDRY");

  // Initialise return object
  RateData result(this->reactant_names);
  result.rate = zeroFrom(n_target);

  // Collision frequencies aren't used in PseudoReaction; set them to zero
  for (auto& coll_freq : result.collision_frequencies) {
    coll_freq = zeroFrom(n_target);
  }

  // Set rate such that the source term will restore the target species to its threshold density over delta_t = [low_n_timescale]
  //   (deficit = source = 0  when n_target > n_thresh)
  BOUT_FOR(i, region) {
    const BoutReal n_deficit =
        std::max<BoutReal>(0.0, this->density_threshold - n_target[i]);
    result.rate[i] = n_deficit / this->low_n_timescale;
  }

  // If this reaction has a gate species, zero the rate when it's depleted.
  if (!this->gate_species.empty()) {
    const Field3D& n_gate =
        state["species"][this->gate_species]["density"].GetRef<Field3D>();
    BOUT_FOR(i, region) {
      if (n_gate[i] < this->gate_species_threshold) {
        result.rate[i] = 0.0;
      }
    }
  }

  // If this reaction has a boost species, modify the rate based on the level of its depletion.
  // Disable reactions entirely if a boost species deficit is required (i.e. for the external ion source)
  if (!this->boost_species.empty()) {
    const Field3D& n_boost =
        state["species"][this->boost_species]["density"].GetRef<Field3D>();
    BOUT_FOR(i, region) {
      const BoutReal thresh_ratio = n_boost[i] / this->boost_species_threshold;
      if (thresh_ratio < 1.0) {
        // The more deficient the boost species is, the further the rate is boosted
        // (Avoid divide-by-zero by enforcing inverse boost factor >= 0.01)
        const BoutReal inv_boost_fac = std::clamp(thresh_ratio, 0.01, 1.0);
        result.rate[i] /= inv_boost_fac;
      } else if (this->boost_requires_deficit) {
        // Some reactions (e.g. an external ion source) should only run at all while
        // the boost species is deficient, rather than merely running un-boosted.
        result.rate[i] = 0.0;
      }
    }
  }

  return result;
}

RateData PseudoReaction::get_rate_uedge([[maybe_unused]] GuardedOptions& state) {
  throw BoutException("PseudoReaction::get_rate_uedge() not yet implemented");
}

RateData PseudoReaction::get_rate(GuardedOptions& state) {
  switch (this->strategy) {
  case LowNSrcsStrategy::standard:
    return get_rate_standard(state);
  case LowNSrcsStrategy::uedge:
    return get_rate_uedge(state);
  default:
    throw std::runtime_error("Unknown low-n sources strategy [" + toString(this->strategy)
                             + "]");
  }
}

} // namespace hermes
