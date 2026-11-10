#include "low_n_sources.hxx"

#include <algorithm>
#include <string>

#include <bout/boutexception.hxx>
#include <bout/options.hxx>
#include <bout/output.hxx>
#include <bout/utils.hxx>

#include "guarded_options.hxx"
#include "pseudo_reaction.hxx"
#include "species_parser.hxx"

// ======================= Helper struct and functions ========================
namespace {

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
 * @brief Assemble the information needed to define a pseudo-reaction.
 *
 * @param type The type of pseudo-reaction
 * @param species_group Group containing related ion, neutral, and molecule species. Some or all of
 * the species may be involved in the reaction.
 * @return PseudoReactionSpec The reaction specification
 */
hermes::PseudoReactionSpec
make_reaction_spec(PseudoReactionType type,
                   const hermes::PseudoReactionSpeciesGroup& species_group) {
  hermes::PseudoReactionSpec spec;
  spec.type = type;

  switch (type) {
  case PseudoReactionType::association: {
    // Allow for molecules made up of any number of atoms (of the same element)
    const SpeciesParser molecule(species_group.molecule);
    const int atoms_per_molecule =
        stringToInt(molecule.get_base_species().substr(molecule.get_element().size()));
    // N.B. SpeciesParser constructor rejects molecule strings with <1 atoms, so we can assume
    // atoms_per_molecule is always +ve
    spec.reaction_str = fmt::format("{}{} -> {}", atoms_per_molecule,
                                    species_group.neutral, species_group.molecule);
    spec.target_species = species_group.molecule;
    return spec;
  }
  case PseudoReactionType::recombination:
    spec.reaction_str =
        fmt::format("{} + e -> {}", species_group.ion, species_group.neutral);
    spec.target_species = species_group.neutral;
    spec.boost_species = species_group.molecule;
    return spec;
  case PseudoReactionType::ionization:
    spec.reaction_str =
        fmt::format("{} -> {} + e", species_group.neutral, species_group.ion);
    spec.target_species = species_group.ion;
    spec.gate_species = species_group.neutral;
    return spec;
  case PseudoReactionType::external_src:
    // N.B. Since no species in this reaction has a -ve pop. change, Reaction::transform_impl doesn't add momentum and energy sources.
    // These are added in PseudoReaction::transform_additional instead, unless turned off with external_src_mom_energy_enabled=false.
    spec.reaction_str = fmt::format("{} -> 2{}", species_group.ion, species_group.ion);
    spec.target_species = species_group.ion;
    if (!species_group.neutral.empty()) {
      spec.boost_species = species_group.neutral;
      spec.boost_requires_deficit = true;
    }
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
LowNSources::~LowNSources() = default;

/// Constructor used by the Component-factory; solver arg is discarded.
LowNSources::LowNSources(const std::string& name, Options& options,
                         [[maybe_unused]] Solver* solver)
    : LowNSources(name, options) {}

/**
 * @brief Main constructor: Read options and set up pseudo-reactions.
 */
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

/**
 * @brief Add species to their respective groups, rejecting duplicate entries in each category.
 */
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
      this->species_groups.try_emplace(element_str);
  if (element_group_inserted) {
    this->species_groups_order.push_back(element_str);
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
    auto [mol_group_it, mol_group_inserted] = this->species_groups.try_emplace(mol_str);
    if (mol_group_inserted) {
      this->species_groups_order.push_back(mol_str);
    }
    auto& mol_group = this->species_groups.at(mol_str);
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
 * @return A copy of the cached reaction definitions
 */
std::vector<LowNSources::ReactionDef> LowNSources::get_reaction_defs() const {
  populate_reaction_defs();
  return this->reaction_defs;
}

void LowNSources::instantiate_reactions() {
  this->reaction_options_store.clear();
  this->reactions.clear();

  populate_reaction_defs();
  this->reaction_options_store.reserve(this->reaction_defs.size());
  this->reactions.reserve(this->reaction_defs.size());

  /*
    Loop over the reaction definitions, instantiating a PseudoReaction instance for each one.
    Note that the existence of a definition at this point guarantees that all required heavy species
    are present in the simulation.
    */
  for (const auto& [species1, species2, type] : this->reaction_defs) {
    const auto& species =
        this->species_groups.at(SpeciesParser(species1).get_base_species());
    PseudoReactionSpec spec = make_reaction_spec(type, species);

    // Settings required by the Reaction base class need to be in an Options object
    this->reaction_options_store.emplace_back(this->options.copy());
    Options& reaction_options = this->reaction_options_store.back();
    reaction_options[PseudoReaction::component_name]["data_srcs"] =
        toString(ReactionDataTypes::None);
    reaction_options[PseudoReaction::component_name]["data_ids"] =
        toString(ReactionDataTypes::None);
    reaction_options[PseudoReaction::component_name]["type"] = spec.reaction_str;
    reaction_options[PseudoReaction::component_name]["diagnose"] = this->diagnose;

    // Configure all other settings in a PseudoReactionSpec struct
    spec.density_threshold = this->density_thresholds.at(spec.target_species);
    spec.timescale = this->timescale;
    spec.gate_species_threshold =
        spec.gate_species.empty() ? 0.0 : this->density_thresholds.at(spec.gate_species);
    spec.boost_species_threshold = spec.boost_species.empty()
                                       ? 0.0
                                       : this->density_thresholds.at(spec.boost_species);
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
  }
}

/**
 * @brief Rebuild the cached reaction definitions if a new species grouping has invalidated them.
 */
void LowNSources::populate_reaction_defs() const {
  if (this->reaction_defs_populated) {
    return;
  }

  this->reaction_defs.clear();

  for (const auto& base_species : this->species_groups_order) {
    const auto& groups = this->species_groups.at(base_species);

    // Add reaction definitions.  If a required species isn't present, no definition is generated

    // Recombination restores neutrals; for an atomic product, the rate is boosted further when molecules are also depleted.
    append_reaction_def_if_present(groups.ion, groups.neutral,
                                   PseudoReactionType::recombination,
                                   this->reaction_defs);
    // Ionization restores the ion population, with smooth suppression when neutrals are depleted.
    append_reaction_def_if_present(groups.neutral, groups.ion,
                                   PseudoReactionType::ionization, this->reaction_defs);
    // Association restores the molecule population.
    append_reaction_def_if_present(groups.neutral, groups.molecule,
                                   PseudoReactionType::association, this->reaction_defs);
    if (!groups.neutral.empty()) {
      // If both neutrals and ions are present, the external ion source is enabled only when BOTH are depleted.
      // It is boosted by the level of neutral depletion.
      append_reaction_def_if_present(groups.ion, groups.neutral,
                                     PseudoReactionType::external_src,
                                     this->reaction_defs);
    } else {
      append_reaction_def_if_present(
          // With no neutral partner, the external ion source only responds to the ion deficit.
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

} // namespace hermes
