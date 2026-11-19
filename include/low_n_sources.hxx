#pragma once

#ifndef LOW_N_SOURCES_H
#define LOW_N_SOURCES_H

#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <bout/bout_enum_class.hxx>
#include <bout/field3d.hxx>

#include "component.hxx"
#include "species_parser.hxx"

/// An enum representing the different possible pseudo-reaction types.
BOUT_ENUM_CLASS(PseudoReactionType, association, recombination, ionization, external_src)
namespace hermes {

// Forward declare PseudoReaction to avoid exposing its full interface in this header.
class PseudoReaction;

/**
 * @brief Struct that groups together related ion, neutral, and molecule species for use in
 * constructing pseudo-reactions.
 *
 * @details Groups can be either element-based (ion, atom, molecule) or molecule-based
 * (molecular-ion, neutral molecule, empty). Only species present in the simulation are populated.
 */
struct PseudoReactionSpeciesGroup {
  /// Ion species (atomic or molecular)
  std::string ion;
  /// Neutral species (atom or molecule)
  std::string neutral;
  /// Molecule species (only populated for element-based groups)
  std::string molecule;
};

/**
 * @brief Component to handle addition of sources/sinks that prevent low density for species
 * in a simulation.
 */
struct LowNSources : public NamedComponent<LowNSources> {
public:
  static constexpr auto type = "low_n_sources";

  /// Species label(s) and reaction type used to define a pseudo-reaction.
  using ReactionDef = std::tuple<std::string, std::string, PseudoReactionType>;

  /// Read options and assemble reactions
  LowNSources(const std::string& name, Options& options);
  /// Constructor used by the Component-factory; solver arg is discarded.
  LowNSources(const std::string& name, Options& options, Solver* solver);
  ~LowNSources() override;

  // Don't allow copying or moving of LowNSources instances
  LowNSources(const LowNSources&) = delete;
  LowNSources& operator=(const LowNSources&) = delete;
  LowNSources(LowNSources&&) = delete;
  LowNSources& operator=(LowNSources&&) = delete;

  std::vector<ReactionDef> get_reaction_defs() const;

  void outputVars(Options& state) final;

private:
  Options& options;

  // Whether or not to output diagnostics for the low-n sources
  bool diagnose;

  // Control params
  /// Whether to enable momentum and energy sources for external ion sources
  bool external_src_mom_energy_enabled;
  /// Reactions kick in when density is below (species density_floor)*floor_factor
  BoutReal floor_factor;
  /// Deficit-relaxation timescale in normalised units, before gates and boosts
  BoutReal timescale;

  /// Per-species density thresholds, used to enable/disable sources
  std::map<std::string, BoutReal> density_thresholds;

  /// Species names extracted from the component list
  std::vector<std::string> all_species;

  /// Map species to their respective groups
  std::map<std::string, PseudoReactionSpeciesGroup> species_groups;
  /// First-seen order of species groups to guarantee deterministic reaction ordering
  std::vector<std::string> species_groups_order;

  /// Reaction definitions, populated by populate_reaction_defs() and used to instantiate PseudoReaction objects
  mutable std::vector<ReactionDef> reaction_defs;
  mutable bool reaction_defs_populated{false};

  /// Own the options objects passed to each PseudoReaction so referenced unit
  /// sections remain valid for the reaction lifetime.
  std::vector<Options> reaction_options_store;
  std::vector<std::unique_ptr<PseudoReaction>> reactions;

  /// Strategy to use when generating sources
  std::string strategy;

  /// Add species to their respective groups, rejecting duplicate entries in each category.
  void append_species_by_type(const SpeciesParser& species, SpeciesType type);

  void instantiate_reactions();
  /// Rebuild the cached definitions if a new species grouping has invalidated them.
  void populate_reaction_defs() const;
  void transform_impl(GuardedOptions& state) final;
};

} // namespace hermes
namespace {
RegisterComponent<hermes::LowNSources> registercomponentlownsource;
}

#endif // LOW_N_SOURCES_H
