#pragma once

#ifndef LOW_N_SOURCES_H
#define LOW_N_SOURCES_H

#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <bout/bout_enum_class.hxx>
#include <bout/field3d.hxx>

#include "component.hxx"
#include "reaction.hxx"
#include "species_parser.hxx"

BOUT_ENUM_CLASS(LowNSrcsStrategy, uedge, standard)
BOUT_ENUM_CLASS(PseudoReactionType, association, recombination, ionization, external_src)
namespace hermes {

// Forward declare PseudoReaction so we can keep the LowNSources interface at the top of this file
class PseudoReaction;

/**
 * @brief Component to handle addition of sources/sinks that prevent low density for species
 * in a simulation.
 */
struct LowNSources : public NamedComponent<LowNSources> {
public:
  static constexpr auto type = "low_n_sources";

  using ReactionDef = std::tuple<std::string, std::string, PseudoReactionType>;

  LowNSources(std::string name, Options& options);
  LowNSources(std::string name, Options& options, Solver* solver);

  std::vector<ReactionDef> get_reaction_defs() const;

private:
  struct ElementSpeciesGroups {
    std::string ion;
    std::string atom;
    std::string molecule;
  };

  Options& options;

  // Whether or not to output diagnostics for the low-n sources
  bool diagnose;

  // Control params
  /// Reactions kick in when density is at or below (species density_floor)*low_n_floor_fac
  BoutReal low_n_floor_fac;
  /// Timescale over which sources aim to return density to the threshold value
  BoutReal low_n_timescale;

  /// Per-species density thresholds, used to enable/disable sources
  std::map<std::string, BoutReal> density_thresholds;

  /// Species names extracted from the component list
  std::vector<std::string> atom_species;
  std::vector<std::string> ion_species;
  std::vector<std::string> molecule_species;

  /// Species grouped by element
  std::map<std::string, ElementSpeciesGroups> species_by_element;
  /// First-seen order of elements for deterministic reaction ordering
  std::vector<std::string> element_order;

  /// Reaction definitions, populated by populate_reaction_defs() and used to instantiate PseudoReaction objects
  mutable std::vector<ReactionDef> reaction_defs;
  mutable bool reaction_defs_populated{false};

  /// Own the options objects passed to each PseudoReaction so referenced unit
  /// sections remain valid for the reaction lifetime.
  std::vector<Options> reaction_options_store;
  std::vector<std::unique_ptr<PseudoReaction>> reactions;

  /// Strategy to use when generating sources
  std::string strategy;

  void append_species_by_type(const SpeciesParser& species, SpeciesTypeAlt type);

  void instantiate_reactions();
  void outputVars(Options& state) final;
  void populate_reaction_defs() const;
  void transform_impl(GuardedOptions& state) final;
};

const std::string PSEUDO_REACTION_COMPONENT_NAME = "pseudo_reaction";

struct PseudoReactionSpec {
  std::string reaction_str;
  PseudoReactionType type;
  std::string target_species;
  BoutReal density_threshold{0.0};
  BoutReal low_n_timescale{0.0};
  std::string gate_species;
  BoutReal gate_species_threshold{0.0};
  std::string boost_species;
  BoutReal boost_species_threshold{0.0};
  bool boost_requires_deficit{false};
  std::string strategy;
};

class PseudoReaction : public Reaction {
public:
  /**
   * @brief Construct a new Pseudo Reaction object
   *
   * @param options Options object used by Reaction base class
   * @param spec Typed pseudo-reaction configuration
   */
  PseudoReaction(Options& options, const PseudoReactionSpec& spec);

protected:
  /// Override Reaction::get_rate to bypass RateHelper and use custom rates
  RateData get_rate(GuardedOptions& state) override;

private:
  /// Type of the pseudo-reaction
  PseudoReactionType type;
  /// Threshold density below which the source term is activated
  BoutReal density_threshold;
  /// Timescale over which the source term would restore density to the threshold value
  BoutReal low_n_timescale;
  /// Species whose density may be boosted by the reaction
  std::string target_species;
  /// Species whose depletion will disable this reaction entirely (empty string means no gating)
  std::string gate_species;
  /// Threshold density below which gate_species disables the reaction
  BoutReal gate_species_threshold;
  /// Species whose depletion boosts this reaction's rate (empty string means no boost)
  std::string boost_species;
  /// Threshold density below which boost_species enhances the reaction rate
  BoutReal boost_species_threshold;
  /// If true, the reaction is disabled entirely if boost_species is not deficient
  bool boost_requires_deficit;
  /// Strategy to use when generating sources
  LowNSrcsStrategy strategy;

  /// Calculate rate using 'standard' strategy (default)
  RateData get_rate_standard(GuardedOptions& state);
  /// Calculate rate using 'uedge' strategy
  RateData get_rate_uedge(GuardedOptions& state);
};

} // namespace hermes
namespace {
RegisterComponent<hermes::LowNSources> registercomponentlownsource;
}

#endif // LOW_N_SOURCES_H
