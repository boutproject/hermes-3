#pragma once

#ifndef PSEUDO_REACTION_H
#define PSEUDO_REACTION_H

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include <bout/bout_enum_class.hxx>

#include "low_n_sources.hxx"
#include "reaction.hxx"

// An enum representing different potential strategies for low-n sources; only standard is implemented.
// (BOUT_ENUM_CLASS specialises Options::as(), so must be in the global scope)
BOUT_ENUM_CLASS(LowNSrcsStrategy, uedge, standard)

namespace hermes {

/// A struct to aggregate the information needed to define a pseudo-reaction
struct PseudoReactionSpec {
  std::string reaction_str;
  PseudoReactionType type;
  std::string target_species;
  BoutReal density_threshold{0.0};
  BoutReal timescale{0.0};
  std::string gate_species;
  BoutReal gate_species_threshold{0.0};
  std::string boost_species;
  BoutReal boost_species_threshold{0.0};
  bool boost_requires_deficit{false};
  bool external_src_mom_energy_enabled{true};
  std::string strategy;
};

/**
 * @brief Class representing a pseudo-reaction used in `LowNSources`.
 *
 * @details Uses the `Reaction` class to create sources that restore target species towards a threshold density.
 * Rates depend on the deficit of the target species and can be modified by 'gate' and 'boost' species.
 *
 */
class PseudoReaction : public Reaction {
public:
  static constexpr auto component_name = "pseudo_reaction";

  PseudoReaction(Options& options, const PseudoReactionSpec& spec);

protected:
  ///
  RateData get_rate(GuardedOptions& state) override;

  ///
  void transform_additional(GuardedOptions& state,
                            const RateData& rate_calc_results) override;

private:
  /// Type of the pseudo-reaction
  PseudoReactionType type;
  /// Whether external ion sources add momentum and energy sources as well as density
  bool external_src_mom_energy_enabled;
  /// Threshold density below which the source term is activated
  BoutReal density_threshold;
  /// Deficit-reduction timescale in normalised units, before gates and boosts
  BoutReal timescale;
  /// Species whose density is restored by the reaction
  std::string target_species;
  /// Species whose depletion smoothly suppresses this reaction (empty string means no gating)
  std::string gate_species;
  /// Density at and below which the gate weight is zero
  BoutReal gate_species_threshold;
  /// Species whose depletion boosts this reaction's rate (empty string means no boost)
  std::string boost_species;
  /// Threshold density below which boost_species enhances the reaction rate
  BoutReal boost_species_threshold;
  /// If true, the reaction is disabled entirely if boost_species is not deficient
  bool boost_requires_deficit;
  /// Strategy to use when generating sources
  LowNSrcsStrategy strategy;

  static BoutReal smooth_clamp(BoutReal x, BoutReal lo, BoutReal hi, BoutReal width);
  static BoutReal smooth_gate(BoutReal ratio, BoutReal width);

  ///
  RateData get_rate_standard(GuardedOptions& state);
};

/**
 * @brief Map `x` to a smooth sigmoid between `lo` and `hi`.
 * The transition is centred on (`lo` +`hi`) / 2.
 *
 * @param x The value to map
 * @param lo Lower bound, approached smoothly as x -> -inf
 * @param hi Upper bound, approached smoothly as x -> +inf
 * @param width Transition scale; smaller is closer to a hard step
 */
inline BoutReal PseudoReaction::smooth_clamp(BoutReal x, BoutReal lo, BoutReal hi,
                                             BoutReal width) {
  return lo + (hi - lo) / (1.0 + exp(-(x - 0.5 * (lo + hi)) / width));
}

/**
 * @brief Smoothly transition from 0 (ratio <= 1) to 1 (ratio >= 1 + width),
 * with continuous first and second derivatives.
 *
 * @param ratio A value expected to vary around 1 (e.g. density / density_threshold)
 * @param width Positive transition width (above ratio = 1)
 */
inline BoutReal PseudoReaction::smooth_gate(BoutReal ratio, BoutReal width) {
  const BoutReal x = std::clamp((ratio - 1.0) / width, 0.0, 1.0);
  return x * x * x * (10.0 + x * (-15.0 + 6.0 * x));
}

/**
 * @brief Construct a new Pseudo Reaction object
 *
 * @param options Options object
 * @param spec Typed pseudo-reaction configuration
 */
inline PseudoReaction::PseudoReaction(Options& options, const PseudoReactionSpec& spec)
    : Reaction(component_name, options), type(spec.type),
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

/**
 * @brief Calculate a rate using the 'standard' strategy.
 */
inline RateData PseudoReaction::get_rate_standard(GuardedOptions& state) {

  // Get target species density
  // Local variable to workaround gcc dangling ref false positive
  const GuardedOptions n_target_opt = state["species"][this->target_species]["density"];
  const Field3D& n_target = n_target_opt.GetRef<Field3D>();
  auto region = n_target.getRegion("RGN_NOBNDRY");

  // Initialise return object. Pass an empty vector of reactants to avoid constructing unnecessary
  // collisions frequency fields, which aren't set by pseudo-reactions.
  RateData result({});
  result.rate = zeroFrom(n_target);

  // Relax the target deficit with e-folding time [timescale], before gates and boosts.
  //   (deficit = source = 0 when n_target >= n_thresh)
  BOUT_FOR(i, region) {
    const BoutReal n_deficit =
        std::max<BoutReal>(0.0, this->density_threshold - n_target[i]);
    result.rate[i] = n_deficit / this->timescale;
  }

  // If this reaction has a gate species, smoothly reduce the rate to zero at its threshold
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
        // (Use a sigmoid between 0.01 and 1.0 for the inverse boost factor,
        // bounding the boost between 1 and 100 and avoiding division by zero)
        const BoutReal inv_boost_fac = smooth_clamp(thresh_ratio, 0.01, 1.0, 0.01);
        result.rate[i] /= inv_boost_fac;
      }
      if (this->boost_requires_deficit) {
        // Some reactions (e.g. an external ion source) should only run at all while
        // the boost species is deficient, rather than merely running un-boosted.
        result.rate[i] *= smooth_gate(2.0 - thresh_ratio, 0.1);
      }
    }
  }

  // Set flag if all (rank-local) rates are zero so that Reaction transform can skip source term construction
  result.is_zero = std::all_of(region.begin(), region.end(),
                               [&](const auto& i) { return result.rate[i] == 0.0; });
  return result;
}

/**
 * @brief Get the rate for this pseudo-reaction based on the current simulation state.
 *
 * @param state GuardedOptions object containing the current state
 * @return RateData Calculated rate data for the pseudo-reaction
 */
inline RateData PseudoReaction::get_rate(GuardedOptions& state) {
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
 * Since the only reactant has a +ve population change, Reaction::transform_impl won't generate momentum
 * or energy sources. Add them here instead, unless they've been explicitly disabled.
 *
 */
inline void PseudoReaction::transform_additional(GuardedOptions& state,
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

#endif // PSEUDO_REACTION_H
