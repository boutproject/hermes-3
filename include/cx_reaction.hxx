#pragma once
#ifndef CX_REACTION_H
#define CX_REACTION_H

#include "reaction.hxx"
#include <string>

namespace hermes {

/**
 * @brief `Reaction` subclass that handles charge exchange reactions.
 *
 * @details  See `CXReaction::transform_additional` for the CX-specific sources.
 */
struct CXReaction : public Reaction {
  /**
   * @brief Main constructor for CXReaction.
   *
   * @details Adds appropriate permissions and diagnostics, checks that reaction string is
   * a valid CX reaction.
   *
   * @param name
   * @param alloptions  The options object
   */
  CXReaction(std::string name, Options& alloptions);

  /**
   * @brief CXReaction constructor used by component factory.
   *
   * @param name
   * @param alloptions  The options object
   * @param solver  The solver object for the simulation (discarded by this class)
   */
  CXReaction(std::string name, Options& alloptions, Solver*);

  /**
   * @brief Characterise the roles of the species involved in this reaction and check
   * that the species are valid for CX. Intended to be term-order-independent.
   *
   * @details Use the reaction parser to categorise reactants and products. See `r1`,
   * `r2`, `p1` and `p2` docstrings for details.
   */
  void set_species_and_validate();

  /**
   * @brief Perform additional transform tasks specific to CX reactions.
   *
   * @param state  The current simulation state to be modified
   * @param rate_data  The rate calculation results from parent class
   */
  void transform_additional(GuardedOptions& state, const RateData& rate_data) final;

private:
  /// True if one reactant has zero charge, false otherwise.
  bool has_neutral_reactant;
  /// Make neutral-ion CX behave as in diffusive neutrals
  bool no_neutral_cx_mom_gain;

  /// Reactant in the lower charge state
  ///   (e.g. 'h' in 'ne+3 + h -> ne+2 + h+' or 'd' in 'd + t+ -> d+ + t')
  std::string r1;
  /// Reactant in the higher charge state
  ///   (e.g. 'ne+3' in 'ne+3 + h -> ne+2 + h+' or 't+' in 'd + t+ -> d+ + t')
  std::string r2;
  /// Product formed by r1 losing an electron
  ///   (e.g. 'h+' in 'ne+3 + h -> ne+2 + h+' or 'd+' in 'd + t+ -> d+ + t')
  std::string p1;
  /// Product formed by r2 gaining an electron
  ///   (e.g. 'ne+2' in 'ne+3 + h -> ne+2 + h+' or 't' in 'd + t+ -> d+ + t')
  std::string p2;
};

} // namespace hermes
namespace {
/// Register components for symmetric HCX, one per isotope
RegisterComponent<hermes::CXReaction> register_cx_hh("h + h+ -> h+ + h");
RegisterComponent<hermes::CXReaction> register_cx_dd("d + d+ -> d+ + d");
RegisterComponent<hermes::CXReaction> register_cx_tt("t + t+ -> t+ + t");

/// Register components for non-symmetric HCX; one per atom,ion combination
RegisterComponent<hermes::CXReaction> register_cx_hd("h + d+ -> h+ + d");
RegisterComponent<hermes::CXReaction> register_cx_dh("d + h+ -> d+ + h");

RegisterComponent<hermes::CXReaction> register_cx_ht("h + t+ -> h+ + t");
RegisterComponent<hermes::CXReaction> register_cx_th("t + h+ -> t+ + h");

RegisterComponent<hermes::CXReaction> register_cx_dt("d + t+ -> d+ + t");
RegisterComponent<hermes::CXReaction> register_cx_td("t + d+ -> t+ + d");
} // namespace

#endif // CX_REACTION_H
