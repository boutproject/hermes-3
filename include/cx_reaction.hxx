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
// hermes::CXReaction couldn't inherit from NamedComponent while
// keeping the function implementations in the .cxx file. Instead
// we apply CRTP here.
template <typename T>
struct CXReaction : public hermes::CXReaction {
  using hermes::CXReaction::CXReaction;
  std::string typeName() const final { return T::type; }
};

/// Register components for symmetric HCX, one per isotope
struct CXHH : CXReaction<CXHH> {
  using CXReaction<CXHH>::CXReaction;
  static constexpr auto type = "h + h+ -> h+ + h";
};
struct CXDD : CXReaction<CXDD> {
  using CXReaction<CXDD>::CXReaction;
  static constexpr auto type = "d + d+ -> d+ + d";
};
struct CXTT : CXReaction<CXTT> {
  using CXReaction<CXTT>::CXReaction;
  static constexpr auto type = "t + t+ -> t+ + t";
};
RegisterComponent<CXHH> register_cx_hh;
RegisterComponent<CXDD> register_cx_dd;
RegisterComponent<CXTT> register_cx_tt;

/// Register components for non-symmetric HCX; one per atom,ion combination
struct CXHD : CXReaction<CXHD> {
  using CXReaction<CXHD>::CXReaction;
  static constexpr auto type = "h + d+ -> h+ + d";
};
struct CXDH : CXReaction<CXDH> {
  using CXReaction<CXDH>::CXReaction;
  static constexpr auto type = "d + h+ -> d+ + h";
};
RegisterComponent<CXHD> register_cx_hd;
RegisterComponent<CXDH> register_cx_dh;

struct CXHT : CXReaction<CXHT> {
  using CXReaction<CXHT>::CXReaction;
  static constexpr auto type = "h + t+ -> h+ + t";
};
struct CXTH : CXReaction<CXTH> {
  using CXReaction<CXTH>::CXReaction;
  static constexpr auto type = "t + h+ -> t+ + h";
};
RegisterComponent<CXHT> register_cx_ht;
RegisterComponent<CXTH> register_cx_th;

struct CXDT : CXReaction<CXDT> {
  using CXReaction<CXDT>::CXReaction;
  static constexpr auto type = "d + t+ -> d+ + t";
};
struct CXTD : CXReaction<CXTD> {
  using CXReaction<CXTD>::CXReaction;
  static constexpr auto type = "t + d+ -> t+ + d";
};
RegisterComponent<CXDT> register_cx_dt;
RegisterComponent<CXTD> register_cx_td;
} // namespace

#endif // CX_REACTION_H
