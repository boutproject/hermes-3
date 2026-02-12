#pragma once
#ifndef REACTION_H
#define REACTION_H

#include "component.hxx"
#include "rate_helper.hxx"
#include "reaction_diagnostic.hxx"
#include "reaction_parser.hxx"

using OPTYPE = GuardedOptions && (GuardedOptions&&, Field3D);

/**
 * @brief Temporary struct to use as a base class for all reactions components.
 *
 * @details Ensures reaction strings are paired up correctly with component classes.
 *
 * @todo Remove this when all reaction classes have been refactored to inherit from
 * Reaction.
 */
struct ReactionBase : public Component {
  ReactionBase(Permissions&& permissions)
      : Component(std::move(permissions)), inst_num(get_instance_num() + 1) {}
  static int get_instance_num() {
    static int instance_num{0};
    return instance_num++;
  }

protected:
  int inst_num;
};

/**
 * @brief Struct intended to act as a base for all reactions.
 *
 * @details Note that several parameterisations are possible when providing reaction rate
 * data; subclasses are expected to (re-)implement one or more of eval_sigma_v_ET,
 * eval_sigma_v_nT and eval_sigma_v_T. Subclasses that have an associated electron energy
 * loss rate are expected to implement eval_sigma_vE_nT. Subclasses that DO NOT have an
 * associated electron energy loss rate should set includes_sigma_v_e=false.
 *
 */
struct Reaction : public ReactionBase {
  Reaction(std::string name, Options& alloptions);

  void outputVars(Options& state) final;

protected:
  /// Reaction string parser
  std::unique_ptr<ReactionParser> parser;

  /// Units and normalisations extracted to member vars for convenience
  Options& units;
  BoutReal Tnorm, Nnorm, FreqNorm;

  /// Rate multipliers, extracted from input options
  BoutReal rate_multiplier, radiation_multiplier;

  /// Output diagnostics?
  bool diagnose;

  /// map of (species_name,diagnostic_type)->diagnostic_object
  std::multimap<std::pair<std::string, ReactionDiagnosticType>, ReactionDiagnostic>
      diagnostics;

  /// Whether to use parallel averaging when calculating reaction rates
  /// (Default to true)
  bool do_parallel_averaging = true;

  /// Whether or not reaction data includes <sigma v E>
  /// (Default to true as a reminder to override eval_sigma_vE_nT)
  bool includes_sigma_v_e = true;

  /**
   * @brief Get the type of parameters expected by the rate calculation.
   * @details Allow the rate parameter types to be provided after construction
   * (e.g. after subclasses have read that information from file) by forcing subclasses to
   * implement this function, rather than storing the types as a Reaction member variable.
   * @return RateParamsTypes the rate params type
   */
  virtual RateParamsTypes get_rate_params_type() const = 0;

  /**
   * @brief Add a new entry in this Reaction's diagnostic (multi)map. The (non-unique) Key
   * is < \p sp_name, \p type >
   *
   * @param sp_name name of the species with which to associate the diagnostic
   * @param diag_name name of the diagnostic (also the key used when updating it in the
   * state object)
   * @param long_diag_name doc string to use as the diagnostic description
   * @param type an enum that associates the diagnostic with density, momentum or energy
   * sources
   * @param data_source name of the associated data source (e.g. 'Amjuel H.x.y')
   * @param transformer optional transformer function to call on field data when the
   * diagnostic is updated
   * @param standard_name optional 'standard_name' to use in the output file
   */
  void add_diagnostic(const std::string& sp_name, const std::string& diag_name,
                      const std::string& long_diag_name, ReactionDiagnosticType type,
                      const std::string& data_source,
                      DiagnosticTransformerType transformer = negate,
                      const std::string& standard_name = "");

  /**
   * @brief Calculate weightsums used in transform(). Can't be done at construction
   * because the species masses may not be set.
   *
   *
   * @param state Current sim state
   */
  void init_channel_weights(GuardedOptions& state);

  /**
   * @brief Evaluate <sigma . v . E> at a particular density and temperature
   * (Subclasses MAY define)
   *
   * @param T a temperature
   * @param n a density
   * @return BoutReal the electron energy loss rate
   */
  virtual BoutReal eval_sigma_vE_nT([[maybe_unused]] BoutReal T,
                                    [[maybe_unused]] BoutReal n) {
    if (this->includes_sigma_v_e) {
      throw BoutException(
          "eval_sigma_vE_nT() needs to be implemented by Reaction instances "
          "which set includes_sigma_v_e=true");
    } else {
      throw BoutException(
          "eval_sigma_vE_nT() was called despite having set includes_sigma_v_e=false!");
    }
    return -1;
  };

  /**
   * @brief Evaluate <sigma.v> at a particular energy and temperature. Rate args are only
   * known at runtime, so need to provide a default implementation that throws and let
   * subclasses decide whether to override.
   *
   * @param E a density
   * @param T a temperature
   * @return BoutReal <sigma.v>(E,T)
   */
  virtual BoutReal eval_sigma_v_ET([[maybe_unused]] BoutReal E,
                                   [[maybe_unused]] BoutReal T) {
    throw BoutException("Trying to call eval_sigma_v_nT but the Reaction subclass hasn't "
                        "implemented it!");
  }
  /**
   * @brief Evaluate <sigma.v> at a particular density and temperature. Rate args are only
   * known at runtime, so need to provide a default implementation that throws and let
   * subclasses decide whether to override.
   *
   * @param T a temperature
   * @return BoutReal <sigma.v>(T)
   */
  virtual BoutReal eval_sigma_v_T([[maybe_unused]] BoutReal T) {
    throw BoutException("Trying to call eval_sigma_v_T but the Reaction subclass hasn't "
                        "implemented it!");
  }

  /**
   * @brief Evaluate <sigma.v> at a particular density and temperature. Rate args are only
   * known at runtime, so need to provide a default implementation that throws and let
   * subclasses decide whether to override.
   *
   * @param T a temperature
   * @param n a density
   * @return BoutReal <sigma.v>(n,T)
   */
  virtual BoutReal eval_sigma_v_nT([[maybe_unused]] BoutReal T,
                                   [[maybe_unused]] BoutReal n) {
    throw BoutException("Trying to call eval_sigma_v_nT but the Reaction subclass hasn't "
                        "implemented it!");
  }

  /**
   * @brief A hook with which subclasses can perform additional transform tasks, over and
   * above those implemented in Reaction::transform. (Subclasses MAY define)
   *
   * @param state
   * @param rate_calc_results
   */
  virtual void transform_additional([[maybe_unused]] GuardedOptions& state,
                                    [[maybe_unused]] const RateData& rate_calc_results) {}

  /**
   * @brief Update both a species source term and the corresponding diagnostics (if any
   * exist and if diagnostics are enabled), determining the key in the state from the
   * type. See alternative form of update_source for further details.
   */
  template <OPTYPE operation>
  void update_source(GuardedOptions& state, const std::string& sp_name,
                     ReactionDiagnosticType type, const Field3D& fld) {

    update_source<operation>(state, sp_name, type, state_labels.at(type), fld);
  }

  /**
   * @brief Update both a species source term and the corresponding diagnostics (if any
   * exist and if diagnostics are enabled)
   *
   * @tparam operation function to call on the state to update the source term and
   * the diagnostic. Either Component::add, Component::subtract or Component::set
   * @param state the state to update
   * @param sp_name the species to update
   * @param type the type of source/diagnostic to update
   * @param lbl label/key for the field in the state object, i.e.
   * state["species"][sp_name][lbl]
   * @param fld the field used in the update
   */
  template <OPTYPE operation>
  void update_source(GuardedOptions& state, const std::string& sp_name,
                     ReactionDiagnosticType type, const std::string& lbl, const Field3D& fld) {
    // Update species data
    operation(state["species"][sp_name][lbl], fld);

    if (this->diagnose) {
      // Update corresponding diagnostic(s) (if any exist)
      auto matches = this->diagnostics.equal_range(std::make_pair(sp_name, type));
      for (auto match = matches.first; match != matches.second; match++) {
        Field3D diag_src_fld = match->second.transform(fld);
        // Apply the update to the diagnostic field in the state, then copy it to the
        // diagnostic
        operation(state[match->second.get_name()], diag_src_fld);
        match->second.set_data(getNonFinal<Field3D>(state[match->second.get_name()]));
      }
    }
  }

  /**
   * @brief Specify what fraction of a reactant's energy is transferred to a particular
   * product.
   *
   * @param reactant_name Name of the reactant species.
   * @param product_name Name of the product species.
   * @param weight Fraction of the energy to transfer.
   */
  void set_energy_channel_weight(const std::string& reactant_name,
                                 const std::string& product_name, BoutReal weight);

  /**
   * @brief Specify what fraction of a reactant's momentum is transferred to a particular
   * product.
   *
   * @param reactant_name Name of the reactant species.
   * @param product_name Name of the product species.
   * @param weight Fraction of the momentum to transfer.
   */
  void set_momentum_channel_weight(const std::string& reactant_name,
                                   const std::string& product_name, BoutReal weight);

private:
  // Channels to determine how momentum and energy are distributed to product species
  std::map<std::string, std::map<std::string, BoutReal>> energy_channels;
  std::map<std::string, std::map<std::string, BoutReal>> momentum_channels;

  /// Label to use for this reaction in a state / Options object
  const std::string name;

  /// Participation factors of all species
  std::map<std::string, BoutReal> pfactors;

  void zero_diagnostics(GuardedOptions& state);

  void transform_impl(GuardedOptions& state) override final;
};
#endif
