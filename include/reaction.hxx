#pragma once
#ifndef REACTION_H
#define REACTION_H

#include <map>
#include <memory>

#include "component.hxx"
#include "rate_helper.hxx"
#include "reaction_data.hxx"
#include "reaction_diagnostic.hxx"
#include "reaction_parser.hxx"

namespace hermes {

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
      : Component(std::move(permissions)), inst_num(incremented_instance_num() + 1) {}

  static std::size_t incremented_instance_num() { return instance_num_ref()++; }

  // Reset the instance counter; needed to avoid unit tests affecting each other!
  static void reset_instance_counter() { instance_num_ref() = 0; }

private:
  static std::size_t& instance_num_ref() {
    static std::size_t instance_num{0};
    return instance_num;
  }

protected:
  std::size_t inst_num;
};

/**
 * @brief Struct intended to act as a base for all reactions.
 *
 * @details Stores a ReactionParser for manipulating the reaction string and reads data
 * via (a subclass of) ReactionData. Also computes generic population change source terms
 * in `transfrom_impl` and records reaction 'channels' that can be used to override the
 * default strategy for distributing momentum and energy between reaction products.
 *
 */
struct Reaction : public ReactionBase {
  /**
   * @brief Construct a new Reaction object
   *
   * @details Extract reaction data options, parse the reaction string, set some
   * permissions and options that apply to all reactions.
   *
   * @param name
   * @param options Options object
   */
  Reaction(std::string name, Options& options);

  /**
   * @brief Copy all diagnostics into the output, setting the appropriate metadata at the
   * same time. Subclasses can't override - instead, they should make add_diagnostic()
   * calls in their constructor.
   *
   * @param state the output state to update
   */
  void outputVars(Options& state) final;

protected:
  /// Reaction string parser
  std::unique_ptr<ReactionParser> parser;

  /// Reaction data
  std::unique_ptr<ReactionData> rate_data;

  /// Units and normalisations extracted to member vars for convenience
  Options& units;
  BoutReal Tnorm, Nnorm, FreqNorm;

  /// Rate multiplier. Defaults to 1, but subclasses may read an option and overwrite.
  BoutReal rate_multiplier;

  /// Output diagnostics?
  bool diagnose;

  /// map of (species_name,diagnostic_type)->diagnostic_object
  std::multimap<std::pair<std::string, ReactionDiagnosticType>, ReactionDiagnostic>
      diagnostics;

  /**
   * Whether to use parallel averaging when calculating reaction rates
   * (Defaults to true)
   */
  bool do_parallel_averaging = true;

  /**
   * @brief Add a new reaction diagnostic.
   *
   * @param sp_name Species with which the diagnostic will be associated
   * @param diag_name Label used in the output (and to store it temporarily in the state)
   * @param diag_desc Description to use as the 'long_name' output attribute
   * @param diag_type enum identifying the diagnostic type, also used to determine source
   * name
   * @param data_source Name to use as the 'source' output attribute
   * @param standard_name Optional string to use as the 'standard_name' output attribute.
   * Defaults to diag_name.
   * @param transformer Optional transformer function to use when modifying the diagnostic
   * (default is 'negate', i.e. the diagnostic has the opposite sign to the source)
   *
   * @details Adds a new entry in the diagnostic (multi)map. The (non-unique)
   * key is < \p sp_name, \p type >
   */
  void add_diagnostic(const std::string& sp_name, const std::string& diag_name,
                      const std::string& long_diag_name, ReactionDiagnosticType type,
                      const std::string& data_source,
                      DiagnosticTransformerType transformer = negate,
                      const std::string& standard_name = "");

  /**
   * @brief Set weights for any reactant => product momentum / energy channel that hasn't
   * already been specified via set_energy_channel_weight and set_momentum_channel_weight.
   *
   * @note Can't be done at construction because the species masses may not be set.
   *
   * @param state Current sim state
   */
  void init_channel_weights(GuardedOptions& state);

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
                     ReactionDiagnosticType type, const Field3D& update_with_field) {

    update_source<operation>(state, sp_name, type, sp_data_keys.at(type),
                             update_with_field);
  }

  /**
   * @brief Update both a species source term and the corresponding diagnostics (if any
   * exist and if diagnostics are enabled)
   *
   * @tparam operation function to call on the state to update the source term and
   * the diagnostic. Either Component::add, Component::subtract or Component::set
   * @param state the state to update
   * @param sp_name the name of the species to update
   * @param type the type of source/diagnostic to update
   * @param sp_data_key label/key for the field in the state object, i.e.
   * state["species"][sp_name][sp_data_key]
   * @param update_with_field the field used in the update
   */
  template <OPTYPE operation>
  void update_source(GuardedOptions& state, const std::string& sp_name,
                     ReactionDiagnosticType type, const std::string& sp_data_key,
                     const Field3D& update_with_field) {
    // Update species data
    operation(state["species"][sp_name][sp_data_key], update_with_field);

    if (this->diagnose) {
      // Update corresponding diagnostic(s) (if any exist)
      auto matches = this->diagnostics.equal_range(std::make_pair(sp_name, type));
      for (auto match = matches.first; match != matches.second; match++) {
        Field3D diag_src_fld = match->second.transform(update_with_field);
        // Apply the update to the diagnostic field in the state, then copy it to the
        // diagnostic
        operation(state[match->second.get_name()], diag_src_fld);
        match->second.set_data(getNonFinal<Field3D>(state[match->second.get_name()]));
      }
    }
  }

private:
  // Channels to determine how momentum and energy are distributed to product species
  std::map<std::string, std::map<std::string, BoutReal>> energy_channels;
  std::map<std::string, std::map<std::string, BoutReal>> momentum_channels;

  /// Label used in the state for reaction configuration.
  const std::string name;

  /// Participation factors of all species - currently set to 1!
  std::map<std::string, BoutReal> pfactors;

  /**
   * @brief Extract reaction string and data type and data ID for this reaction from the
   * input options, or set suitable defaults for the type and ID if they aren't specified.
   *
   * @param options
   * @param name
   * @param[out] reaction_str the extracted reaction string
   * @param[out] data_type the extracted data type enum or a default if no type specified
   * @param[out] data_id the extracted data id or a default if no id specified
   *
   * @details The current input file format (all reactions in a single, comma-separated
   string) is a bit awkward, but is being preserved for now. ReactionBase sets
   this->inst_num according to the order of instantiation for each reaction object, then
   this function extracts the reaction string, data type and data ID using inst_num as an
   index, setting suitable type and ID defaults if either are omitted.
   */
  void get_reaction_settings(Options& options, std::string& reaction_str,
                             ReactionDataTypes& data_type, std::string& data_id);

  /**
   * @brief Add density, momentum and energy sources that apply to all reactions (e.g.
   * those driven by species population changes), then call transform_additional() to
   * allow subclasses to add other terms.
   *
   * @param state
   */
  void transform_impl(GuardedOptions& state) override final;

  /**
   * @brief Reset the temporary values of the diagnostics stored in the state.
   *
   * @param state
   */
  void zero_diagnostics(GuardedOptions& state);
};

} // namespace hermes

#endif
