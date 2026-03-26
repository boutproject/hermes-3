#pragma once
#ifndef IZN_REC_REACTION_H
#define IZN_REC_REACTION_H

#include "reaction.hxx"
#include <string>

/**
 * @brief `Reaction` subclass that acts as a common base for ionisation and
 * recombination reactions.
 *
 * @details See `IznRecReaction::transform_additional` for the
 * ionisation/recombination-specific sources.
 */
struct IznRecReaction : public Reaction {
protected:
  /**
   * @brief Main constructor for ionisation/recombination reaction base class.
   *
   * @details Adds appropriate permissions and diagnostics for ionisation/recombination
   * reactions, checks that the reaction string includes at least one heavy reactant and
   * one heavy product.
   *
   * @param name
   * @param options  The options object containing configuration parameters
   */
  IznRecReaction(std::string short_reaction_type, std::string name, Options& options);

  /**
   * @brief Perform additional transform tasks specific to ionisation/recombination
   * reactions.
   *
   * @details This function adds:
   * - Sources to handle the kinetic=>thermal energy transfer for electrons, and for the
   * heavy product.
   * - An electron energy source that capture the losses associated with the ionisation
   * potential energy cost as well as the photon emission during excitation and
   * de-excitation.
   *
   * @param state  The current simulation state to be modified
   * @param rate_data  The rate calculation results from parent class
   */
  void transform_additional(GuardedOptions& state, const RateData& rate_data) final;

  /// Name of the heavy reactant species
  std::string heavy_reactant;
  /// Name of the heavy product species
  std::string heavy_product;
  /// Name of the (heavy) species with which collision freqs. are associated in the state
  std::string heavy_collfreq_species;

private:
  /// Short reaction type string used in diagnostic names ("iz" or "rec")
  const std::string short_reaction_type;

  /// Pointer to reaction data used to compute electron energy loss
  std::unique_ptr<ReactionData> e_energy_loss_data;
};

/**
 * @brief `Reaction` subclass for ionisation reactions.
 *
 * @details This class only exists to add ionisation-specific options; see
 * `IznRecReaction::transform_additional` for related sources.
 */
struct IznReaction : public IznRecReaction {
  /**
   * @brief Main constructor for ionisation reaction objects.
   *
   * @param name
   * @param options The options object
   */
  IznReaction(std::string name, Options& options);
  /**
   * @brief Constructor used by the component factory.
   *
   * @param name
   * @param options  The options object
   * @param solver  The solver object for the simulation (discarded by this class)
   */
  IznReaction(std::string name, Options& options, [[maybe_unused]] Solver* solver)
      : IznReaction(name, options){};
};

/**
 * @brief `Reaction` subclass for recombination reactions.
 *
 * @details This class only exists to add recombination-specific options; see
 * `IznRecReaction::transform_additional` for related sources.
 */
struct RecReaction : public IznRecReaction {
  /**
   * @brief Main constructor for recombination reaction objects.
   *
   * @param name
   * @param options
   */
  RecReaction(std::string name, Options& options);
  /**
   * @brief Constructor used by the component factory.
   *
   * @param name
   * @param options  The options object
   * @param solver  The solver object for the simulation (discarded by this class)
   */
  RecReaction(std::string name, Options& options, [[maybe_unused]] Solver* solver)
      : RecReaction(name, options){};
};

namespace {
/// Register components for Hydrogen isotope ionisation and recombination
RegisterComponent<IznReaction> register_izn_h("h + e -> h+ + 2e");
RegisterComponent<IznReaction> register_izn_d("d + e -> d+ + 2e");
RegisterComponent<IznReaction> register_izn_t("t + e -> t+ + 2e");
RegisterComponent<RecReaction> register_rec_h("h+ + e -> h");
RegisterComponent<RecReaction> register_rec_d("d+ + e -> d");
RegisterComponent<RecReaction> register_rec_t("t+ + e -> t");

/// Register components for Helium ionisation and recombination
RegisterComponent<IznReaction> register_izn_he("he + e -> he+ + 2e");
RegisterComponent<RecReaction> register_rec_he("he+ + e -> he");

/*
 He+ ionisation (Amjuel data would be 2.2C, page 189) is not included yet
 Currently missing energy loss / radiation data
*/
// RegisterComponent<IznReaction> register_izn_hep("e + he+ -> he+2 + 2e");

/*
 He+2 recombination is not included yet
*/
// RegisterComponent<RecReaction>register_rec_hep2("he+2 + e -> he+");

} // namespace

#endif // IZN_REC_REACTION_H
