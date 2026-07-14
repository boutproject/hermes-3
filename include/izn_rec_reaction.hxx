#pragma once
#ifndef IZN_REC_REACTION_H
#define IZN_REC_REACTION_H

#include "reaction.hxx"
#include <string>

namespace hermes {

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
      : IznReaction(name, options) {};
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
      : RecReaction(name, options) {};
};

} // namespace hermes

namespace {

/// Register components for Hydrogen isotope ionisation and recombination
struct IznH : public hermes::IznReaction {
  using IznReaction::IznReaction;
  static constexpr auto type = "h + e -> h+ + 2e";
};
struct IznD : public hermes::IznReaction {
  using IznReaction::IznReaction;
  static constexpr auto type = "d + e -> d+ + 2e";
};
struct IznT : public hermes::IznReaction {
  using IznReaction::IznReaction;
  static constexpr auto type = "t + e -> t+ + 2e";
};
struct RecH : public hermes::RecReaction {
  using RecReaction::RecReaction;
  static constexpr auto type = "h+ + e -> h";
};
struct RecD : public hermes::RecReaction {
  using RecReaction::RecReaction;
  static constexpr auto type = "d+ + e -> d";
};
struct RecT : public hermes::RecReaction {
  using RecReaction::RecReaction;
  static constexpr auto type = "t+ + e -> t";
};
RegisterComponent<IznH> register_izn_h;
RegisterComponent<IznD> register_izn_d;
RegisterComponent<IznT> register_izn_t;
RegisterComponent<RecH> register_rec_h;
RegisterComponent<RecD> register_rec_d;
RegisterComponent<RecT> register_rec_t;

/// Register components for Helium ionisation and recombination
struct IznHe : public hermes::IznReaction {
  using IznReaction::IznReaction;
  static constexpr auto type = "he + e -> he+ + 2e";
};
struct RecHe : public hermes::RecReaction {
  using RecReaction::RecReaction;
  static constexpr auto type = "he+ + e -> he";
};
RegisterComponent<IznHe> register_izn_he;
RegisterComponent<RecHe> register_rec_he;

/*
 He+ ionisation (Amjuel data would be 2.2C, page 189) is not included yet
 Currently missing energy loss / radiation data
*/
/// Register components for Helium ionisation and recombination
// struct IznHe2 : public hermes::IznReaction {
//   using IznReaction::IznReaction;
//   static constexpr auto type = "e + h2+ -> he+2 + 2e";
// };
// RegisterComponent<IxnHe2> register_izn_hep;

/*
 He+2 recombination is not included yet
*/
// struct RecHe2 : public hermes::RecReaction {
//   using RecReaction::RecReaction;
//   static constexpr auto type = "he+2 + e -> he+";
// };
// RegisterComponent<RecHe2>register_rec_hep2;

} // namespace

#endif // IZN_REC_REACTION_H
