#pragma once
#ifndef ELASTIC_COLLISIONS_H
#define ELASTIC_COLLISIONS_H

#include "reaction.hxx"
#include <string>

namespace hermes {

/**
 * @brief Base class for elastic collision 'reactions'.
 */
struct ElasticCollision : public Reaction {
  /**
   * @brief Main constructor for ElasticCollision.
   *
   * @details Passes `add_pop_change_sources`=`false` to the Reaction
   * constructor to disable standard density, momentum, energy sources due to population
   * changes.
   *
   * @param name
   * @param options  The options object
   */
  ElasticCollision(std::string name, Options& options);

  /**
   * @brief ElasticCollision constructor used by component factory.
   *
   * @param name
   * @param options  The options object
   * @param solver  The solver object for the simulation (discarded by this class)
   */
  ElasticCollision(std::string name, Options& options, Solver*);

  /**
   * @brief Perform additional transformations specific to elastic collisions.
   *
   * @param state
   * @param rate_calc_results
   */
  void transform_additional([[maybe_unused]] GuardedOptions& state,
                            [[maybe_unused]] const RateData& rate_calc_results) final;

private:
  /// @brief Names of the two reactant species (== the product species)
  std::string r1, r2;
};

} // namespace hermes

#endif // ELASTIC_COLLISIONS_H
