#pragma once
#ifndef MOLECULAR_REACTIONS_H
#define MOLECULAR_REACTIONS_H

#include "cx_reaction.hxx"
#include "reaction.hxx"
#include <string>

namespace hermes {

/**
 * @brief Base class for dissociation reactions involving molecules.
 *
 * @details Inherits from Reaction to handle dissociation of molecular species.
 */
struct Dissociation : public Reaction {
  /**
   * @brief Main constructor for Dissociation.
   *
   * @param name
   * @param options The options object
   */
  Dissociation(std::string name, Options& options) : Reaction(name, options) {};

  /**
   * @brief Constructor used by component factory.
   *
   * @param name
   * @param options The options object
   * @param solver The solver object for the simulation
   */
  Dissociation(std::string name, Options& options, [[maybe_unused]] Solver* solver)
      : Dissociation(name, options) {};
};

/**
 * @brief Base class for dissociative excitation reactions involving molecules.
 */
struct DissociativeExc : public Reaction {
  /**
   * @brief Main constructor for DissociativeExc.
   *
   * @param name
   * @param options The options object
   */
  DissociativeExc(std::string name, Options& options) : Reaction(name, options) {};

  /**
   * @brief Constructor used by component factory.
   *
   * @param name
   * @param options The options object
   * @param solver The solver object for the simulation
   */
  DissociativeExc(std::string name, Options& options, [[maybe_unused]] Solver* solver)
      : DissociativeExc(name, options) {};
};

/**
 * @brief Base class for dissociative ionization reactions involving molecules.
 */
struct DissociativeIzn : public Reaction {
  /**
   * @brief Main constructor for DissociativeIzn.
   *
   * @param name
   * @param options The options object
   */
  DissociativeIzn(std::string name, Options& options) : Reaction(name, options) {};

  /**
   * @brief Constructor used by component factory.
   *
   * @param name
   * @param options The options object
   * @param solver The solver object for the simulation
   */
  DissociativeIzn(std::string name, Options& options, [[maybe_unused]] Solver* solver)
      : DissociativeIzn(name, options) {};
};

/**
 * @brief Base class for non-dissociative ionization reactions involving molecules.
 */
struct NonDissociativeIzn : public Reaction {
  /**
   * @brief Main constructor for NonDissociativeIzn.
   *
   * @param name
   * @param options The options object
   */
  NonDissociativeIzn(std::string name, Options& options) : Reaction(name, options) {};

  /**
   * @brief Constructor used by component factory.
   *
   * @param name
   * @param options The options object
   * @param solver The solver object for the simulation
   */
  NonDissociativeIzn(std::string name, Options& options, [[maybe_unused]] Solver* solver)
      : NonDissociativeIzn(name, options) {};
};

/**
 * @brief Base class for dissociative recombination reactions involving molecules.
 */
struct DissociativeRec : public Reaction {
  /**
   * @brief Main constructor for DissociativeRec.
   *
   * @param name
   * @param options The options object
   */
  DissociativeRec(std::string name, Options& options) : Reaction(name, options) {};

  /**
   * @brief Constructor used by component factory.
   *
   * @param name
   * @param options The options object
   * @param solver The solver object for the simulation
   */
  DissociativeRec(std::string name, Options& options, [[maybe_unused]] Solver* solver)
      : DissociativeRec(name, options) {};
};

} // namespace hermes

#endif // MOLECULAR_REACTIONS_H
