#pragma once
#ifndef REACTION_H
#define REACTION_H

#include "component.hxx"
#include "reaction_parser.hxx"
struct Reaction : public Component {
  Reaction(std::string name, Options& alloptions);

  static int get_instance_num() {
    static int instance_num{0};
    return instance_num++;
  }
  void transform(Options& state) override final;

protected:
  // Reaction string parser
  std::unique_ptr<ReactionParser> parser;

  /// Normalisations
  BoutReal Tnorm, Nnorm, FreqNorm;

  // Stoichiometric table (species_name -> population_change)
  std::map<std::string, int> S;

  // For diagnostics
  bool diagnose;                                  ///< Outputting diagnostics?
  BoutReal rate_multiplier, radiation_multiplier; ///< Scaling factor on reaction rate
  //   Field3D S;                                      ///< Particle exchange
  //   Field3D F;                                      ///< Momentum exchange
  //   Field3D E;                                      ///< Energy exchange
  //   Field3D R;                                      ///< Radiation loss

  /**
   * @brief Evaluates electron energy loss rate coefficients at a particular density and
   * temperature (Subclasses MUST define)
   *
   * @param T a temperature
   * @param n a density
   * @return BoutReal the electron energy loss rate
   */
  virtual BoutReal eval_electron_energy_loss_rate(BoutReal T, BoutReal n) = 0;

  //
  //
  /**
   * @brief Evaluate reaction rate coefficients at a particular density and temperature
   * (Subclasses MUST define)
   *
   * @param T a temperature
   * @param n a density
   * @return BoutReal the reaction rate
   */
  virtual BoutReal eval_reaction_rate(BoutReal T, BoutReal n) = 0;

  /**
   * @brief Set the diagnostic fields object (Subclasses MAY define)
   *
   * @param reaction_rate
   * @param momentum_exchange
   * @param energy_exchange
   * @param energy_loss
   */
  virtual void set_diagnostic_fields(Field3D& reaction_rate, Field3D& momentum_exchange,
                                     Field3D& energy_exchange, Field3D& energy_loss){};

  /**
   * @brief A hook with with subclasses can perform additional transform tasks, if
   * necessary. (Subclasses MAY define)
   *
   * @param state
   * @param reaction_rate
   * @param momentum_exchange
   * @param energy_exchange
   * @param energy_loss
   */
  virtual void transform_additional(Options& state, Field3D& reaction_rate,
                                    Field3D& momentum_exchange, Field3D& energy_exchange,
                                    Field3D& energy_loss) {}

private:
  /// Label to use for this reaction in a state / Options object
  const std::string name;
};
#endif