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

  /// Normalisations, extracted from input options
  BoutReal Tnorm, Nnorm, FreqNorm;

  // Rate multipliers, extracted from input options
  BoutReal rate_multiplier, radiation_multiplier;

  // Output diagnostics?
  bool diagnose;

  /**
   * @brief Calculate weightsums used in transform(). Can't be done at construction
   * because the species masses may not be set.
   *
   *
   * @param state Current sim state
   */
  void calc_weightsums(Options& state);

  /**
   * @brief Evaluates electron energy loss rate coefficients at a particular density and
   * temperature.
   * (Subclasses MUST define)
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
   * @brief Set the diagnostic fields object.
   * (Subclasses MAY define)
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
   * necessary.
   * (Subclasses MAY define)
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

  /// Participation factors of all species
  std::map<std::string, BoutReal> pfactors;

  /// Sum of weights to use when calculating energy source due to population change
  BoutReal energy_weightsum;
  /// Sum of weights to use when calculating momentum source due to population change
  BoutReal momentum_weightsum;
};

/**
 * @brief Struct to simplify cell-averaging of the reaction rate, particularly the
 * calculation of the mass action factor.
 *
 */
typedef std::function<BoutReal(BoutReal, BoutReal, BoutReal)> RateFunctionType;
template <typename LimiterType = hermes::Limiter, typename IdxType = Ind3D>
struct RateHelper {
  RateHelper(const Options& state, const std::vector<std::string>& reactant_species,
             RateFunctionType rate_calc_func, const Region<IdxType> region);

  Field3D calc_rate();

private:
  const Region<IdxType> region;
  /// Function to calculate reaction rate as a function of n_e, T_e
  RateFunctionType rate_calc_func;
  /// Electron density and temperature
  Field3D n_e;
  Field3D T_e;
  // Reactant densities
  std::vector<Field3D> n_reactants;

  BoutReal mass_action(IdxType i);

  BoutReal mass_action_left(IdxType i, IdxType ym, IdxType yp);

  BoutReal mass_action_right(IdxType i, IdxType ym, IdxType yp);
};
#endif