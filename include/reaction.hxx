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
  /// Normalisations
  BoutReal Tnorm, Nnorm, FreqNorm;

  // Stoich vector
  std::map<std::string, int> S;

  // For diagnostics
  bool diagnose;                                  ///< Outputting diagnostics?
  BoutReal rate_multiplier, radiation_multiplier; ///< Scaling factor on reaction rate
  //   Field3D S;                                      ///< Particle exchange
  //   Field3D F;                                      ///< Momentum exchange
  //   Field3D E;                                      ///< Energy exchange
  //   Field3D R;                                      ///< Radiation loss

  // Evaluate electron energy loss rate coefficients at a particular density and
  // temperature (delegated to subclasses)
  virtual BoutReal eval_radiation_rate(BoutReal T, BoutReal n) = 0;

  // Evaluate reaction rate coefficients at a particular density and temperature
  // (delegated to subclasses)
  virtual BoutReal eval_reaction_rate(BoutReal T, BoutReal n) = 0;

  virtual void set_diagnostic_fields(Field3D& reaction_rate, Field3D& momentum_exchange,
                                     Field3D& energy_exchange, Field3D& energy_loss){};

private:
  const std::string name;
  std::unique_ptr<ReactionParser> parser;
};
#endif