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

  // Stoich vector
  std::map<std::string, int> S;

  //   void transform(Options& state) override {
  //     Options& electron = state["species"]["e"];
  //     Options& from = state["species"][this->from_species]; // e.g. "h"
  //     Options& to = state["species"][this->to_species];     // e.g. "h+"
  //     Field3D reaction_rate, momentum_exchange, energy_exchange, energy_loss;

  //     electron_reaction(electron, from, to, get_rate_coeffs(), get_rad_coeffs(),
  //                       get_electron_heating(), reaction_rate, momentum_exchange,
  //                       energy_exchange, energy_loss, this->rate_multiplier,
  //                       this->radiation_multiplier);

  //     if (this->diagnose) {
  //       set_diagnostic_fields(reaction_rate, momentum_exchange, energy_exchange,
  //                             energy_loss);
  //     }
  //   }

  /**
   * @brief Set values in Field3D objects associated with diagnostics.
   *
   * @todo Make pure virtual to force implementation?
   *
   * @param reaction_rate
   * @param momentum_exchange
   * @param energy_exchange
   * @param energy_loss
   */

protected:
  /// Normalisations
  BoutReal Tnorm, Nnorm, FreqNorm;

  // For diagnostics
  bool diagnose;                                  ///< Outputting diagnostics?
  BoutReal rate_multiplier, radiation_multiplier; ///< Scaling factor on reaction rate
  //   Field3D S;                                      ///< Particle exchange
  //   Field3D F;                                      ///< Momentum exchange
  //   Field3D E;                                      ///< Energy exchange
  //   Field3D R;                                      ///< Radiation loss

  virtual void set_diagnostic_fields(Field3D& reaction_rate, Field3D& momentum_exchange,
                                     Field3D& energy_exchange, Field3D& energy_loss){};

private:
  const std::string name;
  std::unique_ptr<ReactionParser> parser;
};
#endif