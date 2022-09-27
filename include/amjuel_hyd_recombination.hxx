#pragma once
#ifndef AMJUEL_HYD_RECOMBINATION_H
#define AMJUEL_HYD_RECOMBINATION_H

#include "amjuel_reaction.hxx"

/// Hydrogen recombination, Amjuel rates
///
/// Includes both radiative and 3-body recombination
struct AmjuelHydRecombination : public AmjuelReaction {
  AmjuelHydRecombination(std::string name, Options& alloptions, Solver* solver)
      : AmjuelReaction(name, alloptions, solver) {}

  void calculate_rates(Options& electron, Options& atom, Options& ion,
                        Field3D &reaction_rate, Field3D &momentum_exchange,
                        Field3D &energy_exchange, Field3D &energy_loss);
};

/// Hydrogen recombination
/// Templated on a char to allow 'h', 'd' and 't' species to be treated with the same code
template <char Isotope>
struct AmjuelHydRecombinationIsotope : public AmjuelHydRecombination {
  AmjuelHydRecombinationIsotope(std::string name, Options& alloptions, Solver* solver)
      : AmjuelHydRecombination(name, alloptions, solver) {

    diagnose = alloptions[name]["diagnose"]
      .doc("Output additional diagnostics?")
      .withDefault<bool>(false);

    if (diagnose) {
      // Save particle, momentum and energy channels

      bout::globals::dump.addRepeat(S, {'S', Isotope, '+', '_', 'r', 'e', 'c'}); // Particle source
      bout::globals::dump.addRepeat(F, {'F', Isotope, '+', '_', 'r', 'e', 'c'}); // Momentum exchange
      bout::globals::dump.addRepeat(E, {'E', Isotope, '+', '_', 'r', 'e', 'c'}); // Energy exchange
      bout::globals::dump.addRepeat(R, {'R', Isotope, '+', '_', 'r', 'e', 'c'}); // Radiation loss
      }
    
    
    }
  
  void transform(Options& state) override {
    Options& electron = state["species"]["e"];
    Options& atom = state["species"][{Isotope}];     // e.g. "h"
    Options& ion = state["species"][{Isotope, '+'}]; // e.g. "h+"
    Field3D reaction_rate, momentum_exchange, energy_exchange, energy_loss;

    calculate_rates(electron, atom, ion, reaction_rate, momentum_exchange, energy_exchange, energy_loss);

    if (diagnose) {
      S = -reaction_rate;
      F = -momentum_exchange;
      E = -energy_exchange;
      R = -energy_loss;
    }
  }
private:
  bool diagnose; ///< Outputting diagnostics?
  Field3D S; ///< Particle exchange
  Field3D F; ///< Momentum exchange
  Field3D E; ///< Energy exchange
  Field3D R; ///< Radiation loss
};

namespace {
/// Register three components, one for each hydrogen isotope
/// so no isotope dependence included.
RegisterComponent<AmjuelHydRecombinationIsotope<'h'>>
    register_recombination_h("h+ + e -> h");
RegisterComponent<AmjuelHydRecombinationIsotope<'d'>>
    register_recombination_d("d+ + e -> d");
RegisterComponent<AmjuelHydRecombinationIsotope<'t'>>
    register_recombination_t("t+ + e -> t");
} // namespace

#endif // AMJUEL_HYD_RECOMBINATION_H
