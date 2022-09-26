#pragma once
#ifndef AMJUEL_HYD_IONISATION_H
#define AMJUEL_HYD_IONISATION_H

#include "amjuel_reaction.hxx"

/// Hydrogen ionisation, Amjuel rates
struct AmjuelHydIonisation : public AmjuelReaction {
  AmjuelHydIonisation(std::string name, Options& alloptions, Solver* solver)
      : AmjuelReaction(name, alloptions, solver) {}

  void calculate_rates(Options& electron, Options& atom, Options& ion, 
                        Field3D &reaction_rate, Field3D &momentum_exchange,
                        Field3D &energy_exchange, Field3D &energy_loss);
};

/// Hydrogen ionisation
/// Templated on a char to allow 'h', 'd' and 't' species to be treated with the same code
template <char Isotope>
struct AmjuelHydIonisationIsotope : public AmjuelHydIonisation {
  AmjuelHydIonisationIsotope(std::string name, Options& alloptions, Solver* solver)
      : AmjuelHydIonisation(name, alloptions, solver) {

    diagnose = alloptions[name]["diagnose"]
      .doc("Output additional diagnostics?")
      .withDefault<bool>(false);

    if (diagnose) {
      // Save particle, momentum and energy channels

      bout::globals::dump.addRepeat(S, {"Sd+_iz"}); // Particle source
      bout::globals::dump.addRepeat(F, {"Fd+_iz"}); // Momentum exchange
      bout::globals::dump.addRepeat(E, {"Ed+_iz"}); // Energy exchange
      bout::globals::dump.addRepeat(R, {"Rd+_ex"}); // Radiation loss
    }


  }

  void transform(Options& state) override {
    Options& electron = state["species"]["e"];
    Options& atom = state["species"][{Isotope}];     // e.g. "h"
    Options& ion = state["species"][{Isotope, '+'}]; // e.g. "h+"
    Field3D reaction_rate, momentum_exchange, energy_exchange, energy_loss;

    calculate_rates(electron, atom, ion, reaction_rate, momentum_exchange, energy_exchange, energy_loss);

    if (diagnose) {
      S = reaction_rate;
      F = momentum_exchange;
      E = energy_exchange;
      R = energy_loss;
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
RegisterComponent<AmjuelHydIonisationIsotope<'h'>>
    registerionisation_h("h + e -> h+ + 2e");
RegisterComponent<AmjuelHydIonisationIsotope<'d'>>
    registerionisation_d("d + e -> d+ + 2e");
RegisterComponent<AmjuelHydIonisationIsotope<'t'>>
    registerionisation_t("t + e -> t+ + 2e");
} // namespace

#endif // AMJUEL_HYD_IONISATION_H
