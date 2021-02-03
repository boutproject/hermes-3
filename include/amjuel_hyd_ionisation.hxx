#pragma once
#ifndef AMJUEL_HYD_IONISATION_H
#define AMJUEL_HYD_IONISATION_H

#include "amjuel_reaction.hxx"

/// Hydrogen ionisation, Amjuel rates
struct AmjuelHydIonisation : public AmjuelReaction {
  AmjuelHydIonisation(std::string name, Options& alloptions, Solver* solver)
      : AmjuelReaction(name, alloptions, solver) {}

  void calculate_rates(Options& electron, Options& atom, Options& ion);
};

/// Hydrogen ionisation
/// Templated on a char to allow 'h', 'd' and 't' species to be treated with the same code
template <char Isotope>
struct AmjuelHydIonisationIsotope : public AmjuelHydIonisation {
  AmjuelHydIonisationIsotope(std::string name, Options& alloptions, Solver* solver)
      : AmjuelHydIonisation(name, alloptions, solver) {}

  void transform(Options& state) override {
    Options& electron = state["species"]["e"];
    Options& atom = state["species"][{Isotope}];     // e.g. "h"
    Options& ion = state["species"][{Isotope, '+'}]; // e.g. "h+"

    calculate_rates(electron, atom, ion);
  }
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
