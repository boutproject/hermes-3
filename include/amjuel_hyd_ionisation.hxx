#pragma once
#ifndef AMJUEL_HYD_IONISATION_H
#define AMJUEL_HYD_IONISATION_H

#include <bout/constants.hxx>

#include "amjuel_reaction.hxx"

/// Hydrogen ionisation, Amjuel rates
struct AmjuelHydIonisation : public AmjuelReaction {
  AmjuelHydIonisation(std::string name, Options& alloptions, Solver* solver)
      : AmjuelReaction(name, alloptions, solver) {}

  void calculate_rates(Options& electron, Options& atom, Options& ion,
                       Field3D& reaction_rate, Field3D& momentum_exchange,
                       Field3D& energy_exchange, Field3D& energy_loss, BoutReal& multiplier);
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

    multiplier = alloptions[name]["scale_ionisation"]
                           .doc("Scale the ionisation rate by this factor")
                           .withDefault<BoutReal>(1.0);
  }

  void transform(Options& state) override {
    Options& electron = state["species"]["e"];
    Options& atom = state["species"][{Isotope}];     // e.g. "h"
    Options& ion = state["species"][{Isotope, '+'}]; // e.g. "h+"
    Field3D reaction_rate, momentum_exchange, energy_exchange, energy_loss;

    calculate_rates(electron, atom, ion, reaction_rate, momentum_exchange,
                    energy_exchange, energy_loss, multiplier);

    if (diagnose) {
      S = reaction_rate;
      F = momentum_exchange;
      E = energy_exchange;
      R = -energy_loss;
    }
  }

  void outputVars(Options& state) override {
    AUTO_TRACE();
    // Normalisations
    auto Nnorm = get<BoutReal>(state["Nnorm"]);
    auto Tnorm = get<BoutReal>(state["Tnorm"]);
    BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation
    auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
    auto Cs0 = get<BoutReal>(state["Cs0"]);

    if (diagnose) {
      // Save particle, momentum and energy channels

      std::string atom{Isotope};
      std::string ion{Isotope, '+'};

      set_with_attrs(state[{'S', Isotope, '+', '_', 'i', 'z'}], S,
                     {{"time_dimension", "t"},
                      {"units", "m^-3 s^-1"},
                      {"conversion", Nnorm * Omega_ci},
                      {"standard_name", "particle source"},
                      {"long_name", std::string("Ionisation of ") + atom + " to " + ion},
                      {"source", "amjuel_hyd_ionisation"}});

      set_with_attrs(
          state[{'F', Isotope, '+', '_', 'i', 'z'}], F,
          {{"time_dimension", "t"},
           {"units", "kg m^-2 s^-2"},
           {"conversion", SI::Mp * Nnorm * Cs0 * Omega_ci},
           {"standard_name", "momentum transfer"},
           {"long_name", (std::string("Momentum transfer due to ionisation of ") + atom
                          + " to " + ion)},
           {"source", "amjuel_hyd_ionisation"}});

      set_with_attrs(state[{'E', Isotope, '+', '_', 'i', 'z'}], E,
                     {{"time_dimension", "t"},
                      {"units", "W / m^3"},
                      {"conversion", Pnorm * Omega_ci},
                      {"standard_name", "energy transfer"},
                      {"long_name", (std::string("Energy transfer due to ionisation of ")
                                     + atom + " to " + ion)},
                      {"source", "amjuel_hyd_ionisation"}});

      set_with_attrs(state[{'R', Isotope, '+', '_', 'e', 'x'}], R,
                     {{"time_dimension", "t"},
                      {"units", "W / m^3"},
                      {"conversion", Pnorm * Omega_ci},
                      {"standard_name", "radiation loss"},
                      {"long_name", (std::string("Radiation loss due to ionisation of ")
                                     + atom + " to " + ion)},
                      {"source", "amjuel_hyd_ionisation"}});
    }
  }

private:
  bool diagnose; ///< Outputting diagnostics?
  BoutReal multiplier; ///< Scaling factor on reaction rate
  Field3D S;     ///< Particle exchange
  Field3D F;     ///< Momentum exchange
  Field3D E;     ///< Energy exchange
  Field3D R;     ///< Radiation loss
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
