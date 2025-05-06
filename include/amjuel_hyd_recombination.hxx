#pragma once
#ifndef AMJUEL_HYD_RECOMBINATION_H
#define AMJUEL_HYD_RECOMBINATION_H

#include <bout/constants.hxx>

#include "amjuel_reaction.hxx"

/// Hydrogen recombination, Amjuel rates
///
/// Includes both radiative and 3-body recombination

/// Hydrogen recombination
/// Templated on a char to allow 'h', 'd' and 't' species to be treated with the same code
template <char Isotope>
struct AmjuelHydRecombinationIsotope : public AmjuelReaction {
  AmjuelHydRecombinationIsotope(std::string name, Options& alloptions, Solver* solver)
      : AmjuelReaction(name, "hyd_recombination", {Isotope, '+'}, {Isotope}, alloptions,
                       solver) {

    rate_multiplier = alloptions[{Isotope}]["K_rec_multiplier"]
                          .doc("Scale the recombination rate by this factor")
                          .withDefault<BoutReal>(1.0);

    radiation_multiplier =
        alloptions[{Isotope}]["R_rec_multiplier"]
            .doc("Scale the recombination radiation (incl. 3 body) rate by this factor")
            .withDefault<BoutReal>(1.0);
  }

  void set_diagnostic_fields(Field3D& reaction_rate, Field3D& momentum_exchange,
                             Field3D& energy_exchange, Field3D& energy_loss) override final {
    S = -reaction_rate;
    F = -momentum_exchange;
    E = -energy_exchange;
    R = -energy_loss;
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
      set_with_attrs(state[{'S', Isotope, '+', '_', 'r', 'e', 'c'}], S,
                     {{"time_dimension", "t"},
                      {"units", "m^-3 s^-1"},
                      {"conversion", Nnorm * Omega_ci},
                      {"standard_name", "particle source"},
                      {"long_name", std::string("Particle source due to recombnation of ")
                                        + from_species + " to " + to_species},
                      {"source", "amjuel_hyd_recombination"}});

      set_with_attrs(
          state[{'F', Isotope, '+', '_', 'r', 'e', 'c'}], F,
          {{"time_dimension", "t"},
           {"units", "kg m^-2 s^-2"},
           {"conversion", SI::Mp * Nnorm * Cs0 * Omega_ci},
           {"standard_name", "momentum transfer"},
           {"long_name", (std::string("Momentum transfer due to recombination of ")
                          + from_species + " to " + to_species)},
           {"source", "amjuel_hyd_recombination"}});

      set_with_attrs(
          state[{'E', Isotope, '+', '_', 'r', 'e', 'c'}], E,
          {{"time_dimension", "t"},
           {"units", "W / m^3"},
           {"conversion", Pnorm * Omega_ci},
           {"standard_name", "energy transfer"},
           {"long_name", (std::string("Energy transfer due to recombination of ")
                          + from_species + " to " + to_species)},
           {"source", "amjuel_hyd_recombination"}});

      set_with_attrs(
          state[{'R', Isotope, '+', '_', 'r', 'e', 'c'}], R,
          {{"time_dimension", "t"},
           {"units", "W / m^3"},
           {"conversion", Pnorm * Omega_ci},
           {"standard_name", "radiation loss"},
           {"long_name", (std::string("Radiation loss due to recombination of ")
                          + from_species + " to " + to_species)},
           {"source", "amjuel_hyd_recombination"}});
    }
  }
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
