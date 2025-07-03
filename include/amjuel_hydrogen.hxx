#pragma once
#ifndef AMJUEL_HYD_RECOMBINATION_H
#define AMJUEL_HYD_RECOMBINATION_H

#include <bout/constants.hxx>

#include "amjuel_reaction.hxx"

/**
 * @brief Base class for ionisation and recombination reaction Components for Hydrogen
 * Isotopes.
 *
 * @tparam Isotope char representing the isotope type; 'h', 'd' and 't'. Allows all three
 * species to be treated with the same code.
 */
template <char Isotope>
struct AmjuelHydIsotopeReaction : public AmjuelReaction {
  AmjuelHydIsotopeReaction(std::string name, std::string reaction_type,
                           std::string from_species, std::string to_species,
                           Options& alloptions)
      : AmjuelReaction(name, std::string("hyd_") + reaction_type, from_species,
                       to_species, alloptions),
        reaction_type(reaction_type) {}

  void outputVars(Options& state) override {
    AUTO_TRACE();
    // Normalisations
    BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation
    auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
    auto Cs0 = get<BoutReal>(state["Cs0"]);

    if (this->diagnose) {
      // Save particle, momentum and energy channels
      set_with_attrs(
          state[std::string("S") + this->reaction_suffix], this->S,
          {{"time_dimension", "t"},
           {"units", "m^-3 s^-1"},
           {"conversion", Nnorm * Omega_ci},
           {"standard_name", "particle source"},
           {"long_name", std::string("Particle source due to ") + this->reaction_type
                             + " of " + this->from_species + " to " + this->to_species},
           {"source", this->amjuel_src}});

      set_with_attrs(
          state[std::string("F") + this->reaction_suffix], this->F,
          {{"time_dimension", "t"},
           {"units", "kg m^-2 s^-2"},
           {"conversion", SI::Mp * Nnorm * Cs0 * Omega_ci},
           {"standard_name", "momentum transfer"},
           {"long_name", (std::string("Momentum transfer due to ") + this->reaction_type
                          + " of " + this->from_species + " to " + this->to_species)},
           {"source", this->amjuel_src}});

      set_with_attrs(
          state[std::string("E") + this->reaction_suffix], this->E,
          {{"time_dimension", "t"},
           {"units", "W / m^3"},
           {"conversion", Pnorm * Omega_ci},
           {"standard_name", "energy transfer"},
           {"long_name", (std::string("Energy transfer due to ") + this->reaction_type
                          + " of " + this->from_species + " to " + this->to_species)},
           {"source", this->amjuel_src}});

      set_with_attrs(
          state[std::string("R") + this->radiation_suffix], this->R,
          {{"time_dimension", "t"},
           {"units", "W / m^3"},
           {"conversion", Pnorm * Omega_ci},
           {"standard_name", "radiation loss"},
           {"long_name", (std::string("Radiation loss due to ") + this->reaction_type
                          + " of " + this->from_species + " to " + this->to_species)},
           {"source", this->amjuel_src}});
    }
  }

protected:
  // Strings used in outputVars to avoid duplicating code for ionisation/recombination
  const std::string reaction_type;
  std::string reaction_suffix;
  std::string radiation_suffix;
};

/**
 * @brief Component for Hydrogen recombination based on Amjuel rates. Includes both
 * radiative and 3-body recombination.
 *
 * @tparam Isotope char representing the isotope type; 'h', 'd' and 't'.
 */
template <char Isotope>
struct AmjuelHydRecombinationIsotope : public AmjuelHydIsotopeReaction<Isotope> {
  AmjuelHydRecombinationIsotope(std::string name, Options& alloptions, Solver*)
      : AmjuelHydIsotopeReaction<Isotope>(name, "recombination", {Isotope, '+'},
                                          {Isotope}, alloptions) {

    this->reaction_suffix = std::string({Isotope}) + "+_rec";
    this->radiation_suffix = std::string({Isotope}) + "+_rec";

    this->rate_multiplier = alloptions[{Isotope}]["K_rec_multiplier"]
                                .doc("Scale the recombination rate by this factor")
                                .withDefault<BoutReal>(1.0);

    this->radiation_multiplier =
        alloptions[{Isotope}]["R_rec_multiplier"]
            .doc("Scale the recombination radiation (incl. 3 body) rate by this factor")
            .withDefault<BoutReal>(1.0);
  }

  void set_diagnostic_fields(Field3D& reaction_rate, Field3D& momentum_exchange,
                             Field3D& energy_exchange,
                             Field3D& energy_loss) override final {
    this->S = -reaction_rate;
    this->F = -momentum_exchange;
    this->E = -energy_exchange;
    this->R = -energy_loss;
  }
};

/**
 * @brief Component for Hydrogen ionisation based on Amjuel rates.
 *
 * @tparam Isotope char representing the isotope type; 'h', 'd' and 't'.
 */
template <char Isotope>
struct AmjuelHydIonisationIsotope : public AmjuelHydIsotopeReaction<Isotope> {
  AmjuelHydIonisationIsotope(std::string name, Options& alloptions, Solver*)
      : AmjuelHydIsotopeReaction<Isotope>(name, "ionisation", {Isotope}, {Isotope, '+'},
                                          alloptions) {

    this->reaction_suffix = std::string({Isotope}) + "+_iz";
    this->radiation_suffix = std::string({Isotope}) + "+_ex";

    this->rate_multiplier = alloptions[{Isotope}]["K_iz_multiplier"]
                                .doc("Scale the ionisation rate by this factor")
                                .withDefault<BoutReal>(1.0);

    this->radiation_multiplier = alloptions[{Isotope}]["R_ex_multiplier"]
                                     .doc("Scale the ionisation excitation/de-excitation "
                                          "radiation rate by this factor")
                                     .withDefault<BoutReal>(1.0);
  }

  void set_diagnostic_fields(Field3D& reaction_rate, Field3D& momentum_exchange,
                             Field3D& energy_exchange,
                             Field3D& energy_loss) override final {
    this->S = reaction_rate;
    this->F = momentum_exchange;
    this->E = energy_exchange;
    this->R = -energy_loss;
  }
};

namespace {
/// Register three ionisation components, one for each hydrogen isotope
/// so no isotope dependence included.
RegisterComponent<AmjuelHydIonisationIsotope<'h'>>
    registerionisation_h("h + e -> h+ + 2e");
RegisterComponent<AmjuelHydIonisationIsotope<'d'>>
    registerionisation_d("d + e -> d+ + 2e");
RegisterComponent<AmjuelHydIonisationIsotope<'t'>>
    registerionisation_t("t + e -> t+ + 2e");

/// Register three recombination components, one for each hydrogen isotope
/// so no isotope dependence included.
RegisterComponent<AmjuelHydRecombinationIsotope<'h'>>
    register_recombination_h("h+ + e -> h");
RegisterComponent<AmjuelHydRecombinationIsotope<'d'>>
    register_recombination_d("d+ + e -> d");
RegisterComponent<AmjuelHydRecombinationIsotope<'t'>>
    register_recombination_t("t+ + e -> t");

} // namespace

#endif // AMJUEL_HYD_RECOMBINATION_H
