#pragma once
#ifndef AMJUEL_HYDROGEN_H
#define AMJUEL_HYDROGEN_H

#include <bout/constants.hxx>

#include "amjuel_reaction.hxx"

static std::map<std::string, std::string> long_reaction_types_map = {
    {"cx", "charge exchange"}, {"iz", "ionisation"}, {"rec", "recombination"}};

/**
 * @brief Base class for ionisation and recombination reaction Components for Hydrogen
 * Isotopes.
 *
 * @tparam Isotope char representing the isotope type; 'h', 'd' and 't'. Allows all three
 * species to be treated with the same code.
 */
template <char Isotope>
struct AmjuelHydIsotopeReaction : public AmjuelReaction {
  AmjuelHydIsotopeReaction(std::string name, std::string short_reaction_type,
                           std::string amjuel_label, const std::string & reaction_str, std::string from_species,
                           std::string to_species, Options& alloptions)
    : AmjuelReaction(name, short_reaction_type, amjuel_label, reaction_str, from_species, to_species,
                       alloptions) {
    if (this->diagnose) {
      // Set up diagnostics

      // N.B. The first three diagnostics are named using the ion for both izn and rec,
      // but the sign is changed accordingly
      std::string ion_sp = fmt::format("{:s}+", std::string({Isotope}));
      std::string default_diag_suffix =
          fmt::format("{:s}_{:s}", ion_sp, short_reaction_type);
      std::string rad_diag_suffix = default_diag_suffix;
      DiagnosticTransformerType default_transformer = negate;
      DiagnosticTransformerType rad_transformer = identity;

      // For izn, tweak R diagnostic name and reverse the signs of the S, F, E diagnostics
      if (short_reaction_type == "iz") {
        rad_diag_suffix = fmt::format("{:s}_ex", ion_sp);
        default_transformer = identity;
      }

      auto heavy_products =
          this->parser->get_species(species_filter::heavy, species_filter::products);
      ASSERT1(heavy_products.size() == 1);
      std::string heavy_product = heavy_products[0];
      std::string long_reaction_type = long_reaction_types_map.at(short_reaction_type);
      add_diagnostic(
          heavy_product, fmt::format("S{:s}", default_diag_suffix),
          fmt::format("Particle source due to {:s} of {:s} to {:s}", long_reaction_type,
                      this->from_species, this->to_species),
          ReactionDiagnosticType::density_src, this->amjuel_src, default_transformer);

      add_diagnostic(
          heavy_product, fmt::format("F{:s}", default_diag_suffix),
          fmt::format("Momentum transfer due to {:s} of {:s} to {:s}", long_reaction_type,
                      this->from_species, this->to_species),
          ReactionDiagnosticType::momentum_src, this->amjuel_src, default_transformer);

      add_diagnostic(
          heavy_product, fmt::format("E{:s}", default_diag_suffix),
          fmt::format("Energy transfer due to {:s} of {:s} to {:s}", long_reaction_type,
                      this->from_species, this->to_species),
          ReactionDiagnosticType::energy_src, this->amjuel_src, default_transformer);

      add_diagnostic(
          "e", fmt::format("R{:s}", rad_diag_suffix),
          fmt::format("Radiation loss due to {:s} of {:s} to {:s}", long_reaction_type,
                      this->from_species, this->to_species),
          ReactionDiagnosticType::energy_loss, this->amjuel_src, rad_transformer);
    }
  }
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
    : AmjuelHydIsotopeReaction<Isotope>(name, "rec", "H.x_2.1.8", fmt::format("{}+ + e -> {}", Isotope, Isotope), {Isotope, '+'},
                                          {Isotope}, alloptions) {

    this->rate_multiplier = alloptions[{Isotope}]["K_rec_multiplier"]
                                .doc("Scale the recombination rate by this factor")
                                .withDefault<BoutReal>(1.0);

    this->radiation_multiplier =
        alloptions[{Isotope}]["R_rec_multiplier"]
            .doc("Scale the recombination radiation (incl. 3 body) rate by this factor")
            .withDefault<BoutReal>(1.0);
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
      : AmjuelHydIsotopeReaction<Isotope>(name, "iz", "H.x_2.1.5", fmt::format("{} + e -> {}+ + 2e", Isotope, Isotope), {Isotope},
                                          {Isotope, '+'}, alloptions) {
    this->rate_multiplier = alloptions[{Isotope}]["K_iz_multiplier"]
                                .doc("Scale the ionisation rate by this factor")
                                .withDefault<BoutReal>(1.0);

    this->radiation_multiplier = alloptions[{Isotope}]["R_ex_multiplier"]
                                     .doc("Scale the ionisation excitation/de-excitation "
                                          "radiation rate by this factor")
                                     .withDefault<BoutReal>(1.0);
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
