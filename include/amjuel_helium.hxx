#pragma once
#ifndef AMJUEL_HELIUM_H
#define AMJUEL_HELIUM_H

#include "amjuel_reaction.hxx"

/**
 * @brief Component for Helium ionisation (e + he -> he+ + 2e) based on Amjuel
 * reaction 2.3.9a, page 161. Not resolving metastables, only transporting ground state.
 *
 */
struct AmjuelHeIonisation01 : public AmjuelReaction {
  AmjuelHeIonisation01(std::string name, Options& alloptions, Solver*)
    : AmjuelReaction(name, "iz", "H.x_2.3.9a", "he + e -> he+ + 2e", "he", "he+", alloptions) {

    rate_multiplier = alloptions[std::string("he")]["ionisation_rate_multiplier"]
                          .doc("Scale the ionisation rate by this factor")
                          .withDefault<BoutReal>(1.0);

    radiation_multiplier =
        alloptions[std::string("he")]["ionisation_radiation_rate_multiplier"]
            .doc("Scale the ionisation radiation rate by this factor")
            .withDefault<BoutReal>(1.0);
  }
};

/**
 * @brief Component for Helium recombination (e + he+ -> he) based on Amjuel
 * reaction 2.3.13a, page 181 and 295. Not resolving metastables. Includes radiative +
 * threebody + dielectronic. Fujimoto Formulation II.
 *
 */
struct AmjuelHeRecombination10 : public AmjuelReaction {
  AmjuelHeRecombination10(std::string name, Options& alloptions, Solver*)
    : AmjuelReaction(name, "rec", "H.x_2.3.13a", "he+ + e -> he", "he+", "he", alloptions) {

    rate_multiplier = alloptions[name]["recombination_rate_multiplier"]
                          .doc("Scale the recombination rate by this factor")
                          .withDefault<BoutReal>(1.0);

    radiation_multiplier =
        alloptions[name]["recombination_radiation_multiplier"]
            .doc("Scale the recombination radiation rate by this factor")
            .withDefault<BoutReal>(1.0);
  }
};

// /**
//  * @brief Component for He+ ionisation (e + he+ -> he+2 + 2e) based on Amjuel
//  * reaction 2.2C, page 189. NOTE: Currently missing  energy loss / radiation data
//  *
//  */
// struct AmjuelHeIonisation12 : public AmjuelReaction {
//   AmjuelHeIonisation12(std::string name, Options& alloptions, Solver*)
//       : AmjuelReaction(name, "iz", "H.x_2.2C", "he+", "he+2", alloptions) {}
// };

namespace {
RegisterComponent<AmjuelHeIonisation01> register_ionisation_01("he + e -> he+ + 2e");
RegisterComponent<AmjuelHeRecombination10> register_recombination_10("he+ + e -> he");
//-RegisterComponent<AmjuelHeIonisation12> register_ionisation_12("e + he+ -> he+2 + 2e");
} // namespace

#endif // AMJUEL_HELIUM_H
