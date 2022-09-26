#pragma once
#ifndef AMJUEL_HELIUM_H
#define AMJUEL_HELIUM_H

#include "amjuel_reaction.hxx"


/// e + he -> he+ + 2e
/// Amjuel reaction 2.3.9a, page 161
/// Not resolving metastables, only transporting ground state
struct AmjuelHeIonisation01 : public AmjuelReaction {
  AmjuelHeIonisation01(std::string name, Options& alloptions, Solver* solver)
      : AmjuelReaction(name, alloptions, solver) {}

  void calculate_rates(Options& state, Field3D &reaction_rate);

  void transform(Options& state) override{
    Field3D reaction_rate;
    calculate_rates(state, reaction_rate);
  };
};


/// e + he+ -> he
/// Amjuel reaction 2.3.13a
/// Not resolving metastables. Includes radiative + threebody + dielectronic.
/// Fujimoto Formulation II
struct AmjuelHeRecombination10 : public AmjuelReaction {
  AmjuelHeRecombination10(std::string name, Options& alloptions, Solver* solver)
      : AmjuelReaction(name, alloptions, solver) {}

  void calculate_rates(Options& state, Field3D &reaction_rate);

  void transform(Options& state) override{
    Field3D reaction_rate;
    calculate_rates(state, reaction_rate);
  };
};


/// e + he+ -> he+2 + 2e
/// Amjuel reaction 2.2C, page 189
// NOTE: Currently missing  energy loss / radiation data
// struct AmjuelHeIonisation12 : public AmjuelReaction {
//   AmjuelHeIonisation12(std::string name, Options& alloptions, Solver* solver)
//       : AmjuelReaction(name, alloptions, solver) {}

//   void transform(Options& state) override;
// };


namespace {
RegisterComponent<AmjuelHeIonisation01> register_ionisation_01("he + e -> he+ + 2e");
RegisterComponent<AmjuelHeRecombination10> register_recombination_10("he+ + e -> he");
//-RegisterComponent<AmjuelHeIonisation12> register_ionisation_12("e + he+ -> he+2 + 2e");
} // namespace

#endif // AMJUEL_HELIUM_H
