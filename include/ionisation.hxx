#pragma once

#ifndef IONISATION_H
#define IONISATION_H

#include "component.hxx"

#include "radiation.hxx"

class Ionisation : public Component {
public:
  Ionisation(std::string name, Options &options, Solver *);
  void transform(Options &state) override;
  
private:
  UpdatedRadiatedPower atomic_rates {}; // Atomic rates (H.Willett)
  
  BoutReal Eionize;   // Energy loss per ionisation [eV]

  BoutReal Tnorm, Nnorm, FreqNorm; // Normalisations
};

namespace {
RegisterComponent<Ionisation> registersolverionisation("ionisation");
}

#endif // IONISATION_H
