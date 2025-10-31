#pragma once

#ifndef IONISATION_H
#define IONISATION_H

#include "component.hxx"

#include "radiation.hxx"

class Ionisation : public Component {
public:
  Ionisation(std::string name, Options &options, Solver *);
  
private:
  UpdatedRadiatedPower atomic_rates {}; // Atomic rates (H.Willett)
  
  BoutReal Eionize;   // Energy loss per ionisation [eV]

  BoutReal Tnorm, Nnorm, FreqNorm; // Normalisations

  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<Ionisation> registersolverionisation("ionisation");
}

#endif // IONISATION_H
