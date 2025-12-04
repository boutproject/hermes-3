#pragma once

#ifndef IONISATION_H
#define IONISATION_H

#include "component.hxx"

class Ionisation : public Component {
public:
  Ionisation(std::string name, Options &options, Solver *);
  void transform(Options &state) override;
  
private:
  BoutReal Eionize;   // Energy loss per ionisation [eV]

  BoutReal Tnorm, Nnorm, FreqNorm; // Normalisations
};

namespace {
RegisterComponent<Ionisation> registersolverionisation("ionisation");
}

#endif // IONISATION_H
