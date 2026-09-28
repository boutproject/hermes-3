#pragma once

#ifndef IONISATION_H
#define IONISATION_H

#include "component.hxx"

class Ionisation : public NamedComponent<Ionisation> {
public:
  Ionisation(std::string name, Options& options, Solver*);

  static constexpr auto type = "ionisation";

private:
  BoutReal Eionize; // Energy loss per ionisation [eV]

  BoutReal Tnorm, Nnorm, FreqNorm; // Normalisations

  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<Ionisation> registersolverionisation;
}

#endif // IONISATION_H
