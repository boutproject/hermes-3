#pragma once
#ifndef ELECTROMAGNETIC_H
#define ELECTROMAGNETIC_H

#include "component.hxx"

class Laplacian;

/// Electromagnetic potential A||
///
/// Reinterprets all species' parallel momentum as a combination
/// of a parallel flow and a magnetic contribution, i.e. canonical momentum.
///
///     m n v_{||} + Z e n A_{||}
///
/// Changes the "momentum" of each species so that after this component
/// the momentuum of each species is just
///
///     m n v_{||}
///
/// This component should be run after all species have set their
/// momentum, but before the momentum is used e.g to set boundary
/// conditions.
///
/// Calculates the electromagnetic potential A_{||} using
///
/// Laplace(Apar) - alpha_em * Apar = -Ajpar
///
/// By default outputs Apar every timestep. When `diagnose = true` in
/// also saves alpha_em and Ajpar.
///
struct Electromagnetic : public Component {
  /// Options
  /// - units
  /// - <name>
  ///   - diagnose   Saves Ajpar and alpha_em time-dependent values
  ///
  Electromagnetic(std::string name, Options &options, Solver *solver);

  /// Inputs
  /// - species
  ///   - <..>      All species with charge and parallel momentum
  ///     - charge
  ///     - momentum
  ///     - density
  ///     - AA
  ///
  /// Sets
  /// - species
  ///   - <..>      All species with charge and parallel momentum
  ///     - momentum  (modifies) to m n v||
  /// - fields
  ///   - Apar      Electromagnetic potential
  ///
  void transform(Options &state) override;
private:
  Field3D Apar; // Electromagnetic potential A_||
  Field3D Ajpar; // Total parallel current density
  Field3D alpha_em; // Coefficient
  BoutReal beta_em; // Normalisation coefficient mu_0 e T n / B^2

  std::unique_ptr<Laplacian> aparSolver; // Laplacian solver in X-Z
};

namespace {
RegisterComponent<Electromagnetic> registercomponentelectromagnetic("electromagnetic");
}

#endif // ELECTROMAGNETIC_H
