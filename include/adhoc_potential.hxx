#pragma once
#ifndef ADHOC_POTENTIAL_H
#define ADHOC_POTENTIAL_H

#include "component.hxx"

/// Set temperature to a fixed value
///
struct AdhocPotential : public Component {
  AdhocPotential(std::string name, Options &options, Solver *);

  /// Inputs
  /// - species
  ///   - <name>
  ///     - density (optional)
  ///
  /// Sets in the state
  ///
  /// - species
  ///   - <name>
  ///     - temperature
  ///     - pressure (if density is set)
  ///
  void transform(Options &state) override;

  void outputVars(Options &state) override;
private:
  std::string name; // Species name
  BoutReal lambda;
  Field3D T_e; ///< The normalised temperature
  Field3D phi_adhoc;
};

namespace {
RegisterComponent<AdhocPotential> registercomponentisothermal("adhoc_potential");
}

#endif // ADHOC_POTENTIAL_H
