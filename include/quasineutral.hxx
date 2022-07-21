#pragma once
#ifndef QUASINEUTRAL
#define QUASINEUTRAL

#include "component.hxx"

/// Calculate density from sum of other species densities * charge
/// to ensure that net charge = 0
///
/// This is useful in simulations where multiple species are being evolved.
/// Note that only one species' density can be calculated this way,
/// and it should be calculated last once all other densities are known.
///
/// Saves the density to the output (dump) files as N<name>
///
struct Quasineutral : public Component {
  /// Inputs
  /// ------
  ///
  /// @param name     Short name for species e.g. "e"
  /// @param alloptions  Component configuration options
  ///   - <name>
  ///     - charge   Required to have a particle charge
  ///     - AA       Atomic mass
  ///
  Quasineutral(std::string name, Options &alloptions, Solver *UNUSED(solver));

  /// 
  /// Sets in state
  /// - species
  ///   - <name>
  ///     - density
  ///     - charge
  ///     - AA
  void transform(Options &state) override;

  /// Get the final density for output
  /// including any boundary conditions applied
  void finally(const Options &state) override;
private:
  std::string name; ///< Name of this species
  BoutReal charge;  ///< The charge of this species
  BoutReal AA;      ///< Atomic mass

  Field3D density;  ///< The density (for writing to output)
};

namespace {
RegisterComponent<Quasineutral>
    registercomponentquasineutral("quasineutral");
}

#endif // QUASINEUTRAL
