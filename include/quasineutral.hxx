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
struct Quasineutral : public Component {
  /// Inputs
  /// ------
  ///
  /// name     Short name for species e.g. "e"
  /// options  Component configuration options
  ///   - <name>
  ///     - charge   Required to have a particle charge
  ///
  Quasineutral(std::string name, Options &alloptions, Solver *UNUSED(solver));

  /// 
  /// 
  void transform(Options &state) override;
private:
  std::string name; ///< Name of this species
  BoutReal charge;  ///< The charge of this species
  BoutReal AA;      ///< Atomic mass
};

namespace {
RegisterComponent<Quasineutral>
    registercomponentquasineutral("quasineutral");
}

#endif // QUASINEUTRAL
