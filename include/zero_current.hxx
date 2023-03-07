#pragma once
#ifndef ZERO_CURRENT
#define ZERO_CURRENT

#include "component.hxx"

/// Set the velocity of a species so that
/// there is no net current, by summing the current from
/// other species.
///
/// This is most often used in the electron species, but
/// does not need to be.
struct ZeroCurrent : public Component {
  /// Inputs
  /// ------
  ///
  /// @param name    Short name for species e.g. "e"
  /// @param alloptions   Component configuration options
  ///   - <name>
  ///     - charge   (must not be zero)
  ZeroCurrent(std::string name, Options& alloptions, Solver*);
  
  /// Required inputs
  /// - species
  ///   - <name>
  ///     - density
  ///     - charge
  ///   - <one or more other species>
  ///     - density
  ///     - velocity
  ///     - charge
  ///
  /// Sets in the state
  /// - species
  ///   - <name>
  ///     - velocity
  /// 
  void transform(Options &state) override;

  void finally(const Options &state) override {
    // Get the velocity with boundary condition applied.
    // This is for output only
    velocity = get<Field3D>(state["species"][name]["velocity"]);
  }

  void outputVars(Options &state) override;
private:
  std::string name; ///< Name of this species
  BoutReal charge;  ///< The charge of this species

  Field3D velocity; ///< Species velocity (for writing to output)
};

namespace {
RegisterComponent<ZeroCurrent> registercomponentzerocurrent("zero_current");
}

#endif // ZERO_CURRENT
