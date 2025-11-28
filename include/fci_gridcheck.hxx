#pragma once
#ifndef FCI_GRIDCHECK_H
#define FCI_GRIDCHECK_H

#include "component.hxx"

/// 2D closure, modelling currents through a sheath
///
/// This should only be used where one grid cell is used in y (ny=1).
/// For domains with multiple Y points, use sheath_boundary
struct FCIGridcheck : public Component {
  /// Inputs
  ///  - units
  ///    - meters    Length normalisation
  ///  - <name>
  ///    - connection_length    Parallel connection length in meters
  ///
  FCIGridcheck(std::string name, Options &options, Solver *);

  /// Inputs
  /// - fields
  ///   - phi      Electrostatic potential
  ///
  /// Optional inputs
  /// - species
  ///   - density
  ///   - pressure
  ///
  /// Modifies
  /// - species
  ///   - e
  ///     - density_source   (If density present)
  ///   - density_source and energy_source (If sinks=true)
  /// - fields
  ///   - DivJdia     Divergence of current
  ///
  void transform(Options &state) override;

  void outputVars(Options &state) override;
  
private:

};

namespace {
RegisterComponent<FCIGridcheck>
    registercomponentsheathclosure("fci_gridcheck");
}


#endif // FCI_GRIDCHECK_H
