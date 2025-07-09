#pragma once
#ifndef CENTRIFUGAL_DRIFT_H
#define CENTRIFUGAL_DRIFT_H

#include "component.hxx"

/// Calculate centrifugal drifts and currents
struct CentrifugalDrift : public Component {
  CentrifugalDrift(std::string name, Options& alloptions, Solver* UNUSED(solver));

  /// Inputs:
  /// - species
  ///   - ...  All species with both charge and mass
  ///     - AA
  ///     - charge
  ///     - density
  ///     - pressure
  ///     - momentum
  ///
  /// Modifies:
  /// - species
  ///  - ...
  ///    - density_source
  ///    - energy_source
  ///    - momentum_source
  ///
  /// - fields
  ///   - DivJextra
  void transform(Options& state) override;

  void outputVars(Options& state) override;

private:
  Vector2D acceleration;
  Vector2D axb_B; // (acceleration x b) / B

  bool bndry_flux; ///< Allow fluxes through boundary?
};

namespace {
RegisterComponent<CentrifugalDrift> registercomponentcentrifugal("centrifugal_drift");
}

#endif // CENTRIFUGAL_DRIFT_H
