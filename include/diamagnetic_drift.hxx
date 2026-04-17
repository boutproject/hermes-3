#pragma once
#ifndef DIAMAGNETIC_DRIFT_H
#define DIAMAGNETIC_DRIFT_H

#include "component.hxx"

/// Calculate diamagnetic flows

struct DiamagneticDrift : public Component {
  DiamagneticDrift(std::string name, Options &options, Solver *UNUSED(solver));

  /// For every species, if it has:
  ///  - temperature
  ///  - charge
  ///
  /// Modifies:
  ///  - density_source
  ///  - energy_source
  ///  - momentum_source
  void transform(Options &state) override;

private:
  Vector2D Curlb_B;
  bool bndry_flux;      /// Allow boundary fluxes?
  bool divergence_form; ///< Use divergence form?

  bool average_core; ///< Average around core boundary?

  /// Smooth source around core boundary
  /// Modifies the input field
  void coreAverage(Field3D& f);
  Field2D cell_volume;
  BoutReal core_ring_volume;  

};

namespace {
RegisterComponent<DiamagneticDrift> registercomponentdiamagnetic("diamagnetic_drift");
}

#endif // DIAMAGNETIC_DRIFT_H


