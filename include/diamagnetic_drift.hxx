#pragma once
#ifndef DIAMAGNETIC_DRIFT_H
#define DIAMAGNETIC_DRIFT_H

#include "component.hxx"
#include "guarded_options.hxx"

#include <bout/bout_types.hxx>
#include <bout/field2d.hxx>
#include <bout/options.hxx>
#include <bout/vector2d.hxx>

#include <string>

/// Calculate diamagnetic flows

struct DiamagneticDrift : public Component {
  DiamagneticDrift(std::string name, Options& options, [[maybe_unused]] Solver* solver);

private:
  Vector2D Curlb_B;
  bool bndry_flux;      /// Allow boundary fluxes?
  bool divergence_form; ///< Use divergence form?

  bool average_core; ///< Average around core boundary?

  /// For every species, if it has:
  ///  - temperature
  ///  - charge
  ///
  /// Modifies:
  ///  - density_source
  ///  - energy_source
  ///  - momentum_source
  void transform_impl(GuardedOptions& state) override;

  /// Smooth source around core boundary
  /// Modifies the input field
  void coreAverage(Field3D& f);
  Field2D cell_volume;
  BoutReal core_ring_volume{0.0};
};

namespace {
RegisterComponent<DiamagneticDrift> registercomponentdiamagnetic("diamagnetic_drift");
}

#endif // DIAMAGNETIC_DRIFT_H
