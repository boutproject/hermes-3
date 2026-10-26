#pragma once
#ifndef DIAMAGNETIC_DRIFT_H
#define DIAMAGNETIC_DRIFT_H

#include <bout/vectormetric.hxx>

#include "component.hxx"

/// Calculate diamagnetic flows

struct DiamagneticDrift : public NamedComponent<DiamagneticDrift> {
  DiamagneticDrift(std::string name, Options& options, Solver* UNUSED(solver));

  static constexpr auto type = "diamagnetic_drift";

private:
  VectorMetric Curlb_B;
  bool bndry_flux;
  Coordinates::FieldMetric diamag_form;

  /// For every species, if it has:
  ///  - temperature
  ///  - charge
  ///
  /// Modifies:
  ///  - density_source
  ///  - energy_source
  ///  - momentum_source
  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<DiamagneticDrift> registercomponentdiamagnetic;
}

#endif // DIAMAGNETIC_DRIFT_H
