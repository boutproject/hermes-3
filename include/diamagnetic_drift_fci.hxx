#pragma once
#ifndef DIAMAGNETIC_DRIFT_FCI_H
#define DIAMAGNETIC_DRIFT_FCI_H
#include <bout/vectormetric.hxx>
#include <bout/yboundary_regions.hxx>
#include "component.hxx"
/// Calculate diamagnetic flows
struct DiamagneticDriftFCI : public Component {
  DiamagneticDriftFCI(std::string name, Options &options, Solver *UNUSED(solver));
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
  Coordinates::FieldMetric bracket_factor;
  Field3D logB;
  YBoundary yboundary;
};





namespace {
RegisterComponent<DiamagneticDriftFCI> registercomponentdiamagnetic("diamagnetic_drift_fci");
}





#endif // DIAMAGNETIC_DRIFT_FCI_H
