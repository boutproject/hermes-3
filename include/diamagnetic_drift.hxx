#pragma once
#ifndef DIAMAGNETIC_DRIFT_H
#define DIAMAGNETIC_DRIFT_H

#include <bout/vectormetric.hxx>

#include "component.hxx"
#include "guarded_options.hxx"

#include <bout/bout_enum_class.hxx>
#include <bout/bout_types.hxx>
#include <bout/coordinates.hxx>
#include <bout/options.hxx>
#include <bout/vector2d.hxx>

#include <string>

BOUT_ENUM_CLASS(DiamagneticForm, gradient, divergence, bracket);
BOUT_ENUM_CLASS(DiamagneticBoundary, open, closed, closed_average);

/// Calculate diamagnetic flows
struct DiamagneticDrift : public NamedComponent<DiamagneticDrift> {
  DiamagneticDrift(std::string name, Options& options, [[maybe_unused]] Solver* solver);

  static constexpr auto type = "diamagnetic_drift";

private:
  /// Calculate the diamagnetic sink for a single evolved quantity
  /// using divergence form:
  ///
  ///   factor * Div(quantity * T / q * Curlb_B)
  Field3D calculateDivergenceForm(const Field3D& quantity, const Field3D& temperature,
                                  BoutReal charge, BoutReal factor = 1.0);

  /// Calculate the diamagnetic sink for a single evolved quantity
  /// using gradient form:
  ///
  ///   factor * Curlb_B . Grad(quantity * T / q)
  Field3D calculateGradientForm(const Field3D& quantity, const Field3D& temperature,
                                BoutReal charge, BoutReal factor = 1.0);

  /// Calculate factor * 2 * [log(B), quantity * T / q] using the Arakawa bracket.
  Field3D calculateBracketForm(const Field3D& quantity, const Field3D& temperature,
                               BoutReal charge, BoutReal factor = 1.0);

  /// Calculate the local volume-weighted sum on the inner core ring.
  BoutReal coreRingSum(const Field3D& f) const;

  /// Set the inner core-ring cells to the global volume-weighted average.
  ///
  /// Requires closed_average on the core boundary for this component.
  /// global_sum is the volume-weighted sum over the full core ring.
  /// Modifies the input field in-place.
  void setCoreAverage(Field3D& f, BoutReal global_sum) const;

  VectorMetric Curlb_B;
  DiamagneticForm form;                ///< Form of the diamagnetic sink
  DiamagneticBoundary radial_boundary; ///< Treatment of radial boundary faces
  bool average_core_active;            ///< Average on this processor's core boundary?
  Field3D bracket_factor;              ///< Metric factor for the bracket form
  Field3D logB;                        ///< Log of the normalised magnetic field

  /// For every species, if it has:
  ///  - temperature
  ///  - charge
  ///
  /// Modifies:
  ///  - density_source
  ///  - energy_source
  ///  - momentum_source
  void transform_impl(GuardedOptions& state) override;
  Coordinates::FieldMetric cell_volume;
  BoutReal core_ring_volume{0.0};
};

namespace {
RegisterComponent<DiamagneticDrift> registercomponentdiamagnetic;
}

#endif // DIAMAGNETIC_DRIFT_H
