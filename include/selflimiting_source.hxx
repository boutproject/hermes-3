#pragma once
#ifndef SELFLIMITING_SOURCE_H
#define SELFLIMITING_SOURCE_H

#include "component.hxx"

/// Add selflimiting source for density and pressure evolution
///
struct SelflimitingSource : public NamedComponent<SelflimitingSource> {
  /// # Inputs
  /// See documentation
  /// https://hermes3.readthedocs.io/en/latest/feedback_control.html

  SelflimitingSource(std::string name, Options& alloptions, Solver*);

  void outputVars(Options& state) override;

  static constexpr auto type = "selflimiting_source";

private:
  std::string name; ///< Species name

  bool diagnose;        ///< Outputting diagnostics?
  Field3D densityform;  ///< Tunable form of the density source
  Field3D pressureform; ///< Tunable form of the pressure source

  BoutReal density_setpoint;     ///< Density above the source will be dampened
  BoutReal temperature_setpoint; ///< Temperature above the source will be dampened
  BoutReal density_overshoot;    ///< Density interval in which the source is dampened
  BoutReal
      temperature_overshoot; ///< Temperature interval in which the source is dampened

  Field3D density_source;
  Field3D energy_source;

  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<SelflimitingSource> registercomponentselflimitingsource;
}

#endif // SELFLIMITING_SOURCE_H
