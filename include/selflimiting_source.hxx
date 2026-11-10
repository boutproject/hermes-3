#pragma once
#ifndef SELFLIMITING_SOURCE_H
#define SELFLIMITING_SOURCE_H

#include "component.hxx"

/// Add selflimiting source for density and pressure evolution
///
struct SelflimitingSource : public NamedComponent<SelflimitingSource> {
  /// # Inputs
  ///
  /// - <name>
  ///   - anomalous_D    This overrides D_<name> mesh input
  ///   - anomalous_chi  This overrides chi_<name>
  ///   - anomalous_nu   Overrides nu_<name>
  ///   - anomalous_sheath_flux  Allow anomalous flux into sheath?
  //                             Default false.
  SelflimitingSource(std::string name, Options& alloptions, Solver*);

  void outputVars(Options& state) override;

  static constexpr auto type = "selflimiting_source";

private:
  std::string name; ///< Species name

  bool diagnose;                     ///< Outputting diagnostics?
  Field3D selflimiting_densityform;  ///< Tunable form of the density source
  Field3D selflimiting_pressureform; ///< Tunable form of the pressure source

  BoutReal selflimiting_density;     ///< Density above the source will be dampened
  BoutReal selflimiting_temperature; ///< Temperature above the source will be dampened
  BoutReal
      selflimiting_density_overshoot; ///< Density interval in which the source is dampened
  BoutReal
      selflimiting_temperature_overshoot; ///< Temperature interval in which the source is dampened

  Field3D density_source;
  Field3D energy_source;

  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<SelflimitingSource> registercomponentselflimitingsource;
}

#endif // SELFLIMITING_SOURCE_H
