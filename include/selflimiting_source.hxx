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

  // Selflimiting source function
  // This function uses logistic growth to dampen this over the interval overshoot from 1 to 0
  // k = 5.0/overshoot is chosen so the function goes from 1 to ~0.1 in the inverval of the overshoot
  // This function just returns the correction, which has to be multiplied with the actual source
  Field3D selflimiting_sourceterm(const Field3D& thisfield, const BoutReal maximum,
                                  const BoutReal overshoot) {
    BoutReal k = 5.0 / overshoot;
    BoutReal x0 = maximum + overshoot / 2.0;
    return (1.0 - (1.0 / (1.0 + exp(-k * (thisfield - x0)))));
  }

  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<SelflimitingSource> registercomponentselflimitingsource;
}

#endif // SELFLIMITING_SOURCE_H
