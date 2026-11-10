#pragma once
#ifndef ADAPTIVE_SOURCE_H
#define ADAPTIVE_SOURCE_H

#include "component.hxx"

/// Add adaptive source for density and pressure evolution
///
struct AdaptiveSource : public NamedComponent<AdaptiveSource> {
  /// # Inputs
  ///
  /// - <name>
  ///   - anomalous_D    This overrides D_<name> mesh input
  ///   - anomalous_chi  This overrides chi_<name>
  ///   - anomalous_nu   Overrides nu_<name>
  ///   - anomalous_sheath_flux  Allow anomalous flux into sheath?
  //                             Default false.
  AdaptiveSource(std::string name, Options& alloptions, Solver*);

  void outputVars(Options& state) override;

  static constexpr auto type = "adaptive_source";

private:
  std::string name; ///< Species name

  bool diagnose;                 ///< Outputting diagnostics?
  Field3D adaptive_densityform;  ///< Tunable form of the density source
  Field3D adaptive_pressureform; ///< Tunable form of the pressure source

  BoutReal adaptive_density;     ///< Density above the source will be dampened
  BoutReal adaptive_temperature; ///< Temperature above the source will be dampened
  BoutReal
      adaptive_density_overshoot; ///< Density interval in which the source is dampened
  BoutReal
      adaptive_temperature_overshoot; ///< Temperature interval in which the source is dampened

  // Adaptive source function
  // This function uses logistic growth to dampen this over the interval overshoot from 1 to 0
  // k = 5.0/overshoot is chosen so the function goes from 1 to ~0.1 in the inverval of the overshoot
  // This function just returns the correction, which has to be multiplied with the actual source
  Field3D adaptive_sourceterm(const Field3D& thisfield, const BoutReal maximum,
                              const BoutReal overshoot) {
    BoutReal k = 5.0 / overshoot;
    BoutReal x0 = maximum + overshoot / 2.0;
    return (1.0 - (1.0 / (1.0 + exp(-k * (thisfield - x0)))));
  }

  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<AdaptiveSource> registercomponentadaptivesource;
}

#endif // ADAPTIVE_SOURCE_H
