#pragma once
#ifndef BRAGINSKII_ELECTRON_VISCOSITY_H
#define BRAGINSKII_ELECTRON_VISCOSITY_H

#include <string>

#include <bout/bout_types.hxx>
#include <bout/field3d.hxx>
#include <bout/options.hxx>

#include "component.hxx"

/// Electron viscosity
///
/// Adds Braginskii parallel electron viscosity, with SOLPS-style
/// viscosity flux limiter
///
/// Needs to be calculated after collisions, because collision
/// frequency is used to calculate parallel viscosity
///
/// References
///  - https://farside.ph.utexas.edu/teaching/plasma/lectures1/node56.html
///
struct BraginskiiElectronViscosity : public NamedComponent<BraginskiiElectronViscosity> {
  /// Inputs
  /// - <name>
  ///   - diagnose: bool, default false
  ///     Output diagnostic SNVe_viscosity?
  ///   - eta_limit_alpha: float, default -1.0
  ///     Flux limiter coefficient. < 0 means no limiter
  ///   - eta_free_streaming: bool, default false
  ///     Reduce eta by 1 / (1 + lambda_e / (c L_par)), lambda_e = v_te tau_e,
  ///     1/L_par = 1/L_line (+ |Grad_par Te|/Te if eta_free_streaming_gradient)
  ///   - eta_free_streaming_gradient: bool, default false
  ///   - eta_free_streaming_c: float, default 1
  ///   - eta_free_streaming_length: float [m], L_line. Default: field line
  ///     length from the grid (needs a grid without branch cuts)
  BraginskiiElectronViscosity(const std::string& name, Options& alloptions, Solver*);

  void outputVars(Options& state) override;

  static constexpr auto type = "braginskii_electron_viscosity";

private:
  BoutReal eta_limit_alpha; ///< Flux limit coefficient
  bool eta_free_streaming;  ///< Reduce eta when lambda_e is not << L_par?
  BoutReal eta_free_streaming_c; ///< Coefficient c in 1 / (1 + lambda_e / (c L_par))
  bool eta_free_streaming_gradient; ///< Include |Grad_par Te|/Te in 1/L_par?
  Field2D L_line;           ///< Field line length
  bool diagnose;            ///< Output viscosity diagnostic?
  Field3D viscosity;        ///< The viscosity momentum source

  bout::FieldMetricParallel sqrtB; ///< Magnetic field variables

  /// Inputs
  /// - species
  ///   - e
  ///     - pressure  (skips if not present)
  ///     - velocity  (skips if not present)
  ///     - collision_frequency
  ///
  /// Sets in the state
  /// - species
  ///   - e
  ///     - momentum_source
  ///
  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<BraginskiiElectronViscosity>
    registercomponentbraginskiielectronviscosity;
}

#endif
