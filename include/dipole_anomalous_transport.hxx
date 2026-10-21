#pragma once
#ifndef DIPOLE_ANOMALOUS_DIFFUSION_H
#define DIPOLE_ANOMALOUS_DIFFUSION_H

#include "component.hxx"

/// Anomalous and quasilinear (flute interchange) transport in a dipole,
/// following OSTR-Mod.1 section 4.4
///
/// - dipole_anomalous_D: particle diffusion of N, carrying momentum and
///   energy as in anomalous_diffusion
/// - dipole_anomalous_chi: thermal diffusion of T, as in anomalous_diffusion
/// - dipole_quasilinear_D: quasilinear particle flux F_n (eq. 52) and
///   entropy flux F_s (eq. 53), giving the density, energy (eq. 57) and
///   momentum sources
///
struct DipoleAnomalousDiffusion : public NamedComponent<DipoleAnomalousDiffusion> {
  /// # Inputs
  ///
  /// - <name>
  ///   - dipole_anomalous_D        Anomalous particle diffusion
  ///   - dipole_anomalous_chi      Anomalous thermal diffusion
  ///   - dipole_anomalous_model    Scaling of dipole_anomalous_D/chi. Default 1:
  ///                               D = value/(R^2 B^2), value in [Wb^2 rad^-2 s^-1]
  ///   - dipole_quasilinear_D      Quasilinear diffusion flux function Dpsi
  ///   - dipole_quasilinear_model  Scaling of dipole_quasilinear_D. Default 1:
  ///                               D = Dpsi/(R^2 B^2), Dpsi in [Wb^2 rad^-2 s^-1]
  ///   - dipole_gamma          Adiabatic index. Default 5/3
  ///   - dipole_upwind         Use upwind operators? Default true
  DipoleAnomalousDiffusion(std::string name, Options &alloptions, Solver *);

  static constexpr auto type = "dipole_anomalous_diffusion";

  void outputVars(Options &state) override;

private:
  std::string name; ///< Species name

  bool diagnose; ///< Outputting diagnostics?
  Field2D dipole_anomalous_D; ///< Anomalous density diffusion coefficient
  Field2D dipole_anomalous_chi; ///< Anomalous thermal diffusion coefficient
  Field2D dipole_quasilinear_D; ///< Quasilinear particle diffusion coefficient
  Field2D B2D;
  bool dipole_inner_noflux; ///< No quasilinear flux through the inner boundary?

  BoutReal dipole_gamma; ///< Adiabatic index
  bool dipole_upwind;    ///< Use upwind gradients for dipole fluxes?
  Field3D flow_xlow, flow_ylow; ///< Quasilinear particle flows through cell faces
  Field3D F_n; ///< Quasilinear density source, eq. (52)
  Field3D F_s; ///< Quasilinear entropy density source, eq. (53)

  /// Inputs
  /// - species
  ///   - <name>
  ///     - AA
  ///     - density
  ///     - temperature  (optional)
  ///     - pressure     (optional)
  ///     - velocity     (optional)
  ///
  /// Sets in the state
  ///
  /// - species
  ///   - <name>
  ///     - density_source, particle_flow_xlow/ylow
  ///     - momentum_source, momentum_flow_xlow/ylow
  ///     - energy_source, energy_flow_xlow/ylow
  ///
  void transform_impl(GuardedOptions& state) override;
};
namespace {
RegisterComponent<DipoleAnomalousDiffusion> registercomponentdipoleanomalousdiffusion;
}

#endif // DIPOLE_ANOMALOUS_DIFFUSION_H
