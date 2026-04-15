#pragma once
#ifndef DIPOLE_ANOMALOUS_DIFFUSION_H
#define DIPOLE_ANOMALOUS_DIFFUSION_H

#include "component.hxx"

/// Add anomalous diffusion of density, momentum and energy
///
/// # Mesh inputs
///
/// D_<name>, chi_<name>, nu_<name>
/// e.g `D_e`, `chi_e`, `nu_e`
///
/// in units of m^2/s
///
struct DipoleAnomalousDiffusion : public Component {
  /// # Inputs
  ///
  /// - <name>
  ///   - anomalous_D    This overrides D_<name> mesh input
  ///   - anomalous_chi  This overrides chi_<name>
  ///   - anomalous_nu   Overrides nu_<name>
  ///   - anomalous_sheath_flux  Allow anomalous flux into sheath?
  //                             Default false.
  DipoleAnomalousDiffusion(std::string name, Options &alloptions, Solver *);

  void outputVars(Options &state) override;

private:
  std::string name; ///< Species name

  bool diagnose; ///< Outputting diagnostics?
  bool include_D, include_chi, include_nu; ///< Which terms should be included?
  Field2D dipole_anomalous_D; ///< Anomalous density diffusion coefficient
  Field2D dipole_anomalous_chi; ///< Anomalous thermal diffusion coefficient
  Field2D dipole_quasilinear_chi; ///< Quasilinear thermal diffusion coefficient
  Field2D dipole_quasilinear_D; ///< Quasilinear particle diffusion coefficient
  Field2D U2D;
  Field2D transport_on;
  BoutReal dipole_gamma; // Sheath heat transmission coefficient
  BoutReal density_floor; // Minimum mass density if boussinesq=false
  bool zero_inner_gradient_U; ///< Compute local U?
  bool dipole_upwind;       ///< Use upwind gradients for anomalous fluxes?
  bool entropy_conserving;
  bool dipole_anomalous;
  Field3D flow_xlow, flow_ylow; // Flows through cell faces
  bool   dipole_div_form;
  // Field2D anomalous_nu; ///< Anomalous momentum diffusion coefficient

  // bool dipole_anomalous_sheath_flux; ///< Allow anomalous diffusion into sheath?

  /// Inputs
  /// - species
  ///   - <name>
  ///     - AA
  ///     - density
  ///     - temperature  (optional)
  ///     - velocity     (optional)
  ///
  /// Sets in the state
  ///
  /// - species
  ///   - <name>
  ///     - density_source
  ///     - momentum_source
  ///     - energy_source
  ///
  void transform_impl(GuardedOptions& state) override;
};

const void compute_U2D(Field2D& U, bool local_U);
const Field2D isnegative_grad_perp(const Field2D& P);
const Field3D isnegative_grad_perp(const Field3D& P);
namespace {
RegisterComponent<DipoleAnomalousDiffusion> registercomponentdipoleanomalousdiffusion("dipole_anomalous_diffusion");
}

#endif // DIPOLE_ANOMALOUS_DIFFUSION_H
