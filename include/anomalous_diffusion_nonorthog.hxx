#pragma once
#ifndef ANOMALOUS_DIFFUSION_NONORTHOG_H
#define ANOMALOUS_DIFFUSION_NONORTHOG_H

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
struct AnomalousDiffusionNonorthog : public Component {
  /// # Inputs
  ///
  /// - <name>
  ///   - anomalous_D    This overrides D_<name> mesh input
  ///   - anomalous_chi  This overrides chi_<name>
  ///   - anomalous_nu   Overrides nu_<name>
  ///   - anomalous_sheath_flux  Allow anomalous flux into sheath?
  //                             Default false.
  AnomalousDiffusionNonorthog(std::string name, Options &alloptions, Solver *);

  /// Inputs
  /// - species
  ///   - <name>
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
  void transform(Options &state) override;
  void outputVars(Options &state) override;

private:
  std::string name; ///< Species name

  bool diagnose; ///< Outputting diagnostics?
  bool include_D, include_chi, include_nu; ///< Which terms should be included?
  Field3D anomalous_D; ///< Anomalous density diffusion coefficient
  Field3D anomalous_chi; ///< Anomalous thermal diffusion coefficient
  Field3D anomalous_nu; ///< Anomalous momentum diffusion coefficient

  bool anomalous_sheath_flux; ///< Allow anomalous diffusion into sheath?
};



namespace {
  RegisterComponent<AnomalousDiffusionNonorthog> registercomponentanomalousdiffusionnonorthog("anomalous_diffusion_nonorthog");
}

#endif // ANOMALOUS_DIFFUSION_NONORTHOG_H
