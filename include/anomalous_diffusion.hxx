#pragma once
#ifndef ANOMALOUS_DIFFUSION_H
#define ANOMALOUS_DIFFUSION_H

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
struct AnomalousDiffusion : public Component {
  /// # Inputs
  ///
  /// - <name>
  ///   - anomalous_D    This overrides D_<name> mesh input
  ///   - anomalous_chi  This overrides chi_<name>
  ///   - anomalous_nu   Overrides nu_<name>
  ///   - anomalous_sheath_flux  Allow anomalous flux into sheath?
  //                             Default false.
  AnomalousDiffusion(std::string name, Options &alloptions, Solver *) {
    diagnose = alloptions[name]["diagnose"]
                   .doc("Output additional diagnostics?")
                   .withDefault<bool>(false);
  };

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

  void outputVars(Options& state) override {
    AUTO_TRACE();
    // Normalisations
    auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
    auto rho_s0 = get<BoutReal>(state["rho_s0"]);

    if (diagnose) {
      // Save particle, momentum and energy channels

      set_with_attrs(state[{std::string("anomalous_D_") + name}], anomalous_D,
                      {{"time_dimension", "t"},
                      {"units", "m^2 s^-1"},
                      {"conversion", rho_s0 * rho_s0 * Omega_ci},
                      {"standard_name", "anomalous density diffusion"},
                      {"long_name", std::string("Anomalous density diffusion of ") + name},
                      {"source", "anomalous_diffusion"}});

    }
  }

private:
  std::string name; ///< Species name

  bool diagnose; ///< Outputting diagnostics?
  bool include_D, include_chi, include_nu; ///< Which terms should be included?
  Field2D anomalous_D; ///< Anomalous density diffusion coefficient
  Field2D anomalous_chi; ///< Anomalous thermal diffusion coefficient
  Field2D anomalous_nu; ///< Anomalous momentum diffusion coefficient

  bool anomalous_sheath_flux; ///< Allow anomalous diffusion into sheath?
};



namespace {
RegisterComponent<AnomalousDiffusion> registercomponentanomalousdiffusion("anomalous_diffusion");
}

#endif // ANOMALOUS_DIFFUSION_H
