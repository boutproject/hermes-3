#pragma once
#ifndef ANOMALOUS_DIFFUSION_H
#define ANOMALOUS_DIFFUSION_H

#include "component.hxx"

/// Add anomalous diffusion of density, momentum and energy
///
struct AnomalousDiffusion : public Component {
  AnomalousDiffusion(std::string name, Options &alloptions, Solver *);

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

private:
  std::string name; ///< Species name

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
