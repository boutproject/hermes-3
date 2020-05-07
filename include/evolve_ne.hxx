#pragma once
#ifndef EVOLVE_NE_H
#define EVOLVE_NE_H

#include "component.hxx"

/// Evolve electron density in time
/// 
struct EvolveNe : public Component {
  EvolveNe(std::string name, Options &options, Solver *solver);

  /// This sets in the state
  /// - species
  ///   - e
  ///     - density
  void transform(Options &state) override;

  /// Calculate ddt(Ne).
  ///
  /// Requires state components
  /// - species
  ///   - e
  ///     - density
  ///
  /// Optional components
  /// - species
  ///   - e
  ///     - velocity        If included, requires sound_speed or temperature
  ///     - density_source
  /// - fields
  ///   - phi 
  void finally(const Options &state) override;
private:
  Field3D Ne;  // Electron density (normalised, evolving)

  bool ne_bndry_flux;
  bool poloidal_flows;
  
  BoutReal anomalous_D; // Anomalous diffusion (axisymmetric)
  bool low_n_diffuse;   // Parallel diffusion at low density
  bool low_n_diffuse_perp;  // Perpendicular diffusion at low density
  bool hyper_z;  // Hyper-diffusion in Z
};

namespace {
RegisterComponent<EvolveNe> registercomponentevolvene("evolve_ne");
}


#endif // EVOLVE_NE_H
