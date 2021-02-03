#pragma once
#ifndef EVOLVE_DENSITY_H
#define EVOLVE_DENSITY_H

#include "component.hxx"

/// Evolve species density in time
/// 
struct EvolveDensity : public Component {
  EvolveDensity(std::string name, Options &options, Solver *solver);

  /// This sets in the state
  /// - species
  ///   - <name>
  ///     - AA
  ///     - charge
  ///     - density
  void transform(Options &state) override;

  /// Calculate ddt(N).
  ///
  /// Requires state components
  /// - species
  ///   - <name>
  ///     - density
  ///
  /// Optional components
  /// - species
  ///   - <name>
  ///     - velocity        If included, requires sound_speed or temperature
  ///     - density_source
  /// - fields
  ///   - phi               If included, ExB drift is calculated
  void finally(const Options &state) override;
private:
  std::string name;     ///< Short name of species e.g "e"

  BoutReal charge;      ///< Species charge e.g. electron = -1
  BoutReal AA;          ///< Atomic mass e.g. proton = 1
  
  Field3D N;            ///< Species density (normalised, evolving)

  bool bndry_flux;      // Allow flows through boundaries?
  bool poloidal_flows;  // Include ExB flow in Y direction?

  bool low_n_diffuse;   // Parallel diffusion at low density
  bool low_n_diffuse_perp;  // Perpendicular diffusion at low density
  bool hyper_z;  // Hyper-diffusion in Z

  Field3D Sn; ///< Density source
};

namespace {
RegisterComponent<EvolveDensity> registercomponentevolvedensity("evolve_density");
}


#endif // EVOLVE_DENSITY_H
