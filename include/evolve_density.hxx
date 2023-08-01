#pragma once
#ifndef EVOLVE_DENSITY_H
#define EVOLVE_DENSITY_H

#include "component.hxx"

/// Evolve species density in time
///
/// # Mesh inputs
///
/// N<name>_src  A source of particles, per cubic meter per second.
///              This can be over-ridden by the `source` option setting.
struct EvolveDensity : public Component {
  /// # Inputs
  ///
  /// - <name>
  ///   - charge         Particle charge e.g. hydrogen = 1
  ///   - AA             Atomic mass number e.g. hydrogen = 1
  ///   - bndry_flux     Allow flow through radial boundaries? Default is true.
  ///   - poloidal_flows Include poloidal ExB flows? Default is true.
  ///   - density_floor  Minimum density floor. Default is 1e-5 normalised units
  ///   - low_n_diffuse  Enhance parallel diffusion at low density? Default false
  ///   - hyper_z        Hyper-diffusion in Z. Default off.
  ///   - evolve_log     Evolve logarithm of density? Default false.
  ///   - diagnose       Output additional diagnostics?
  ///
  /// - N<name>   e.g. "Ne", "Nd+"
  ///   - source    Source of particles [/m^3/s]
  ///               NOTE: This overrides mesh input N<name>_src
  ///   - source_only_in_core         Zero the source outside the closed field-line region?
  ///   - neumann_boundary_average_z  Apply Neumann boundaries with Z average?
  ///
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

  void outputVars(Options &state) override;
private:
  std::string name;     ///< Short name of species e.g "e"

  BoutReal charge;      ///< Species charge e.g. electron = -1
  BoutReal AA;          ///< Atomic mass e.g. proton = 1
  
  Field3D N;            ///< Species density (normalised, evolving)

  bool bndry_flux;      ///< Allow flows through boundaries?
  bool poloidal_flows;  ///< Include ExB flow in Y direction?
  bool neumann_boundary_average_z; ///< Apply neumann boundary with Z average?

  BoutReal density_floor;
  bool low_n_diffuse;   ///< Parallel diffusion at low density
  bool low_n_diffuse_perp;  ///< Perpendicular diffusion at low density
  BoutReal hyper_z;    ///< Hyper-diffusion in Z

  bool evolve_log; ///< Evolve logarithm of density?
  Field3D logN;    ///< Logarithm of density (if evolving)

  Field3D source; ///< External input source
  Field3D Sn; ///< Total density source

  bool diagnose; ///< Output additional diagnostics?
  Field3D flow_xlow, flow_ylow; ///< Particle flow diagnostics
};

namespace {
RegisterComponent<EvolveDensity> registercomponentevolvedensity("evolve_density");
}


#endif // EVOLVE_DENSITY_H
