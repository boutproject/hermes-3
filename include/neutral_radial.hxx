#pragma once
#ifndef NEUTRAL_RADIAL_H
#define NEUTRAL_RADIAL_H

#include <string>
#include <vector>

#include "component.hxx"

/// Neutral gas species evolving only its density, at a fixed temperature, and moving
/// only in the radial (x) direction: no parallel flow. Reduced copy of neutral_mixed:
///
///   ddt(Nn) = (1/J) d/dx (J g11 Dnn dNn/dx) + Sn
///
/// with the neutral_mixed diffusion coefficient Dnn = (Tn/AA) / nu.
struct NeutralRadial : public NamedComponent<NeutralRadial> {
  /// # Inputs
  /// - <name>
  ///   - AA                        Atomic mass. Default 1
  ///   - temperature               Fixed neutral temperature [eV]. Default 3
  ///   - diffusion_collisions_mode afn (CX, IZ, NN; default) or multispecies
  ///   - neutral_lmax              Maximum mean free path [m] (collisionality floor)
  ///   - collisionality_override   Fixed collision frequency (normalised). < 0: off
  ///   - diffusion_limit           Upper limit on Dnn [m^2/s]. < 0: off
  ///   - density_floor             Normalised
  ///   - diagnose                  Save Dnn, Sn and the radial flow
  /// - N<name>
  ///   - source                    Density source [m^-3 s^-1]
  ///   - bndry_*                   Boundary conditions on Nn (default neumann)
  /// - Dnn<name>
  ///   - bndry_*                   Default dirichlet: no diffusive flux through x walls
  NeutralRadial(const std::string& name, Options& options, Solver* solver);

  void finally(const Options& state) override;
  void outputVars(Options& state) override;

  static constexpr auto type = "neutral_radial";

private:
  std::string name; ///< Species name

  Field3D Nn;             ///< Evolved neutral density
  Field3D Dnn;            ///< Radial diffusion coefficient
  Field3D density_source; ///< User source
  Field3D Sn;             ///< Total density source (user + other components)
  Field3D flow_xlow;      ///< Radial particle flow through x faces

  BoutReal AA;             ///< Atomic mass
  BoutReal Tn;             ///< Fixed temperature (normalised)
  BoutReal density_floor;  ///< Normalised
  BoutReal neutral_lmax;   ///< Normalised
  BoutReal collisionality_override;
  BoutReal diffusion_limit; ///< Normalised; < 0: off
  std::string diffusion_collisions_mode;
  std::vector<std::string> collision_names; ///< Collisions used for nu
  bool diagnose;

  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<NeutralRadial> registercomponentneutralradial;
}

#endif // NEUTRAL_RADIAL_H
