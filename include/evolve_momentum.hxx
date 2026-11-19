#pragma once
#ifndef EVOLVE_MOMENTUM_H
#define EVOLVE_MOMENTUM_H

#include "component.hxx"

#include <bout/invert_pardiv.hxx>

#include <memory>

/// Evolve parallel momentum
struct EvolveMomentum : public NamedComponent<EvolveMomentum> {
  EvolveMomentum(std::string name, Options& options, Solver* solver);

  /// Calculate ddt(NV).
  ///
  /// Inputs
  /// - species
  ///   - <name>
  ///     - density
  ///     - velocity
  ///     - pressure (optional)
  ///     - momentum_source (optional)
  ///     - sound_speed (optional, used for numerical dissipation)
  ///     - temperature (only needed if sound_speed not provided)
  /// - fields
  ///   - phi (optional)
  void finally(const Options& state) override;

  void outputVars(Options& state) override;

  /// Parallel viscosity preconditioner, used when braginskii_ion_viscosity or
  /// braginskii_electron_viscosity sets `precondition` for this species. Inverts the
  /// viscous force linearised at frozen eta and density, with the sheath-face
  /// closure at the targets.
  void precon(const Options& state, BoutReal gamma) override;

  static constexpr auto type = "evolve_momentum";

private:
  std::string name; ///< Short name of species e.g "e"

  Field3D NV;        ///< Species parallel momentum (normalised, evolving)
  Field3D NV_err;    ///< Difference from momentum as input from solver
  Field3D NV_solver; ///< Momentum as calculated in the solver
  Field3D V;         ///< Species parallel velocity

  Field3D momentum_source; ///< From other components. Stored for diagnostic output

  bool bndry_flux;     // Allow flows through boundaries?
  bool poloidal_flows; // Include ExB flow in Y direction?

  BoutReal density_floor;
  bool low_n_diffuse_perp; ///< Cross-field diffusion at low density?
  BoutReal pressure_floor;
  bool low_p_diffuse_perp;    ///< Cross-field diffusion at low pressure?
  bool use_div_par_fvv{true}; ///< Use Div_par_fvv instead of Div_par

  BoutReal hyper_z; ///< Hyper-diffusion

  bool diagnose;                   ///< Output additional diagnostics?
  bool fix_momentum_boundary_flux; ///< Fix momentum flux to boundary condition?

  std::unique_ptr<InvertParDiv> viscosity_inverse;

  /// Per target guard cell (field-aligned): -1 where the RHS fixes V at the sheath face
  /// (V_guard = 2 v_sheath - V), +1 where it extrapolates with zero gradient
  Field3D sheathClosure(const Field3D& V) const;
  Field3D flow_xlow, flow_ylow; ///< Momentum flow diagnostics

  /// This takes as inputs
  /// - species
  ///   - <name>
  ///     - AA
  ///     - density
  ///
  /// This sets in the state
  /// - species
  ///   - <name>
  ///     - momentum
  ///     - velocity  if density is defined
  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<EvolveMomentum> registercomponentevolvemomentum;
}

#endif // EVOLVE_MOMENTUM_H
