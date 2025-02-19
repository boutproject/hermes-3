
#pragma once
#ifndef NEUTRAL_MIXED_H
#define NEUTRAL_MIXED_H

#include <memory>
#include <string>

#include <bout/invert_laplace.hxx>

#include "component.hxx"

/// Evolve density, parallel momentum and pressure
/// for a neutral gas species with cross-field diffusion
struct NeutralMixed : public Component {
  ///
  /// @param name     The name of the species e.g. "h"
  /// @param options  Top-level options. Settings will be taken from options[name]
  /// @param solver   Time-integration solver to be used
  NeutralMixed(const std::string& name, Options& options, Solver *solver);
  
  /// Modify the given simulation state
  void transform(Options &state) override;
  
  /// Use the final simulation state to update internal state
  /// (e.g. time derivatives)
  void finally(const Options &state) override;

  /// Add extra fields for output, or set attributes e.g docstrings
  void outputVars(Options &state) override;

  /// Preconditioner
  void precon(const Options &state, BoutReal gamma) override;
private:
  std::string name;  ///< Species name
  
  Field3D Nn, Pn, NVn; // Density, pressure and parallel momentum
  Field3D Vn; ///< Neutral parallel velocity
  Field3D Tn; ///< Neutral temperature
  Field3D Nnlim, Pnlim, logPnlim, Vnlim, Tnlim; // Limited in regions of low density

  BoutReal AA; ///< Atomic mass (proton = 1)

  Field3D Dnn; ///< Diffusion coefficient
  Field3D DnnNn, DnnPn, DnnNVn;

  bool sheath_ydown, sheath_yup;

  BoutReal density_floor; ///< Minimum Nn used when dividing NVn by Nn to get Vn.
  BoutReal pressure_floor; ///< Minimum Pn used when dividing Pn by Nn to get Tn.

  BoutReal flux_limit; ///< Diffusive flux limit
  BoutReal diffusion_limit;    ///< Maximum diffusion coefficient

  bool neutral_viscosity; ///< include viscosity?
  bool neutral_conduction; ///< Include heat conduction?
  bool evolve_momentum; ///< Evolve parallel momentum?

  bool precondition {true}; ///< Enable preconditioner?
  bool lax_flux; ///< Use Lax flux for advection terms
  std::unique_ptr<Laplacian> inv; ///< Laplacian inversion used for preconditioning

  Field3D density_source, pressure_source; ///< External input source
  Field3D Sn, Sp, Snv; ///< Particle, pressure and momentum source
  Field3D sound_speed; ///< Sound speed for use with Lax flux
  Field3D perp_nn_adv_src; ///< Source due to perpendicular advection operator
  Field3D par_nn_adv_src; ///< Source due to parallel advection operator

  bool output_ddt; ///< Save time derivatives?
  bool diagnose; ///< Save additional diagnostics?

  // Flow diagnostics
  Field3D pf_adv_perp_xlow, pf_adv_perp_ylow, pf_adv_par_ylow;
  Field3D mf_adv_perp_xlow, mf_adv_perp_ylow, mf_adv_par_ylow;
  Field3D mf_visc_perp_xlow, mf_visc_perp_ylow, mf_visc_par_ylow;
  Field3D ef_adv_perp_xlow, ef_adv_perp_ylow, ef_adv_par_ylow;
  Field3D ef_cond_perp_xlow, ef_cond_perp_ylow, ef_cond_par_ylow;
};

namespace {
RegisterComponent<NeutralMixed> registersolverneutralmixed("neutral_mixed");
}

#endif // NEUTRAL_MIXED_H
