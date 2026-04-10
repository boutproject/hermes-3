
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
  NeutralMixed(const std::string& name, Options& options, Solver* solver);

  /// Use the final simulation state to update internal state
  /// (e.g. time derivatives)
  void finally(const Options& state) override;

  /// Add extra fields for output, or set attributes e.g docstrings
  void outputVars(Options& state) override;

  /// Preconditioner
  void precon(const Options& state, BoutReal gamma) override;

private:
  std::string name; ///< Species name

  Field3D Nn, Pn, NVn;            // Density, pressure and parallel momentum
  Field3D Vn;                     ///< Neutral parallel velocity
  Field3D Tn;                     ///< Neutral temperature
  Field3D Nnlim, Pnlim, logPnlim; // Limited in regions of low density

  Field3D NVn_err;    ///< Difference from momentum as input from solver
  Field3D NVn_solver; ///< Momentum as calculated in the solver

  BoutReal AA; ///< Atomic mass (proton = 1)

  std::vector<std::string> collision_names; ///< Collisions used for collisionality
  std::string
      diffusion_collisions_mode; ///< Collision selection, either afn or multispecies
  Field3D Vnth_pf = 0.0;
  Field3D Vnth_hf = 0.0;

  Field3D nu;                    ///< Collisionality to use for diffusion
  Field3D nu_pseudo_mfp;         ///< Pseudo-collision frequency based on mean free path
  Field3D nu_total; ///< Total collision frequency used for diffusion, including
                    ///< pseudo-collisions
  Field3D Dnn, Dnn_unlimited, Dmax;    ///< Diffusion coefficient
  Field3D DnnNn, DnnPn, DnnTn, DnnNVn; ///< Used for operators
  Field3D kappa_n_unlimited, kappa_n_max_par, kappa_n_max_perp;
  BoutReal flux_limit_adv;       ///< Diffusive flux limit
  BoutReal flux_limit_cond_par;  ///< Limit for parallel conductive flux
  BoutReal flux_limit_cond_perp; ///< Limit for perpendicular conductive flux
  BoutReal flux_limit_visc_par;  ///< Limit for parallel viscous flux
  BoutReal flux_limit_visc_perp; ///< Limit for perpendicular viscous flux

  BoutReal diffusion_limit; ///< Maximum diffusion coefficient
  BoutReal neutral_lmax;
  BoutReal flux_limiter_sharpness; ///< Sharpness of flux limiter transition

  bool sheath_ydown, sheath_yup;

  BoutReal density_floor; ///< Minimum Nn used when dividing NVn by Nn to get Vn.
  BoutReal temperature_floor;
  BoutReal pressure_floor; ///< Minimum Pn used when dividing Pn by Nn to get Tn.
  bool freeze_low_density; ///< Freeze evolution in low density regions?

  BoutReal collisionality_override;     ///< Rnn input for testing
  BoutReal density_norm, pressure_norm; ///< Normalisations
  BoutReal momentum_norm;               ///< Normalisations

  bool neutral_viscosity;  ///< include viscosity?
  bool neutral_conduction; ///< Include heat conduction?
  bool evolve_momentum;    ///< Evolve parallel momentum?
  bool normalise_sources;  ///< Normalise input sources?
  bool perp_ion_coupling;  ///< Include coupling to ion perpendicular velocity?

  // Temporary variables
  Field3D debug;          ///< Debug variable FIXME: remove
  bool double_count_lmax; ///< Include neutral_lmax in Dmax and kappa_max as well as Dnn?
  bool legacy_thermal_speed; ///< Use legacy definition of thermal speed in flux limiter?
  bool legacy_limiter_form;  ///< Use legacy form of flux limiter rather than SOLPS-style

  Field3D kappa_n, eta_n_unlimited;      ///< Neutral conduction and viscosity
  Field3D kappa_n_perp, eta_n_perp;      ///< Neutral conduction and viscosity
  Field3D kappa_n_par, eta_n_par;        ///< Neutral conduction and viscosity
  Field3D eta_n_max_par, eta_n_max_perp; ///< Viscosity reduction factor for flux-limiting

  bool nonorthogonal_operators;   ///< Use nonorthogonal operators for radial transport?
  bool precondition{true};        ///< Enable preconditioner?
  bool precon_laplacexy{false};   ///< Use LaplaceXY?
  bool lax_flux;                  ///< Use Lax flux for advection terms
  std::unique_ptr<Laplacian> inv; ///< Laplacian inversion used for preconditioning

  Field3D density_source, pressure_source, momentum_source; ///< External input source
  Field3D Sn, Sp, Snv; ///< Particle, pressure and momentum source
  Field3D sound_speed; ///< Sound speed for use with Lax flux

  bool zero_timederivs; ///< Set the time derivatives to zero?

  bool output_ddt; ///< Save time derivatives?
  bool diagnose;   ///< Save additional diagnostics?

  //Physic terms
  Field3D ddtN_par_advection = 0.0;
  Field3D ddtN_perp_diffusion = 0.0;

  Field3D ddtPn_par_advection = 0.0;
  Field3D ddtPn_work_done = 0.0;
  Field3D ddtPn_perp_advection = 0.0;
  Field3D ddtPn_par_conduction = 0.0;
  Field3D ddtPn_perp_conduction = 0.0;

  Field3D ddtNVn_par_advection = 0.0;
  Field3D ddtNVn_pressure_gradient = 0.0;
  Field3D ddtNVn_perp_advection = 0.0;

  Field3D par_viscosity_source = 0.0;
  Field3D perp_viscosity_source = 0.0;
  Field3D ddtNVn_viscosity = 0.0; // par_viscosity_source + perp_viscosity_source
  Field3D ddtPn_viscosity = 0.0; 

  Field3D ddtN_anomalous_transport = 0.0;
  Field3D ddtNVn_anomalous_transport = 0.0;
  Field3D ddtPn_anomalous_transport = 0.0;

  // Flow diagnostics
  Field3D pf_adv_perp_xlow, pf_adv_perp_ylow, pf_adv_par_ylow;
  Field3D mf_adv_perp_xlow, mf_adv_perp_ylow, mf_adv_par_ylow;
  Field3D mf_visc_perp_xlow, mf_visc_perp_ylow, mf_visc_par_ylow;
  Field3D ef_adv_perp_xlow, ef_adv_perp_ylow, ef_adv_par_ylow;
  Field3D ef_cond_perp_xlow, ef_cond_perp_ylow, ef_cond_par_ylow;

  /// Sets
  /// - species
  ///   - <name>
  ///     - AA
  ///     - density
  ///     - momentum
  ///     - pressure
  ///     - temperature
  ///     - velocity
  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<NeutralMixed> registersolverneutralmixed("neutral_mixed");
}

#endif // NEUTRAL_MIXED_H
