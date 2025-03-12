
#pragma once
#ifndef NEUTRAL_MIXED_H
#define NEUTRAL_MIXED_H

#include <memory>
#include <string>

#include <bout/invert_laplace.hxx>

#include "component.hxx"
#include <bout/yboundary_regions.hxx>

extern YBoundary yboundary;

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
  
  std::shared_ptr<FCI::dagp_fv> dagp;
  Field3D Nn, Pn, NVn; // Density, pressure and parallel momentum
  Field3D Vn; ///< Neutral parallel velocity
  Field3D Tn; ///< Neutral temperature
  Field3D Nnlim, Pnlim, logPnlim, Vnlim, Tnlim; // Limited in regions of low density

  BoutReal AA; ///< Atomic mass (proton = 1)

  Field3D Dnn; ///< Diffusion coefficient
  Field3D DnnNn, DnnPn, DnnTn, DnnNVn; ///< Used for operators

  bool sheath_ydown, sheath_yup;

  BoutReal nn_floor; ///< Minimum Nn used when dividing NVn by Nn to get Vn.

  BoutReal flux_limit; ///< Diffusive flux limit
  BoutReal diffusion_limit;    ///< Maximum diffusion coefficient

  bool neutral_viscosity; ///< include viscosity?
  bool evolve_momentum; ///< Evolve parallel momentum?

  bool precondition {true}; ///< Enable preconditioner?
  std::unique_ptr<Laplacian> inv; ///< Laplacian inversion used for preconditioning

  Field3D density_source, pressure_source; ///< External input source
  Field3D Sn, Sp, Snv; ///< Particle, pressure and momentum source

  bool output_ddt; ///< Save time derivatives?
  bool diagnose; ///< Save additional diagnostics?

  Field3D Div_a_Grad_perp(Field3D a, Field3D b) {
    if (a.isFci()) {
      return (*dagp)(a, b, false);
    }
    return FV::Div_a_Grad_perp(a, b);
  }

};

namespace {
RegisterComponent<NeutralMixed> registersolverneutralmixed("neutral_mixed");
}

#endif // NEUTRAL_MIXED_H
