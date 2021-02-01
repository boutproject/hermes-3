
#pragma once
#ifndef NEUTRAL_MIXED_H
#define NEUTRAL_MIXED_H

#include <memory>
#include <string>

#include <invert_laplace.hxx>

#include "component.hxx"

struct NeutralMixed : public Component {
  NeutralMixed(const std::string& name, Options& options, Solver *solver);
  
  /// Modify the given simulation state
  void transform(Options &state);
  
  /// Use the final simulation state to update internal state
  /// (e.g. time derivatives)
  void finally(const Options &state) override;

  /// Add extra fields for output, or set attributes e.g docstrings
  void annotate(Options &state) override;

  /// Preconditioner
  void precon(const Options &state, BoutReal gamma) override;
private:
  std::string name;
  
  Field3D Nn, Pn, NVn; // Density, pressure and parallel momentum
  Field3D Vn; // Neutral parallel velocity
  Field3D Tn; // Neutral temperature
  Field3D Nnlim, Pnlim, Vnlim; // Limited in regions of low density

  Field3D Dnn; // Diffusion coefficient
  
  bool sheath_ydown, sheath_yup;
  
  BoutReal neutral_gamma; // Heat transmission for neutrals
  
  BoutReal nn_floor; ///< Minimum Nn used when dividing NVn by Nn to get Vn.
  
  bool precondition {true}; // Enable preconditioner?
  std::unique_ptr<Laplacian> inv; // Laplacian inversion used for preconditioning

  Field3D Sn, Sp, Snv; ///< Particle, pressure and momentum source
};

namespace {
RegisterComponent<NeutralMixed> registersolverneutralmixed("neutral_mixed");
}

#endif // NEUTRAL_MIXED_H
