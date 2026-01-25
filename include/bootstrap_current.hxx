#pragma once
#ifndef BOOTSTRAP_CURRENT_H
#define BOOTSTRAP_CURRENT_H

#include "component.hxx"

///
///
///
struct BootstrapCurrent : public Component {
  /// Options
  /// - units
  /// - <name>
  ///   - diagnose   Saves Ajpar and alpha_em time-dependent values
  ///
  BootstrapCurrent(std::string name, Options& options, Solver* solver);

  void outputVars(Options& state) override;

  // Individual equations for unit testing

  /// Equation 14a
  BoutReal F31(BoutReal X, BoutReal Z);

  /// Equation 14b
  BoutReal L31(BoutReal Z, BoutReal f, BoutReal nu_e);

  // Equation 15a
  BoutReal L32(BoutReal Z, BoutReal f, BoutReal nu_e);

  // Equation 16a
  BoutReal L34(BoutReal Z, BoutReal f, BoutReal nu_e);

  // Equation 17b
  BoutReal alpha(BoutReal f, BoutReal nu_i);

  /// Flux surface average
  /// <f> = integral(f / Bp) / integral(1 / Bp)
  Field2D fluxSurfaceAverage(const Field3D& f);
  Field2D fluxSurfaceAverage(const Field2D& f);

  /// Calculate trapped fraction
  /// Note: This could be a very slow calculation
  Field2D trappedFraction(const Field2D& Bxy);

  /// Return maximum value on each flux surface
  Field2D maxY(const Field2D& f);

private:
  Mesh* mesh;
  BoutReal Nnorm, Tnorm, Pnorm, Bnorm; ///< Normalization factors

  bool diagnose;    ///< Output additional diagnostics?
  Field2D JparB_av; ///< Current value of <J||B>
  Field2D JparB_bs; ///< Steady state <J||B>

  Field2D Bxy; ///< Magnetic field strength
  Field2D RBt; ///< I(psi) = R * Btor
  Field2D J;   ///< Jacobian = hthe / Bpxy

  Field2D averageJ; /// Used to calculate flux surface averages

  Field2D trapped_fraction; ///< f_t, constant on each flux surface

  /// Inputs
  /// - species
  ///   - <..>      All species with charge and parallel momentum
  ///     - charge
  ///     - momentum
  ///     - density
  ///     - AA
  ///
  /// Sets
  /// - species
  ///   - <..>      All species with charge and parallel momentum
  ///     - momentum  (modifies) to m n v||
  ///     - velocity  (modifies) to v||
  ///
  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<BootstrapCurrent>
    registercomponentbootstrapcurrent("bootstrap_current");
}

#endif // BOOTSTRAP_CURRENT_H
