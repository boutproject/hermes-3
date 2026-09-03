#pragma once
#ifndef KMODEL_H
#define KMODEL_H

#include "component.hxx"

/// Evolve the turbulent kinetic energy with the k-model ansatz
/// from H. Bufferand et al 2021 Nucl. Fusion 61 116052
/// DOI: 10.1088/1741-4326/ac2873

struct Kmodel : public NamedComponent<Kmodel> {

  static constexpr auto type = "kmodel";

  Kmodel(std::string name, Options& options, Solver* solver);

  void finally(const Options& state) override;

  void outputVars(Options& state) override;

private:
  Field3D k;

  Field3D D_k;

  Field3D chi_k, nu_k;

  Field3D alpha;

  Field3D gradPgradB_X, gradPgradB_Y, gradPgradB_Z;

  Field3D DDX_P, DDX_B;

  Field3D DDY_P, DDY_B;

  Field3D DDZ_P, DDZ_B;

  Field3D gamma;

  Field3D S_k, P_k;

  Field3D klim;

  Field3D Bxy;

  Field3D Pi_hat, Ni_hat;

  BoutReal R_major, R_minor;

  BoutReal L_par;
  BoutReal lambda_q;
  BoutReal average_AA;
  Field3D SQ_lambda_SOL;

  Field3D flux_k_x, flux_k_y;

  BoutReal chi_factor, nu_factor;

  bool diagnose;

  bool diffusion, advection;

  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<Kmodel> registercomponentkmodel;
}

#endif // KMODEL_H
