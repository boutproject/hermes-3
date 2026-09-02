#pragma once
#ifndef KMODEL_H
#define KMODEL_H

#include "component.hxx"

/// Evolve the turbulent kinetic energy with the k-model ansatz
/// from H. Bufferand et al 2021 Nucl. Fusion 61 116052
/// DOI: 10.1088/1741-4326/ac2873

struct Kmodel : public NamedComponent<Kmodel> {

  Kmodel(std::string name, Options& options, Solver* solver);

  void finally(const Options& state) override;

  void outputVars(Options& state) override;

private:
  Field3D k;

  Field3D D_k;

  Field3D chi_k, nu_k;

  Field3D alpha;

  Field3D gradPgradB_X, gradPgradB_Z;

  Field3D DDX_P, DDX_B;

  Field3D DDZ_P, DDZ_B;

  Field3D gamma;

  Field3D S_k, P_k;

  Field3D Bxy;

  BoutReal R_major, R_minor;

  BoutReal L_par;
  BoutReal lambda_q;

  Field3D SQ_lambda_SOL;

  BoutReal chi_factor, nu_factor;

  bool diagnose;
};

#endif // KMODEL_H
