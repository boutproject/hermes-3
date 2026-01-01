#pragma once
#ifndef KMODEL_H
#define KMODEL_H

#include <bout/vectormetric.hxx>
#include <bout/yboundary_regions.hxx>

#include "component.hxx"
#include "div_ops.hxx"

struct Kmodel : public Component {

  Kmodel(std::string name, Options &options, Solver *solver);


  void transform(Options &state) override;

  void finally(const Options &state) override;

  void outputVars(Options &state) override;

private:
  Field3D k; // Evolving vorticity

  Field3D D_k; // Electrostatic potential

  Field3D chi_k;

  Field3D nu_k;

  Field3D alpha;

  Field3D gradPgradB_X, gradPgradB_Z;
  Field3D DDX_P, DDX_B;
  Field3D DDZ_P, DDZ_B;
  Field3D gamma;

  bool diffusion, advection;
  
  Field3D S_k;

  Field3D P_k;

  bool dissipative;

  bool propagate;
  
  Field3D Bxy;
  
  BoutReal R_major;
  BoutReal R_minor;

  BoutReal L_par;
  BoutReal lambda_q;
  Field3D SQ_lambdda_SOL;

  BoutReal chi_factor, nu_factor;
  
  std::shared_ptr<FCI::dagp_fv> dagp;

  bool diagnose;
  
  Field3D Div_a_Grad_perp(Field3D a, Field3D b) {
    if (k.isFci()) {
      return (*dagp)(a, b, false);
    }
    return FV::Div_a_Grad_perp(a, b);
  }

  Coordinates::FieldMetric bracket_factor; ///< For non-Clebsch coordinate systems (e.g. FCI)
};

namespace {
RegisterComponent<Kmodel> registercomponentkmodel("kmodel");
}

#endif // KMODEL_H
