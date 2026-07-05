#pragma once
#ifndef CLASSICAL_DIFFUSION_H
#define CLASSICAL_DIFFUSION_H

#include "component.hxx"

struct ClassicalDiffusion : public Component {
  ClassicalDiffusion(std::string name, Options& alloptions, Solver*);

  void outputVars(Options &state) override;
private:
  Field2D Bsq; // Magnetic field squared

  bool diagnose; ///< Output additional diagnostics?
  Field3D Dn; ///< Particle diffusion coefficient
  BoutReal custom_D; ///< User-set particle diffusion coefficient override

  // Flow diagnostics
  Field3D cls_pf_perp_xlow, cls_pf_perp_ylow;
  Field3D cls_mf_perp_xlow, cls_mf_perp_ylow;
  Field3D cls_nef_perp_xlow, cls_nef_perp_ylow;
  Field3D cls_tef_perp_xlow, cls_tef_perp_ylow;
  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<ClassicalDiffusion> registercomponentclassicaldiffusion("classical_diffusion");
}

#endif // CLASSICAL_DIFFUSION_H
