#pragma once
#ifndef CLASSICAL_DIFFUSION_H
#define CLASSICAL_DIFFUSION_H

#include "component.hxx"

struct ClassicalDiffusion : public Component {
  ClassicalDiffusion(std::string name, Options& alloptions, Solver*);

  void outputVars(Options &state) override;
private:
  std::string name; ///< Short name of species e.g "e"
  Field2D Bsq; // Magnetic field squared

  bool diagnose; ///< Output additional diagnostics?
  Field3D Dn; ///< Particle diffusion coefficient
  BoutReal custom_D; ///< User-set particle diffusion coefficient override
  Field3D nu; 
  Field3D Kappa_perp;
  bool nonorthogonal_operators;   ///< Use nonorthogonal operators for radial transport?
  bool zero_BC_transport; ///< Set transport to zero in ghost cells?
  // Flow diagnostics
  Field3D cls_pf_perp_xlow, cls_pf_perp_ylow;
  Field3D cls_mf_perp_xlow, cls_mf_perp_ylow;
  Field3D cls_nef_perp_xlow, cls_nef_perp_ylow;
  Field3D cls_tef_perp_xlow, cls_tef_perp_ylow;
  Field3D cls_energy_flow_xlow;   ///< Cell edge fluxes used for calculating
                              ///<
                              ///< fast recycling energy source
  Field3D cls_particle_flow_xlow; ///< Radial wall particle fluxes for recycling calc. No need
                              ///<
                              ///< to get poloidal from here, it's calculated from sheath
                              ///<
                              ///< velocity
  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<ClassicalDiffusion> registercomponentclassicaldiffusion("classical_diffusion");
}

#endif // CLASSICAL_DIFFUSION_H
