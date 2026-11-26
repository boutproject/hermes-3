#pragma once
#ifndef CLASSICAL_DIFFUSION_H
#define CLASSICAL_DIFFUSION_H

#include <bout/coordinates.hxx>

#include "component.hxx"
#include "domain_balance.hxx"

/// Perpendicular diffusion of particles, carried parallel momentum and energy.
///
/// A positive scalar custom_D sets the particle and thermal-energy diffusivity
/// in m^2/s, converted to normalised units. The energy source is
/// Div(1.5 D Grad_perp P), using the species
/// pressure or N*T when only temperature is supplied. No collision frequencies
/// are required in this mode.
///
/// Negative custom_D (default -1) calculates the ambipolar particle diffusivity
/// from electron collisions, transports 1.5 T with the particles, and adds
/// perpendicular temperature conduction for each charged species. A value of
/// zero calculates the same particle diffusivity but omits conduction.
///
/// All physical X and Y boundary face fluxes are closed, including mixed-metric
/// contributions. Periodic and inter-processor faces are retained.
/// Optional low_density_enhance multiplies Dn smoothly below an electron
/// density floor supplied in m^-3, with a bounded enhancement multiplier.
struct ClassicalDiffusion : public NamedComponent<ClassicalDiffusion> {
  /// Options: custom_D selects the transport model; diagnose outputs the
  /// normalised particle diffusivity as D_classical with SI conversion metadata.
  ClassicalDiffusion(std::string name, Options& alloptions, Solver*);

  void outputVars(Options& state) override;

  static constexpr auto type = "classical_diffusion";

private:
  DomainBalance balance;
  Coordinates::FieldMetric Bsq; // Magnetic field squared

  bool diagnose;     ///< Output additional diagnostics?
  Field3D Dn;        ///< Realised normalised perpendicular particle diffusivity
  BoutReal custom_D; ///< Normalised diffusivity override; positive selects custom mode
  bool low_density_enhance;        ///< Enhance Dn below the electron density floor?
  BoutReal low_density_floor{0.0}; ///< Normalised electron density for enhancement onset
  BoutReal low_density_max_multiplier{1.0}; ///< Upper bound on the enhancement

  // Flow diagnostics
  Field3D cls_pf_perp_xlow, cls_pf_perp_ylow;
  Field3D cls_mf_perp_xlow, cls_mf_perp_ylow;
  Field3D cls_nef_perp_xlow, cls_nef_perp_ylow;
  Field3D cls_tef_perp_xlow, cls_tef_perp_ylow;
  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<ClassicalDiffusion> registercomponentclassicaldiffusion;
}

#endif // CLASSICAL_DIFFUSION_H
