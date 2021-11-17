#pragma once
#ifndef RELAX_POTENTIAL_H
#define RELAX_POTENTIAL_H

#include <vector2d.hxx>

#include "component.hxx"

/// Evolve vorticity and potential in time.
///
/// Uses a relaxation method for the potential, which is valid for
/// steady state, but not for timescales shorter than the relaxation
/// timescale.
///
struct RelaxPotential : public Component {
  /// Options
  ///
  /// - <name>
  ///   - diamagnetic
  ///   - diamagnetic_polarisation
  ///   - average_atomic_mass
  ///   - bndry_flux
  ///   - poloidal_flows
  ///   - split_n0
  ///   - laplacian
  ///     Options for the Laplacian phi solver
  ///
  RelaxPotential(std::string name, Options &options, Solver *solver);

  /// Optional inputs
  ///
  /// - species
  ///   - pressure and charge => Calculates diamagnetic terms [if diamagnetic=true]
  ///   - pressure, charge and mass => Calculates polarisation current terms [if diamagnetic_polarisation=true]
  ///
  /// Sets in the state
  /// - species
  ///   - [if has pressure and charge]
  ///     - energy_source
  /// - fields
  ///   - vorticity
  ///   - phi         Electrostatic potential
  ///   - DivJdia     Divergence of diamagnetic current [if diamagnetic=true]
  ///
  /// Note: Diamagnetic current calculated here, but could be moved
  ///       to a component with the diamagnetic drift advection terms
  void transform(Options &state) override;

  /// Optional inputs
  /// - fields
  ///   - DivJextra    Divergence of current, including parallel current
  ///                  Not including diamagnetic or polarisation currents
  ///
  void finally(const Options &state) override;

private:
  Field3D Vort; // Evolving vorticity

  Field3D phi1; // Scaled electrostatic potential, evolving in time ϕ_1 = λ_2 ϕ
  Field3D phi;  // Electrostatic potential

  bool exb_advection; // Include nonlinear ExB advection?
  bool diamagnetic; // Include diamagnetic current?
  bool diamagnetic_polarisation; // Include diamagnetic drift in polarisation current
  BoutReal average_atomic_mass; // Weighted average atomic mass, for polarisaion current (Boussinesq approximation)
  bool poloidal_flows;   ///< Include poloidal ExB flow?
  bool bndry_flux;  ///< Allow flows through radial boundaries?
  bool boussinesq;  ///< Use the Boussinesq approximation?

  bool sheath_boundary; ///< Set outer boundary to j=0?

  Field2D Bsq; // SQ(coord->Bxy)
  Vector2D Curlb_B; // Curvature vector Curl(b/B)

  BoutReal lambda_1, lambda_2;
};

namespace {
RegisterComponent<RelaxPotential> registercomponentrelaxpotential("relax_potential");
}

#endif // RELAX_POTENTIAL_H
