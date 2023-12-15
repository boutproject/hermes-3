#pragma once
#ifndef VORTICITY_H
#define VORTICITY_H

#include <bout/vector2d.hxx>

#include "component.hxx"

class LaplaceXY;
class Laplacian;

/// Evolve electron density in time
/// 
struct Vorticity : public Component {
  /// Options
  ///
  /// - <name>
  ///   - average_atomic_mass: float, default 2.0
  ///     Weighted average ion atomic mass for polarisation current
  ///   - bndry_flux: bool, default true
  ///     Allow flows through radial (X) boundaries?
  ///   - collisional_friction: bool, default false
  ///     Damp vorticity based on mass-weighted collision frequency?
  ///   - diagnose: bool, false
  ///     Output additional diagnostics?
  ///   - diamagnetic: bool, default true
  ///     Include diamagnetic current, using curvature vector?
  ///   - diamagnetic_polarisation: bool, default true
  ///     Include ion diamagnetic drift in polarisation current?
  ///   - exb_advection: bool, default true
  ///     Include ExB advection (nonlinear term)?
  ///   - hyper_z: float, default -1.0
  ///     Hyper-viscosity in Z. < 0 means off
  ///   - laplacian: subsection
  ///     Options for the Laplacian phi solver
  ///   - phi_boundary_relax: bool, default false
  ///     Relax radial phi boundaries towards zero-gradient?
  ///   - phi_boundary_timescale: float, 1e-4
  ///     Timescale for phi boundary relaxation [seconds]
  ///   - phi_dissipation: bool, default true
  ///     Parallel dissipation of potential (Recommended)
  ///   - poloidal_flows: bool, default true
  ///     Include poloidal ExB flow?
  ///   - sheath_boundary: bool, default false
  ///     If phi_boundary_relax is false, set the radial boundary to the sheath potential?
  ///   - split_n0: bool, default false
  ///     Split phi into n=0 and n!=0 components?
  ///   - viscosity: Field2D, default 0.0
  ///     Kinematic viscosity [m^2/s]
  ///   - vort_dissipation: bool, default false
  ///     Parallel dissipation of vorticity?
  ///   - damp_core_vorticity: bool, default false
  ///     Damp axisymmetric component of vorticity in cell next to core boundary
  ///
  Vorticity(std::string name, Options &options, Solver *solver);

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

  void outputVars(Options &state) override;

  // Save and restore potential phi
  void restartVars(Options& state) override {
    AUTO_TRACE();

    // NOTE: This is a hack because we know that the loaded restart file
    //       is passed into restartVars in PhysicsModel::postInit
    // The restart value should be used in init() rather than here
    static bool first = true;
    if (first and state.isSet("phi")) {
      first = false;
      phi = state["phi"].as<Field3D>();
    }

    // Save the potential
    set_with_attrs(state["phi"], phi,
                   {{"long_name", "plasma potential"},
                    {"source", "vorticity"}});
  }
private:
  Field3D Vort; // Evolving vorticity

  Field3D phi; // Electrostatic potential
  std::unique_ptr<Laplacian> phiSolver; // Laplacian solver in X-Z

  Field3D Pi_hat; ///< Contribution from ion pressure, weighted by atomic mass / charge

  bool exb_advection; // Include nonlinear ExB advection?
  bool exb_advection_simplified; // Simplify nonlinear ExB advection form?
  bool diamagnetic; // Include diamagnetic current?
  bool diamagnetic_polarisation; // Include diamagnetic drift in polarisation current
  BoutReal average_atomic_mass; // Weighted average atomic mass, for polarisaion current (Boussinesq approximation)
  bool poloidal_flows;   ///< Include poloidal ExB flow?
  bool bndry_flux;  ///< Allow flows through radial boundaries?

  bool collisional_friction; ///< Damping of vorticity due to collisional friction

  bool sheath_boundary; ///< Set outer boundary to j=0?

  bool vort_dissipation; ///< Parallel dissipation of vorticity
  bool phi_dissipation;  ///< Parallel dissipation of potential
  bool phi_sheath_dissipation; ///< Dissipation at the sheath if phi < 0
  bool damp_core_vorticity; ///< Damp axisymmetric component of vorticity

  bool phi_boundary_relax; ///< Relax boundary to zero-gradient
  BoutReal phi_boundary_timescale; ///< Relaxation timescale [normalised]
  BoutReal phi_boundary_last_update; ///< Time when last updated

  bool split_n0; // Split phi into n=0 and n!=0 components
  LaplaceXY* laplacexy; // Laplacian solver in X-Y (n=0)

  Field2D Bsq; // SQ(coord->Bxy)
  Vector2D Curlb_B; // Curvature vector Curl(b/B)
  BoutReal hyper_z; ///< Hyper-viscosity in Z
  Field2D viscosity; ///< Kinematic viscosity

  // Diagnostic outputs
  Field3D DivJdia, DivJcol; // Divergence of diamagnetic and collisional current

  bool diagnose; ///< Output additional diagnostics?
};

namespace {
RegisterComponent<Vorticity> registercomponentvorticity("vorticity");
}

#endif // VORTICITY_H
