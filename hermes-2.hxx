/*
    Copyright B.Dudson, J.Leddy, University of York, September 2016
              email: benjamin.dudson@york.ac.uk

    This file is part of Hermes-2 (Hot ion).

    Hermes is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Hermes is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Hermes.  If not, see <http://www.gnu.org/licenses/>.
  
*/

class Hermes;

#ifndef __HERMES_H__
#define __HERMES_H__

#include <bout/physicsmodel.hxx>

#include <invert_laplace.hxx>
#include <bout/invert/laplace3d.hxx>
#include <bout/invert/laplacexy.hxx>
#include <bout/invert/laplacexz.hxx>
#include <bout/constants.hxx>

#include "neutral-model.hxx"

class Hermes : public PhysicsModel {
public:
  virtual ~Hermes() {}
protected:
  int init(bool restarting);
  int rhs(BoutReal t);
  
  int precon(BoutReal t, BoutReal gamma, BoutReal delta);
private:
  // Equilibrium current
  Field2D Jpar0;

  // Evolving variables
  Field3D Ne;         // Electron density
  Field3D Pe, Pi;     // Electron and Ion pressures
  Field3D VePsi;      // Combination of Ve and psi
  Field3D Vort;       // Vorticity
  Field3D NVi;        // Parallel momentum

  FieldGroup EvolvingVars;

  // Auxilliary variables
  Field3D Te;         // Electron temperature
  Field3D Ti;         // Ion temperature
  Field3D Ve, Vi, Jpar;  // Electron and ion parallel velocities
  Field3D psi;        // Electromagnetic potential (-A_||)
  Field3D phi;        // Electrostatic potential
  
  // Limited variables
  Field3D Telim, Tilim;

  // Collisional terms
  Field3D nu, kappa_epar, kappa_ipar, Dn;
  BoutReal tau_e0, tau_i0;
  Field3D tau_e, tau_i;          // Collision times for electrons and ions
  Field3D Wi;                    // Energy transfer from electrons to ions
  Field3D Pi_ciperp, Pi_cipar, Pi_ci;   // Ion collisional stress tensor
  BoutReal resistivity_multiply; ///< Factor in front of nu

  BoutReal flux_limit_alpha;  // Flux limiter. < 0 disables
  BoutReal kappa_limit_alpha; // Heat flux limiter from SOLPS
  BoutReal eta_limit_alpha;   // Momentum flux limiter from SOLPS
  
  // Neutral gas model
  NeutralModel *neutrals; // Handles evolution of neutral gas
  bool neutral_friction;
  BoutReal frecycle;  // Recycling fraction
  BoutReal ion_neutral_rate; // Fixed ion-neutral collision rate
  
  // Impurity radiation
  BoutReal carbon_fraction;
  Field3D Rzrad;             // Radiated power
  RadiatedPower *carbon_rad; // Carbon cooling curve
  
  // Switches
  bool evolve_plasma;   // Should plasma be evolved?
  
  bool electromagnetic; // Include magnetic potential psi
  bool FiniteElMass;    // Finite Electron Mass
  
  bool j_diamag;    // Diamagnetic current: Vort <-> Pe
  bool j_par;       // Parallel current:    Vort <-> Psi
  bool parallel_flow;
  bool parallel_flow_p_term; // Vi advection terms in Pe, Pi
  bool pe_par;      // Parallel pressure gradient: Pe <-> Psi
  bool pe_par_p_term; // Includes terms in Pe,Pi equations
  bool resistivity; // Resistivity: Psi -> Pe
  bool thermal_force; // Force due to temperature gradients
  bool electron_viscosity; // Electron parallel viscosity
  bool ion_viscosity;   // Ion viscosity
  bool electron_neutral;   // Include electron-neutral collisions in resistivity
  bool ion_neutral;        // Include ion-neutral collisions in ion collision time
  bool poloidal_flows;  // Include y derivatives in diamagnetic and ExB drifts
  bool thermal_flux;    // Include parallel and perpendicular energy flux from Te gradients
  bool thermal_conduction; // Braginskii electron heat conduction
  bool electron_ion_transfer; // Electron-ion heat transfer
  bool classical_diffusion; // Collisional diffusion, including viscosity
  
  // Anomalous perpendicular diffusion coefficients
  BoutReal anomalous_D;    // Density diffusion
  BoutReal anomalous_chi;  // Electron thermal diffusion
  BoutReal anomalous_nu;   // Momentum diffusion (kinematic viscosity)

  bool anomalous_D_nvi; // Include terms in momentum equation
  bool anomalous_D_pepi; // Include terms in Pe, Pi equations
  
  bool ion_velocity;  // Include Vi terms

  bool phi3d;         // Use a 3D solver for phi
  
  bool staggered;     // Use staggered differencing along B

  bool boussinesq;     // Use a fixed density (Nnorm) in the vorticity equation

  bool sinks; // Sink terms for running 2D drift-plane simulations
  bool sheath_closure; // Sheath closure sink on vorticity (if sinks = true)
  bool drift_wave;     // Drift-wave closure (if sinks=true)

  bool radial_buffers; // Radial buffer regions
  int radial_inner_width; // Number of points in the inner radial buffer
  int radial_outer_width; // Number of points in the outer radial buffer
  BoutReal radial_buffer_D; // Diffusion in buffer region

  Field2D sink_invlpar; // Parallel inverse connection length (1/L_{||}) for
                        // sink terms
  Field2D alpha_dw;

  // Sheath heat transmission factor
  int sheath_model;     // Sets boundary condition model
  BoutReal sheath_gamma_e, sheath_gamma_i;  // Heat transmission
  BoutReal neutral_vwall; // Scale velocity at the wall
  bool sheath_yup, sheath_ydown; 
  bool test_boundaries;

  Field2D wall_flux; // Particle flux to wall (diagnostic)
  Field2D wall_power; // Power flux to wall (diagnostic)
  
  // Fix density in SOL
  bool sol_fix_profiles;
  std::shared_ptr<FieldGenerator> sol_ne, sol_te; // Generating functions
  
  // Output switches for additional information
  bool verbose;    // Outputs additional fields, mainly for debugging
  bool output_ddt; // Output time derivatives
  
  // Numerical dissipation

  BoutReal numdiff, hyper, hyperpar; ///< Numerical dissipation
  int low_pass_z; // Fourier filter in Z 
  BoutReal z_hyper_viscos, x_hyper_viscos, y_hyper_viscos; // 4th-order derivatives
  bool low_n_diffuse; // Diffusion in parallel direction at low density
  bool low_n_diffuse_perp; // Diffusion in perpendicular direction at low density
  BoutReal ne_hyper_z, pe_hyper_z; // Hyper-diffusion
  BoutReal scale_num_cs; // Scale numerical sound speed
  BoutReal floor_num_cs; // Apply a floor to the numerical sound speed
  bool vepsi_dissipation; // Dissipation term in VePsi equation
  bool vort_dissipation; // Dissipation term in Vorticity equation
  
  // Sources and profiles
  
  bool ramp_mesh;   // Use Ne,Pe in the grid file for starting ramp target
  BoutReal ramp_timescale; // Length of time for the initial ramp
  Field2D NeTarget, PeTarget, PiTarget; // For adaptive sources
  
  bool adapt_source; // Use a PI controller to feedback profiles
  bool core_sources; // Sources only in the core
  bool energy_source; // Add the same amount of energy to each particle
  BoutReal source_p, source_i;  // Proportional-Integral controller
  Field2D Sn, Spe, Spi; // Sources in density, Pe and Pi
  Field3D NeSource, PeSource, PiSource; // Actual sources added
  bool density_inflow;  // Does incoming density have momentum?
  
  bool source_vary_g11; // Multiply source by g11
  Field2D g11norm;
  
  // Boundary fluxes
  
  bool pe_bndry_flux;   // Allow flux of pe through radial boundaries
  bool ne_bndry_flux;   // Allow flux of ne through radial boundaries
  bool vort_bndry_flux; // Allow flux of vorticity through radial boundaries
  
  // Normalisation parameters
  BoutReal Tnorm, Nnorm, Bnorm;
  BoutReal AA, Cs0, rho_s0, Omega_ci;
  BoutReal mi_me, beta_e;
  
  // Curvature, Grad-B drift
  Vector3D Curlb_B; // Curl(b/B)
  
  // Perturbed parallel gradient operators
  const Field3D Grad_parP(const Field3D &f);
  const Field3D Div_parP(const Field3D &f);
  
  // Electromagnetic solver for finite electron mass case
  bool split_n0_psi;   // Split the n=0 component of Apar (psi)?
  //Laplacian *aparSolver;
  LaplaceXZ *aparSolver;
  LaplaceXY *aparXY;    // Solves n=0 component
  Field2D psi2D;        // Axisymmetric Psi
  
  // Solvers for the electrostatic potential

  bool split_n0;        // Split solve into n=0 and n~=0?
  LaplaceXY *laplacexy; // Laplacian solver in X-Y (n=0)
  Field2D phi2D;        // Axisymmetric phi
  
  bool newXZsolver; 
  Laplacian *phiSolver; // Old Laplacian in X-Z
  LaplaceXZ *newSolver; // New Laplacian in X-Z

  // Mesh quantities
  Field2D B32, sqrtB;
};

/// Fundamental constants

const BoutReal e0  = 8.854e-12;      // Permittivity of free space
const BoutReal mu0 = 4.e-7*PI;       // Permeability of free space
const BoutReal qe  = 1.602e-19;      // Electron charge
const BoutReal Me  = 9.109e-31;      // Electron mass
const BoutReal Mp  = 1.67262158e-27; // Proton mass

#endif // __HERMES_H__
