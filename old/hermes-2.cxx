/*
  
    Copyright B.Dudson, J.Leddy, University of York, 2016-2019
              email: benjamin.dudson@york.ac.uk

    This file is part of Hermes-2 (Hot ion version)

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
#include "hermes-2.hxx"

#include <derivs.hxx>
#include <field_factory.hxx>
#include <initialprofiles.hxx>

#include <invert_parderiv.hxx>

#include "div_ops.hxx"
#include "loadmetric.hxx"

#include <bout/constants.hxx>
#include <bout/assert.hxx>
#include <bout/fv_ops.hxx>


// OpenADAS interface Atomicpp by T.Body
#include "atomicpp/ImpuritySpecies.hxx"
#include "atomicpp/Prad.hxx"

namespace FV {
  template<typename CellEdges = MC>
  const Field3D Div_par_fvv(const Field3D &f_in, const Field3D &v_in,
                            const Field3D &wave_speed_in, bool fixflux=true) {

    ASSERT1(areFieldsCompatible(f_in, v_in));
    ASSERT1(areFieldsCompatible(f_in, wave_speed_in));

    Mesh* mesh = f_in.getMesh();

    CellEdges cellboundary;

    /// Ensure that f, v and wave_speed are field aligned
    Field3D f = toFieldAligned(f_in, "RGN_NOX");
    Field3D v = toFieldAligned(v_in, "RGN_NOX");
    Field3D wave_speed = toFieldAligned(wave_speed_in, "RGN_NOX");

    Coordinates *coord = f_in.getCoordinates();

    Field3D result{zeroFrom(f)};
    
    // Only need one guard cell, so no need to communicate fluxes
    // Instead calculate in guard cells to preserve fluxes
    int ys = mesh->ystart-1;
    int ye = mesh->yend+1;

    for (int i = mesh->xstart; i <= mesh->xend; i++) {

      if (!mesh->firstY(i) || mesh->periodicY(i)) {
        // Calculate in guard cell to get fluxes consistent between processors
        ys = mesh->ystart - 1;
      } else {
        // Don't include the boundary cell. Note that this implies special
        // handling of boundaries later
        ys = mesh->ystart;
      }

      if (!mesh->lastY(i) || mesh->periodicY(i)) {
        // Calculate in guard cells
        ye = mesh->yend + 1;
      } else {
        // Not in boundary cells
        ye = mesh->yend;
      }

      for (int j = ys; j <= ye; j++) {
        // Pre-calculate factors which multiply fluxes

        // For right cell boundaries
        BoutReal common_factor = (coord->J(i, j) + coord->J(i, j + 1)) /
          (sqrt(coord->g_22(i, j)) + sqrt(coord->g_22(i, j + 1)));
        
        BoutReal flux_factor_rc = common_factor / (coord->dy(i, j) * coord->J(i, j));
        BoutReal flux_factor_rp = common_factor / (coord->dy(i, j + 1) * coord->J(i, j + 1));

        // For left cell boundaries
        common_factor = (coord->J(i, j) + coord->J(i, j - 1)) /
          (sqrt(coord->g_22(i, j)) + sqrt(coord->g_22(i, j - 1)));

        BoutReal flux_factor_lc = common_factor / (coord->dy(i, j) * coord->J(i, j));
        BoutReal flux_factor_lm = common_factor / (coord->dy(i, j - 1) * coord->J(i, j - 1));
        
        for (int k = 0; k < mesh->LocalNz; k++) {

          ////////////////////////////////////////////
          // Reconstruct f at the cell faces
          // This calculates s.R and s.L for the Right and Left
          // face values on this cell
          
          // Reconstruct f at the cell faces
          Stencil1D s;
          s.c = f(i, j, k);
          s.m = f(i, j - 1, k);
          s.p = f(i, j + 1, k);
          
          cellboundary(s); // Calculate s.R and s.L

          // Reconstruct v at the cell faces
          Stencil1D sv;
          sv.c = v(i, j, k);
          sv.m = v(i, j - 1, k);
          sv.p = v(i, j + 1, k);
          
          cellboundary(sv);
          
          ////////////////////////////////////////////
          // Right boundary

          // Calculate velocity at right boundary (y+1/2)
          BoutReal vpar = 0.5 * (v(i, j, k) + v(i, j + 1, k));
          BoutReal flux;

          if (mesh->lastY(i) && (j == mesh->yend) && !mesh->periodicY(i)) {
            // Last point in domain

            BoutReal bndryval = 0.5 * (s.c + s.p);
            if (fixflux) {
              // Use mid-point to be consistent with boundary conditions
              flux = bndryval * vpar * vpar;
            } else {
              // Add flux due to difference in boundary values
              flux = s.R * vpar * sv.R + wave_speed(i, j, k) * (s.R * sv.R - bndryval * vpar);
            }
          } else {
            
            // Maximum wave speed in the two cells
            BoutReal amax = BOUTMAX(wave_speed(i, j, k), wave_speed(i, j + 1, k));

            if (vpar > amax) {
              // Supersonic flow out of this cell
              flux = s.R * vpar * sv.R;
            } else if (vpar < -amax) {
              // Supersonic flow into this cell
              flux = 0.0;
            } else {
              // Subsonic flow, so a mix of right and left fluxes
              flux = s.R * 0.5 * (vpar + amax) * sv.R;
            }
          }
          
          result(i, j, k) += flux * flux_factor_rc;
          result(i, j + 1, k) -= flux * flux_factor_rp;

          ////////////////////////////////////////////
          // Calculate at left boundary
          
          vpar = 0.5 * (v(i, j, k) + v(i, j - 1, k));

          if (mesh->firstY(i) && (j == mesh->ystart) && !mesh->periodicY(i)) {
            // First point in domain
            BoutReal bndryval = 0.5 * (s.c + s.m);
            if (fixflux) {
              // Use mid-point to be consistent with boundary conditions
              flux = bndryval * vpar * vpar;
            } else {
              // Add flux due to difference in boundary values
              flux = s.L * vpar * sv.L - wave_speed(i, j, k) * (s.L * sv.L - bndryval * vpar);
            }
          } else {
            
            // Maximum wave speed in the two cells
            BoutReal amax = BOUTMAX(wave_speed(i, j, k), wave_speed(i, j - 1, k));

            if (vpar < -amax) {
              // Supersonic out of this cell
              flux = s.L * vpar * sv.L;
            } else if (vpar > amax) {
              // Supersonic into this cell
              flux = 0.0;
            } else {
              flux = s.L * 0.5 * (vpar - amax) * sv.L;
            }
          }
          
          result(i, j, k) -= flux * flux_factor_lc;
          result(i, j - 1, k) += flux * flux_factor_lm;
          
        }
      }
    }
    return fromFieldAligned(result, "RGN_NOBNDRY");
  }
  
}

BoutReal floor(const BoutReal &var, const BoutReal &f) {
  if (var < f)
    return f;
  return var;
}

/// Returns a copy of input \p var with all values greater than \p f replaced by
/// \p f.
const Field3D ceil(const Field3D &var, BoutReal f, REGION rgn = RGN_ALL) {
  checkData(var);
  Field3D result = copy(var);

  BOUT_FOR(d, var.getRegion(rgn)) {
    if (result[d] > f) {
      result[d] = f;
    }
  }

  return result;
}

// Square function for vectors
Field3D SQ(const Vector3D &v) { return v * v; }

int Hermes::init(bool restarting) {

  auto& opt = Options::root();

  // Switches in model section
  auto& optsc = opt["Hermes"];

  OPTION(optsc, evolve_plasma, true);

  electromagnetic = optsc["electromagnetic"]
                        .doc("Include vector potential psi in Ohm's law?")
                        .withDefault<bool>(true);

  FiniteElMass = optsc["FiniteElMass"]
                     .doc("Include electron inertia in Ohm's law?")
                     .withDefault<bool>(true);

  j_diamag = optsc["j_diamag"]
                 .doc("Diamagnetic current: Vort <-> Pe")
                 .withDefault<bool>(true);
  
  j_par = optsc["j_par"]
              .doc("Parallel current:    Vort <-> Psi")
              .withDefault<bool>(true);
  
  OPTION(optsc, parallel_flow, true);
  OPTION(optsc, parallel_flow_p_term, parallel_flow);
  OPTION(optsc, pe_par, true);
  OPTION(optsc, pe_par_p_term, pe_par);
  OPTION(optsc, resistivity, true);
  OPTION(optsc, thermal_flux, true);

  thermal_force = optsc["thermal_force"]
                    .doc("Force on electrons due to temperature gradients")
                    .withDefault<bool>(true);
  
  OPTION(optsc, electron_viscosity, true);
  OPTION(optsc, ion_viscosity, true);

  electron_neutral = optsc["electron_neutral"]
                       .doc("Include electron-neutral collisions in resistivity?")
                       .withDefault<bool>(true);
  
  ion_neutral = optsc["ion_neutral"]
                  .doc("Include ion-neutral collisions in tau_i?")
                  .withDefault<bool>(false);

  poloidal_flows = optsc["poloidal_flows"]
                       .doc("Include ExB flows in X-Y plane")
                       .withDefault(true);
  
  OPTION(optsc, ion_velocity, true);

  OPTION(optsc, thermal_conduction, true);
  OPTION(optsc, electron_ion_transfer, true);

  OPTION(optsc, neutral_friction, false);
  OPTION(optsc, frecycle, 0.9);

  OPTION(optsc, phi3d, false);

  OPTION(optsc, ne_bndry_flux, true);
  OPTION(optsc, pe_bndry_flux, true);
  OPTION(optsc, vort_bndry_flux, false);

  OPTION(optsc, ramp_mesh, true);
  OPTION(optsc, ramp_timescale, 1e4);

  OPTION(optsc, energy_source, false);

  ion_neutral_rate = optsc["ion_neutral_rate"]
                     .doc("A fixed ion-neutral collision rate, normalised to ion cyclotron frequency.")
                     .withDefault(0.0);

  OPTION(optsc, staggered, false);

  OPTION(optsc, boussinesq, false);

  OPTION(optsc, sinks, false);
  OPTION(optsc, sheath_closure, true);
  OPTION(optsc, drift_wave, false);

  // Cross-field transport
  classical_diffusion = optsc["classical_diffusion"]
                          .doc("Collisional cross-field diffusion, including viscosity")
                          .withDefault<bool>(false);
  OPTION(optsc, anomalous_D, -1);
  OPTION(optsc, anomalous_chi, -1);
  OPTION(optsc, anomalous_nu, -1);
  OPTION(optsc, anomalous_D_nvi, true);
  OPTION(optsc, anomalous_D_pepi, true);

  // Flux limiters
  OPTION(optsc, flux_limit_alpha, -1);
  OPTION(optsc, kappa_limit_alpha, -1);
  OPTION(optsc, eta_limit_alpha, -1);

  // Numerical dissipation terms
  OPTION(optsc, numdiff, -1.0);
  OPTION(optsc, hyper, -1);
  OPTION(optsc, hyperpar, -1);
  OPTION(optsc, low_pass_z, -1);
  OPTION(optsc, x_hyper_viscos, -1.0);
  OPTION(optsc, y_hyper_viscos, -1.0);
  OPTION(optsc, z_hyper_viscos, -1.0);
  OPTION(optsc, scale_num_cs, 1.0);
  OPTION(optsc, floor_num_cs, -1.0);
  OPTION(optsc, vepsi_dissipation, false);
  OPTION(optsc, vort_dissipation, false);

  OPTION(optsc, ne_hyper_z, -1.0);
  OPTION(optsc, pe_hyper_z, -1.0);

  OPTION(optsc, low_n_diffuse, false);
  OPTION(optsc, low_n_diffuse_perp, false);

  OPTION(optsc, resistivity_multiply, 1.0);

  OPTION(optsc, sheath_model, 0);
  OPTION(optsc, sheath_gamma_e, 5.5);
  OPTION(optsc, sheath_gamma_i, 1.0);

  OPTION(optsc, neutral_vwall, 1. / 3);  // 1/3rd Franck-Condon energy at wall
  OPTION(optsc, sheath_yup, true);       // Apply sheath at yup?
  OPTION(optsc, sheath_ydown, true);     // Apply sheath at ydown?
  OPTION(optsc, test_boundaries, false); // Test boundary conditions

  // Fix profiles in SOL
  OPTION(optsc, sol_fix_profiles, false);
  if (sol_fix_profiles) {
    sol_ne = FieldFactory::get()->parse("sol_ne", &optsc);
    sol_te = FieldFactory::get()->parse("sol_te", &optsc);
  }

  OPTION(optsc, radial_buffers, false);
  OPTION(optsc, radial_inner_width, 4);
  OPTION(optsc, radial_outer_width, 4);
  OPTION(optsc, radial_buffer_D, 1.0);

  resistivity_boundary = optsc["resistivity_boundary"]
    .doc("Normalised resistivity in radial boundary region")
    .withDefault(1.0);

  resistivity_boundary_width = optsc["resistivity_boundary_width"]
    .doc("Number of grid cells in radial (x) direction")
    .withDefault(0);
  
  // Output additional information
  OPTION(optsc, verbose, false);    // Save additional fields
  OPTION(optsc, output_ddt, false); // Save time derivatives

  // Normalisation
  OPTION(optsc, Tnorm, 100);  // Reference temperature [eV]
  OPTION(optsc, Nnorm, 1e19); // Reference density [m^-3]
  OPTION(optsc, Bnorm, 1.0);  // Reference magnetic field [T]

  OPTION(optsc, AA, 2.0); // Ion mass (2 = Deuterium)

  output.write("Normalisation Te=%e, Ne=%e, B=%e\n", Tnorm, Nnorm, Bnorm);
  SAVE_ONCE(Tnorm, Nnorm, Bnorm, AA); // Save

  Cs0 = sqrt(qe * Tnorm / (AA * Mp)); // Reference sound speed [m/s]
  Omega_ci = qe * Bnorm / (AA * Mp);  // Ion cyclotron frequency [1/s]
  rho_s0 = Cs0 / Omega_ci;

  mi_me = AA * Mp / Me;
  beta_e = qe * Tnorm * Nnorm / (SQ(Bnorm) / (2. * mu0));

  output.write("\tmi_me=%e, beta_e=%e\n", mi_me, beta_e);
  SAVE_ONCE(mi_me, beta_e);

  output.write("\t Cs=%e, rho_s=%e, Omega_ci=%e\n", Cs0, rho_s0, Omega_ci);
  SAVE_ONCE(Cs0, rho_s0, Omega_ci);

  // Collision times
  BoutReal lambda_ei = 24. - log(sqrt(Nnorm / 1e6) / Tnorm);
  BoutReal lambda_ii = 23. - log(sqrt(2. * Nnorm / 1e6) / pow(Tnorm, 1.5));

  tau_e0 = 1. / (2.91e-6 * (Nnorm / 1e6) * lambda_ei * pow(Tnorm, -3. / 2));
  tau_i0 =
      sqrt(AA) / (4.78e-8 * (Nnorm / 1e6) * lambda_ii * pow(Tnorm, -3. / 2));

  output.write("\ttau_e0=%e, tau_i0=%e\n", tau_e0, tau_i0);

  if (anomalous_D > 0.0) {
    // Normalise
    anomalous_D /= rho_s0 * rho_s0 * Omega_ci; // m^2/s
    output.write("\tnormalised anomalous D_perp = %e\n", anomalous_D);
  }
  if (anomalous_chi > 0.0) {
    // Normalise
    anomalous_chi /= rho_s0 * rho_s0 * Omega_ci; // m^2/s
    output.write("\tnormalised anomalous chi_perp = %e\n", anomalous_chi);
  }
  if (anomalous_nu > 0.0) {
    // Normalise
    anomalous_nu /= rho_s0 * rho_s0 * Omega_ci; // m^2/s
    output.write("\tnormalised anomalous nu_perp = %e\n", anomalous_nu);
  }

  if (ramp_mesh) {
    Jpar0 = 0.0;
  } else {
    // Read equilibrium current density
    // GRID_LOAD(Jpar0);
    // Jpar0 /= qe*Nnorm*Cs0;
    Jpar0 = 0.0;
  }

  FieldFactory fact(mesh);

  if (sinks) {
    std::string source = optsc["sink_invlpar"].withDefault<std::string>("0.05"); // 20 m
    sink_invlpar = fact.create2D(source);
    sink_invlpar *= rho_s0; // Normalise
    SAVE_ONCE(sink_invlpar);

    if (drift_wave) {
      alpha_dw = fact.create2D("Hermes:alpha_dw");
      SAVE_ONCE(alpha_dw);
    }
  }

  // Get switches from each variable section
  auto& optne = opt["Ne"];
  NeSource = optne["source"].doc("Source term in ddt(Ne)").withDefault(Field3D{0.0});
  NeSource /= Omega_ci;
  Sn = DC(NeSource);

  // Inflowing density carries momentum
  OPTION(optne, density_inflow, false);

  auto& optpe = opt["Pe"];
  PeSource = optpe["source"].withDefault(Field3D{0.0});
  PeSource /= Omega_ci;
  Spe = DC(PeSource);

  auto& optpi = opt["Pi"];
  PiSource = optpi["source"].withDefault(Field3D{0.0});
  PiSource /= Omega_ci;
  Spi = DC(PiSource);

  OPTION(optsc, core_sources, false);
  if (core_sources) {
    for (int x = mesh->xstart; x <= mesh->xend; x++) {
      if (!mesh->periodicY(x)) {
        // Not periodic, so not in core
        for (int y = mesh->ystart; y <= mesh->yend; y++) {
          Sn(x, y) = 0.0;
          Spe(x, y) = 0.0;
        }
      }
    }
  }

  // Mid-plane power flux q_||
  // Midplane power specified in Watts per m^2
  Field2D qfact;
  GRID_LOAD(qfact); // Factor to multiply to get volume source
  Field2D qmid = optpe["midplane_power"].withDefault(Field2D{0.0}) * qfact;
  // Normalise from W/m^3
  qmid /= qe * Tnorm * Nnorm * Omega_ci;
  Spe += (2. / 3) * qmid;

  // Add variables to solver
  SOLVE_FOR(Ne, Pe, Pi);
  EvolvingVars.add(Ne, Pe, Pi);

  if (output_ddt) {
    SAVE_REPEAT(ddt(Ne), ddt(Pe), ddt(Pi));
  }

  if (j_par || j_diamag) {
    // Have a source of vorticity
    solver->add(Vort, "Vort");
    EvolvingVars.add(Vort);
    if (output_ddt) {
      SAVE_REPEAT(ddt(Vort));
    }
  } else {
    Vort = 0.0;
  }

  if (electromagnetic || FiniteElMass) {
    solver->add(VePsi, "VePsi");
    EvolvingVars.add(VePsi);
    if (output_ddt) {
      SAVE_REPEAT(ddt(VePsi));
    }
  } else {
    // If both electrostatic and zero electron mass,
    // then Ohm's law has no time-derivative terms,
    // but is calculated from other evolving quantities
    VePsi = 0.0;
  }

  if (ion_velocity) {
    solver->add(NVi, "NVi");
    EvolvingVars.add(NVi);
    if (output_ddt) {
      SAVE_REPEAT(ddt(NVi));
    }
  } else {
    NVi = 0.0;
  }

  if (verbose) {
    SAVE_REPEAT(Ti);
    if (electron_ion_transfer) {
      SAVE_REPEAT(Wi);
    }
    if (ion_velocity) {
      SAVE_REPEAT(Vi);
    }
  }

  OPTION(optsc, adapt_source, false);
  if (adapt_source) {
    // Adaptive sources to match profiles

    // PI controller, including an integrated difference term
    OPTION(optsc, source_p, 1e-2);
    OPTION(optsc, source_i, 1e-6);

    Field2D Snsave = copy(Sn);
    Field2D Spesave = copy(Spe);
    Field2D Spisave = copy(Spi);
    SOLVE_FOR(Sn, Spe, Spi);
    Sn = Snsave;
    Spe = Spesave;
    Spi = Spisave;
  } else {
    SAVE_ONCE(Sn, Spe, Spi);
  }

  /////////////////////////////////////////////////////////
  // Load metric tensor from the mesh, passing length and B
  // field normalisations
  TRACE("Loading metric tensor");

  Coordinates *coord = mesh->getCoordinates();

  if (optsc["loadmetric"]
          .doc("Load Rxy, Bpxy etc. to create orthogonal metric?")
          .withDefault(true)) {
    LoadMetric(rho_s0, Bnorm);
  } else {
    Coordinates *coord = mesh->getCoordinates();
    // To use non-orthogonal metric
    // Normalise
    coord->dx /= rho_s0 * rho_s0 * Bnorm;
    coord->Bxy /= Bnorm;
    // Metric is in grid file - just need to normalise
    coord->g11 /= (Bnorm * Bnorm * rho_s0 * rho_s0);
    coord->g22 *= (rho_s0 * rho_s0);
    coord->g33 *= (rho_s0 * rho_s0);
    coord->g12 /= Bnorm;
    coord->g13 /= Bnorm;
    coord->g23 *= (rho_s0 * rho_s0);

    coord->J *= Bnorm / rho_s0;

    coord->g_11 *= (Bnorm * Bnorm * rho_s0 * rho_s0);
    coord->g_22 /= (rho_s0 * rho_s0);
    coord->g_33 /= (rho_s0 * rho_s0);
    coord->g_12 *= Bnorm;
    coord->g_13 *= Bnorm;
    coord->g_23 /= (rho_s0 * rho_s0);

    coord->geometry(); // Calculate other metrics
  }

  sqrtB = sqrt(coord->Bxy); // B^(1/2)
  B32 = sqrtB * coord->Bxy; // B^(3/2)

  SAVE_ONCE(B32);
  
  /////////////////////////////////////////////////////////
  // Neutral models

  TRACE("Initialising neutral models");
  neutrals = NeutralModel::create(solver, mesh, Options::root()["neutral"]);

  // Set normalisations
  if (neutrals) {
    neutrals->setNormalisation(Tnorm, Nnorm, Bnorm, rho_s0, Omega_ci);
  }

  /////////////////////////////////////////////////////////
  // Impurities
  TRACE("Impurities");
  
  impurity_adas = optsc["impurity_adas"]
                      .doc("Use Atomic++ interface to ADAS")
                      .withDefault<bool>(false);

  if (impurity_adas) {

    fimp = optsc["impurity_fraction"]
               .doc("Fixed fraction ADAS impurity, multiple of electron density")
               .withDefault(0.0);

    string impurity_species =
        optsc["impurity_species"]
            .doc("Short name of the ADAS species e.g. 'c' or 'ne'")
            .withDefault("c");

    impurity = new ImpuritySpecies(impurity_species);
  }

  carbon_fraction = optsc["carbon_fraction"]
          .doc("Include a fixed fraction carbon impurity. < 0 means none.")
          .withDefault(-1.);
  if (carbon_fraction > 0.0) {
    SAVE_REPEAT(Rzrad);
    SAVE_ONCE(carbon_fraction);
    carbon_rad = new HutchinsonCarbonRadiation();
  }

  if ((carbon_fraction > 0.0) || impurity_adas) {
    // Save impurity radiation
    SAVE_REPEAT(Rzrad);
  }
  
  /////////////////////////////////////////////////////////
  // Read profiles from the mesh
  TRACE("Reading profiles");

  Field2D NeMesh, TeMesh, TiMesh;
  if (mesh->get(NeMesh, "Ne0")) {
    // No Ne0. Try Ni0
    if (mesh->get(NeMesh, "Ni0")) {
      output << "WARNING: Neither Ne0 nor Ni0 found in mesh input\n";
    }
  }
  NeMesh *= 1e20; // Convert to m^-3

  NeMesh /= Nnorm; // Normalise

  if (mesh->get(TeMesh, "Te0")) {
    // No Te0
    output << "WARNING: Te0 not found in mesh\n";
    // Try to read Ti0
    if (mesh->get(TeMesh, "Ti0")) {
      // No Ti0 either
      output << "WARNING: No Te0 or Ti0. Setting TeMesh to 0.0\n";
      TeMesh = 0.0;
    }
  }

  TeMesh /= Tnorm; // Normalise

  if (mesh->get(TiMesh, "Ti0")) {
    // No Ti0
    output << "WARNING: Ti0 not found in mesh. Setting to TeMesh\n";
    TiMesh = TeMesh;
  }
  TiMesh /= Tnorm; // Normalise
  PiTarget = NeMesh * TiMesh;

  NeTarget = NeMesh;
  PeTarget = NeMesh * TeMesh;

  if (!restarting && !ramp_mesh) {
    if (optsc["startprofiles"].withDefault(true)) {
      Ne += NeMesh; // Add profiles in the mesh file

      Pe += NeMesh * TeMesh;

      Pi += NeMesh * TiMesh;
      // Check for negatives
      if (min(Pi, true) < 0.0) {
        throw BoutException("Starting ion pressure is negative");
      }
      if (max(Pi, true) < 1e-5) {
        throw BoutException("Starting ion pressure is too small");
      }
      mesh->communicateXZ(Pi);
    }

    // Check for negatives
    if (min(Ne, true) < 0.0) {
      throw BoutException("Starting density is negative");
    }
    if (max(Ne, true) < 1e-5) {
      throw BoutException("Starting density is too small");
    }

    if (min(Pe, true) < 0.0) {
      throw BoutException("Starting pressure is negative");
    }
    if (max(Pe, true) < 1e-5) {
      throw BoutException("Starting pressure is too small");
    }

    mesh->communicateXZ(Ne, Pe);
  }

  /////////////////////////////////////////////////////////
  // Sources (after metric)

  // Multiply sources by g11
  OPTION(optsc, source_vary_g11, false);
  if (source_vary_g11) {
    // Average metric tensor component
    g11norm = coord->g11 / averageY(coord->g11);

    NeSource *= g11norm;
    PeSource *= g11norm;
    PiSource *= g11norm;
  }

  /////////////////////////////////////////////////////////
  // Read curvature components
  TRACE("Reading curvature");

  try {
    Curlb_B.covariant = false; // Contravariant
    mesh->get(Curlb_B, "bxcv");
  
  } catch (BoutException &e) {
    try {
      // May be 2D, reading as 3D
      Vector2D curv2d;
      curv2d.covariant = false;
      mesh->get(curv2d, "bxcv");
      Curlb_B = curv2d;
    } catch (BoutException &e) {
      if (j_diamag) {
        // Need curvature
        throw;
      } else {
        output_warn.write("No curvature vector in input grid");
        Curlb_B = 0.0;
      }
    }
  }
  
  if (Options::root()["mesh"]["paralleltransform"].as<std::string>() == "shifted") {
    Field2D I;
    mesh->get(I, "sinty");
    Curlb_B.z += I * Curlb_B.x;
  }

  Curlb_B.x /= Bnorm;
  Curlb_B.y *= rho_s0 * rho_s0;
  Curlb_B.z *= rho_s0 * rho_s0;

  Curlb_B *= 2. / coord->Bxy;

  SAVE_REPEAT(phi);

  if (j_par) {
    SAVE_REPEAT(Ve);

    if (electromagnetic)
      SAVE_REPEAT(psi);
  }

  OPTION(optsc, split_n0, false); // Split into n=0 and n~=0
  OPTION(optsc, split_n0_psi, split_n0);
  // Phi solver
  if (phi3d) {
#ifdef PHISOLVER
    phiSolver3D = Laplace3D::create();
#endif
  } else {
    if (split_n0) {
      // Create an XY solver for n=0 component
      laplacexy = new LaplaceXY(mesh);
      // Set coefficients for Boussinesq solve
      laplacexy->setCoefs(1. / SQ(coord->Bxy), 0.0);
      phi2D = 0.0; // Starting guess
    }

    // Create an XZ solver
    OPTION(optsc, newXZsolver, false);
    if (newXZsolver) {
      // Test new LaplaceXZ solver
      newSolver = LaplaceXZ::create(mesh);
      // Set coefficients for Boussinesq solve
      newSolver->setCoefs(1. / SQ(coord->Bxy), 0.0);
    } else {
      // Use older Laplacian solver
      phiSolver = Laplacian::create(&opt["phiSolver"]);
      // Set coefficients for Boussinesq solve
      phiSolver->setCoefC(1. / SQ(coord->Bxy));
    }
    phi = 0.0;
  }

  // Apar (Psi) solver
  aparSolver = LaplaceXZ::create(mesh, &opt["aparSolver"], CELL_CENTRE);
  if (split_n0_psi) {
    // Use another XY solver for n=0 psi component
    aparXY = new LaplaceXY(mesh);
    psi2D = 0.0;
  }

  Ve.setBoundary("Ve");
  nu.setBoundary("nu");
  phi.setBoundary("phi"); // For y boundaries
  Jpar.setBoundary("Jpar");

  nu = 0.0;
  kappa_epar = 0.0;
  Dn = 0.0;

  SAVE_REPEAT(Telim, Tilim);

  if (verbose) {
    // Save additional fields
    SAVE_REPEAT(Jpar); // Parallel current

    SAVE_REPEAT(tau_e, tau_i);

    SAVE_REPEAT(kappa_epar); // Parallel electron heat conductivity
    SAVE_REPEAT(kappa_ipar); // Parallel ion heat conductivity

    if (resistivity) {
      SAVE_REPEAT(nu); // Parallel resistivity
    }

    // SAVE_REPEAT2(wall_flux, wall_power);

    if (ion_viscosity) {
      // Ion parallel stress tensor
      SAVE_REPEAT(Pi_ci, Pi_ciperp, Pi_cipar);
    }

    // Sources added to Ne, Pe and Pi equations
    SAVE_REPEAT(NeSource, PeSource, PiSource);
  }

  psi = phi = 0.0;

  // Preconditioner
  setPrecon((preconfunc)&Hermes::precon);
  
  return 0;
}

int Hermes::rhs(BoutReal t) {
  Coordinates *coord = mesh->getCoordinates();
  
  if (!evolve_plasma) {
    Ne = 0.0;
    Pe = 0.0;
    Pi = 0.0;
    Vort = 0.0;
    VePsi = 0.0;
    NVi = 0.0;
    sheath_model = 0;
  }

  // Communicate evolving variables
  // Note: Parallel slices are not calculated because parallel derivatives
  // are calculated using field aligned quantities
  
  mesh->communicateXZ(EvolvingVars);
  for (auto* f : EvolvingVars.field3d()) {
    f->clearParallelSlices(); // Make sure no parallel slices
  }

  Field3D Nelim = floor(Ne, 1e-5);

  Te = Pe / Nelim;
  Vi = NVi / Nelim;

  Telim = floor(Te, 0.1 / Tnorm);

  Field3D Pelim = Telim * Nelim;

  Field3D logPelim = log(Pelim);
  logPelim.applyBoundary("neumann");

  Ti = Pi / Nelim;
  Tilim = floor(Ti, 0.1 / Tnorm);
  Field3D Pilim = Tilim * Nelim;

  // Set radial boundary conditions on Te, Ti, Vi
  //
  if (mesh->firstX()) {
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      for (int k = 0; k < mesh->LocalNz; k++) {
        BoutReal ne_bndry = 0.5 * (Ne(1, j, k) + Ne(2, j, k));
        if (ne_bndry < 1e-5)
          ne_bndry = 1e-5;
        BoutReal pe_bndry = 0.5 * (Pe(1, j, k) + Pe(2, j, k));
        BoutReal pi_bndry = 0.5 * (Pi(1, j, k) + Pi(2, j, k));

        BoutReal te_bndry = pe_bndry / ne_bndry;
        BoutReal ti_bndry = pi_bndry / ne_bndry;

        Te(1, j, k) = 2. * te_bndry - Te(2, j, k);
        Ti(1, j, k) = 2. * ti_bndry - Ti(2, j, k);
        Vi(0, j, k) = Vi(1, j, k) = Vi(2, j, k);

        if (te_bndry < 0.1 / Tnorm)
          te_bndry = 0.1 / Tnorm;
        if (ti_bndry < 0.1 / Tnorm)
          ti_bndry = 0.1 / Tnorm;

        Telim(1, j, k) = 2. * te_bndry - Telim(2, j, k);
        Tilim(1, j, k) = 2. * ti_bndry - Tilim(2, j, k);
      }
    }
  }
  if (mesh->lastX()) {
    int n = mesh->LocalNx;
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      for (int k = 0; k < mesh->LocalNz; k++) {
        BoutReal ne_bndry = 0.5 * (Ne(n - 1, j, k) + Ne(n - 2, j, k));
        if (ne_bndry < 1e-5)
          ne_bndry = 1e-5;
        BoutReal pe_bndry = 0.5 * (Pe(n - 1, j, k) + Pe(n - 2, j, k));
        BoutReal pi_bndry = 0.5 * (Pi(n - 1, j, k) + Pi(n - 2, j, k));

        BoutReal te_bndry = pe_bndry / ne_bndry;
        BoutReal ti_bndry = pi_bndry / ne_bndry;

        Te(n - 1, j, k) = 2. * te_bndry - Te(n - 2, j, k);
        Ti(n - 1, j, k) = 2. * ti_bndry - Ti(n - 2, j, k);
        Vi(n - 1, j, k) = Vi(n - 2, j, k);

        if (te_bndry < 0.1 / Tnorm)
          te_bndry = 0.1 / Tnorm;
        if (ti_bndry < 0.1 / Tnorm)
          ti_bndry = 0.1 / Tnorm;

        Telim(n - 1, j, k) = 2. * te_bndry - Telim(n - 2, j, k);
        Tilim(n - 1, j, k) = 2. * ti_bndry - Tilim(n - 2, j, k);
      }
    }
  }

  // Are there any currents? If not, then there is no source
  // for vorticity, phi = 0 and jpar = 0
  bool currents = j_par | j_diamag;

  // Local sound speed. Used for parallel advection operator
  // Assumes isothermal electrons, adiabatic ions
  // The factor scale_num_cs can be used to test sensitity
  Field3D sound_speed = scale_num_cs * sqrt(Telim + Tilim * (5. / 3));
  if (floor_num_cs > 0.0) {
    // Apply a floor function to the sound speed
    sound_speed = floor(sound_speed, floor_num_cs);
  }
  sound_speed.applyBoundary("neumann");
  
  //////////////////////////////////////////////////////////////
  // Calculate electrostatic potential phi
  //
  //

  TRACE("Electrostatic potential");

  if (!currents) {
    // Disabling electric fields
    // phi = 0.0; // Already set in initialisation
  } else {
    // Solve phi from Vorticity
    if (phi3d) {
#ifdef PHISOLVER
      phiSolver3D->setCoefC(Ne / SQ(coord->Bxy));
      // phi.setBoundaryTo(3.*Te);
      if (mesh->lastX()) {
        for (int i = mesh->xend + 1; i < mesh->LocalNx; i++)
          for (int j = mesh->ystart; j <= mesh->yend; j++)
            for (int k = 0; k < mesh->LocalNz; k++) {
              phi(i, j, k) = 3. * Te(i, j, k);
            }
      }
      phi = phiSolver3D->solve(Vort, phi);
#endif
    } else {
      // Phi flags should be set in BOUT.inp
      // phiSolver->setInnerBoundaryFlags(INVERT_DC_GRAD);
      // phiSolver->setOuterBoundaryFlags(INVERT_SET);

      // Sheath multiplier Te -> phi (2.84522 for Deuterium)
      BoutReal sheathmult = log(0.5 * sqrt(mi_me / PI));

      if (boussinesq) {

        if (split_n0) {
          ////////////////////////////////////////////
          // Boussinesq, split
          // Split into axisymmetric and non-axisymmetric components
          Field2D Vort2D = DC(Vort); // n=0 component

          // Set the boundary to 2.8*Te
          // phi2D.setBoundaryTo(sheathmult * floor(Te.DC(), 0.0));
          phi2D.setBoundaryTo(sheathmult * DC(Telim));

          phi2D = laplacexy->solve(Vort2D, phi2D);

          // Solve non-axisymmetric part using X-Z solver
          if (newXZsolver) {
            newSolver->setCoefs(1. / SQ(coord->Bxy), 0.0);
            phi = newSolver->solve(Vort - Vort2D, phi);
          } else {
            phiSolver->setCoefC(1. / SQ(coord->Bxy));
            // phi = phiSolver->solve((Vort-Vort2D)*SQ(coord->Bxy), phi);
            phi = phiSolver->solve((Vort - Vort2D) * SQ(coord->Bxy),
                                   sheathmult * (Telim - DC(Telim)));
          }
          phi += phi2D; // Add axisymmetric part
        } else {
          ////////////////////////////////////////////
          // Boussinesq, non-split
          // Solve all components using X-Z solver

          if (newXZsolver) {
            // Use the new LaplaceXY solver
            // newSolver->setCoefs(1./SQ(coord->Bxy), 0.0); // Set when
            // initialised
            phi = newSolver->solve(Vort, phi);
          } else {
            // Use older Laplacian solver
            // phiSolver->setCoefC(1./SQ(coord->Bxy)); // Set when initialised
            phi = phiSolver->solve(Vort * SQ(coord->Bxy), (sheathmult + log(sqrt(Telim / (Telim + Tilim)))) * Telim );
          }
        }
        
        // Hot ion term in vorticity
        phi -= Pi;

      } else {
        ////////////////////////////////////////////
        // Non-Boussinesq
        //
        throw BoutException("Non-Boussinesq not implemented yet");

        phiSolver->setCoefC(Nelim / SQ(coord->Bxy));
        phi = phiSolver->solve(Vort * SQ(coord->Bxy) / Nelim, sheathmult * Telim);
      }
    }
    phi.applyBoundary(t);
    mesh->communicate(phi);
    phi.clearParallelSlices();
  }

  //////////////////////////////////////////////////////////////
  // Calculate perturbed magnetic field psi
  TRACE("Calculating psi");

  if (!currents) {
    // No magnetic fields or currents
    psi = 0.0;
    Jpar = 0.0;
    Ve = 0.0; // Ve will be set after the sheath boundaries below
  } else {
    // Calculate electomagnetic potential psi from VePsi
    // VePsi = Ve - Vi + 0.5 * mi_me * beta_e * psi
    // where the first term comes from finite electron mass, and the second
    // from the parallel component of the electric field
    // Note that psi is -A_|| so Jpar = Delp2(psi)
    if (electromagnetic) {
      if (FiniteElMass) {
        // Solve Helmholtz equation for psi

        // aparSolver->setCoefA(-Ne*0.5*mi_me*beta_e);
        Field2D NDC = DC(Ne);
        aparSolver->setCoefs(1.0, -NDC * 0.5 * mi_me * beta_e);
        // aparSolver->setCoefs(1.0, -Ne*0.5*mi_me*beta_e);
        if (split_n0_psi) {
          // Solve n=0 component separately

          aparXY->setCoefs(1.0, -NDC * 0.5 * mi_me * beta_e);

          Field2D JDC = -NDC * DC(VePsi);
          aparXY->solve(JDC, psi2D);

          psi = aparSolver->solve(-NDC * VePsi - JDC, psi - psi2D) + psi2D;
        } else {
          psi = aparSolver->solve(-NDC * VePsi, psi);
          // psi = aparSolver->solve(-Ne*VePsi, psi);
        }

        Ve = VePsi - 0.5 * mi_me * beta_e * psi + Vi;

        Ve.applyBoundary(t);
        mesh->communicate(Ve, psi);

        Jpar = Ne * (Vi - Ve);
        Jpar.applyBoundary(t);

      } else {
        // Zero electron mass
        // No Ve term in VePsi, only electromagnetic term

        psi = VePsi / (0.5 * mi_me * beta_e);

        // Ve = (NVi - Delp2(psi)) / Nelim;
        Jpar = FV::Div_a_Laplace_perp(1.0, psi);
        mesh->communicate(Jpar);

        Jpar.applyBoundary(t);
        Ve = (NVi - Jpar) / Nelim;
      }

      // psi -= psi.DC(); // Remove toroidal average, only keep fluctuations
    } else {
      // Electrostatic
      psi = 0.0;
      if (FiniteElMass) {
        // No psi contribution to VePsi
        Ve = VePsi + Vi;
      } else {
        // Zero electron mass and electrostatic.
        // Special case where Ohm's law has no time-derivatives
        mesh->communicate(phi);

        Ve = Vi + (Grad_parP(phi) - Grad_parP(Pe) / Ne) / nu;

        if (thermal_force) {
          Ve -= 0.71 * Grad_parP(Te) / nu;
        }
      }

      Ve.applyBoundary(t);
      // Communicate auxilliary variables
      mesh->communicate(Ve);

      Jpar = NVi - Ne * Ve;
    }
    // Ve -= Jpar0 / Ne; // Equilibrium current density
  }

  //////////////////////////////////////////////////////////////
  // Sheath boundary conditions on Y up and Y down
  //
  // NOTE: Have to apply parallel boundary conditions in field aligned coordinates
  // so shift to and then from field aligned

  TRACE("Sheath boundaries");

  // Manually setting boundary conditions
  Ne.clearParallelSlices();
  Pe.clearParallelSlices();
  Pi.clearParallelSlices();
  Pelim.clearParallelSlices();
  Pilim.clearParallelSlices();
  Te.clearParallelSlices();
  Ti.clearParallelSlices();
  Ve.clearParallelSlices();
  Vi.clearParallelSlices();
  Vi.clearParallelSlices();
  NVi.clearParallelSlices();
  phi.clearParallelSlices();
  Vort.clearParallelSlices();
  Jpar.clearParallelSlices();
  
  // Shift to field aligned

  Ne = toFieldAligned(Ne);
  
  Te = toFieldAligned(Te);
  Ti = toFieldAligned(Ti);
  Telim = toFieldAligned(Telim);
  Tilim = toFieldAligned(Tilim);
  
  Pe = toFieldAligned(Pe);
  Pi = toFieldAligned(Pi);
  Pelim = toFieldAligned(Pelim);
  Pilim = toFieldAligned(Pilim);
  
  Ve = toFieldAligned(Ve);
  Vi = toFieldAligned(Vi);
  NVi = toFieldAligned(NVi);

  phi = toFieldAligned(phi);
  Jpar = toFieldAligned(Jpar);
  Vort = toFieldAligned(Vort);
  
  if (sheath_ydown) {
    switch (sheath_model) {
    case 0: { // Normal Bohm sheath
      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          // Zero-gradient density
          BoutReal nesheath = floor(Ne(r.ind, mesh->ystart, jz), 0.0);

          // Temperature at the sheath entrance
          BoutReal tesheath = floor(Te(r.ind, mesh->ystart, jz), 0.0);
          BoutReal tisheath = floor(Ti(r.ind, mesh->ystart, jz), 0.0);

          // Zero-gradient potential
          BoutReal phisheath = phi(r.ind, mesh->ystart, jz);

          // Ion velocity goes to the sound speed. Note negative since out of the domain
          BoutReal visheath = -sqrt(tesheath + tisheath);

          if (Vi(r.ind, mesh->ystart, jz) < visheath) {
            // If plasma is faster, go to plasma velocity
            visheath = Vi(r.ind, mesh->ystart, jz);
          }

          // Sheath current
          // Note that phi/Te >= 0.0 since for phi < 0
          // vesheath is the electron saturation current
          BoutReal phi_te =
              floor(phisheath / Telim(r.ind, mesh->ystart, jz), 0.0);
          BoutReal vesheath =
              -sqrt(tesheath) * (sqrt(mi_me) / (2. * sqrt(PI))) * exp(-phi_te);

          // J = n*(Vi - Ve)
          BoutReal jsheath = nesheath * (visheath - vesheath);
          if (nesheath < 1e-10) {
            vesheath = visheath;
            jsheath = 0.0;
          }

          // Apply boundary condition half-way between cells
          for (int jy = mesh->ystart - 1; jy >= 0; jy--) {
            // Neumann conditions
            Ne(r.ind, jy, jz) = nesheath;
            phi(r.ind, jy, jz) = phisheath;
            Vort(r.ind, jy, jz) = Vort(r.ind, mesh->ystart, jz);

            // Here zero-gradient Te, heat flux applied later
            Te(r.ind, jy, jz) = Te(r.ind, mesh->ystart, jz);
            Ti(r.ind, jy, jz) = Ti(r.ind, mesh->ystart, jz);

            Pe(r.ind, jy, jz) = Pe(r.ind, mesh->ystart, jz);
            Pelim(r.ind, jy, jz) = Pelim(r.ind, mesh->ystart, jz);
            Pi(r.ind, jy, jz) = Pi(r.ind, mesh->ystart, jz);
            Pilim(r.ind, jy, jz) = Pilim(r.ind, mesh->ystart, jz);

            // Dirichlet conditions
            Vi(r.ind, jy, jz) = 2. * visheath - Vi(r.ind, mesh->ystart, jz);
            Ve(r.ind, jy, jz) = 2. * vesheath - Ve(r.ind, mesh->ystart, jz);
            Jpar(r.ind, jy, jz) = 2. * jsheath - Jpar(r.ind, mesh->ystart, jz);
            NVi(r.ind, jy, jz) =
                2. * nesheath * visheath - NVi(r.ind, mesh->ystart, jz);
          }
        }
      }
      break;
    }
    case 1: {
      /*
        Loizu boundary conditions

        Temperature

        Grad_par(Te) = 0

        Density equation

        Grad_par(n) = -(n/Cs) Grad_par(Vi)
        -> n_p - n_m = - (n_p + n_m)/(2Cs) (Vi_p - Vi_m)

        Pressure

        Grad_par(Pe) = Te Grad_par(n)

        Potential
        Grad_par(phi) = -Cs Grad_par(Vi)
        -> phi_p - phi_m = -Cs (Vi_p - Vi_m)
       */

      throw BoutException("Sorry, no Loizu boundary for hot ions yet");

      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          // Temperature at the sheath entrance
          BoutReal tesheath = floor(Te(r.ind, mesh->ystart, jz), 0.0);

          // Zero gradient Te
          Te(r.ind, mesh->ystart - 1, jz) = Te(r.ind, mesh->ystart, jz);
          BoutReal Cs = sqrt(tesheath); // Sound speed

          // Ion velocity goes to the sound speed
          // Dirichlet boundary condition
          Vi(r.ind, mesh->ystart - 1, jz) =
              -2. * Cs - Vi(r.ind, mesh->ystart, jz);

          BoutReal g = 0.0;
          if (tesheath > 0.1 / Tnorm) {
            // Only divide by Cs if the temperature is greater than 0.1eV
            // to avoid divide-by-zero errors
            g = (Vi(r.ind, mesh->ystart - 1, jz) -
                 Vi(r.ind, mesh->ystart, jz)) /
                (-2. * Cs);
          }

          // Mixed boundary condition on n
          Ne(r.ind, mesh->ystart - 1, jz) =
              Ne(r.ind, mesh->ystart, jz) * (1 - g) / (1 + g);

          // Make sure nesheath doesn't go negative
          Ne(r.ind, mesh->ystart - 1, jz) = floor(
              Ne(r.ind, mesh->ystart - 1, jz), -Ne(r.ind, mesh->ystart, jz));

          // Density at the sheath
          BoutReal nesheath = 0.5 * (Ne(r.ind, mesh->ystart, jz) +
                                     Ne(r.ind, mesh->ystart - 1, jz));

          // Momentum
          NVi(r.ind, mesh->ystart - 1, jz) = NVi(r.ind, mesh->ystart, jz);
          if (NVi(r.ind, mesh->ystart - 1, jz) > 0.0) {
            // Limit flux to be <= 0
            NVi(r.ind, mesh->ystart - 1, jz) = -NVi(r.ind, mesh->ystart, jz);
          }

          // Pressure
          Pe(r.ind, mesh->ystart - 1, jz) =
              Pe(r.ind, mesh->ystart, jz) +
              tesheath * (Ne(r.ind, mesh->ystart - 1, jz) -
                          Ne(r.ind, mesh->ystart, jz));

          // Potential
          phi(r.ind, mesh->ystart - 1, jz) =
              phi(r.ind, mesh->ystart, jz) -
              Cs * (Vi(r.ind, mesh->ystart - 1, jz) -
                    Vi(r.ind, mesh->ystart, jz));
          BoutReal phisheath = 0.5 * (phi(r.ind, mesh->ystart, jz) +
                                      phi(r.ind, mesh->ystart - 1, jz));

          // Sheath current
          BoutReal phi_te = floor(phisheath / tesheath, 0.0);
          BoutReal vesheath =
              -Cs * (sqrt(mi_me) / (2. * sqrt(PI))) * exp(-phi_te);

          BoutReal jsheath = nesheath * (-Cs - vesheath);

          // Electron velocity
          Ve(r.ind, mesh->ystart - 1, jz) =
              2. * vesheath - Ve(r.ind, mesh->ystart, jz);

          // Parallel velocity
          Jpar(r.ind, mesh->ystart - 1, jz) =
              2. * jsheath - Jpar(r.ind, mesh->ystart, jz);

          if (currents && !finite(Jpar(r.ind, mesh->ystart - 1, jz))) {
            output.write("JPAR: %d, %d: %e, %e, %e, %e\n", r.ind, jz, jsheath,
                         vesheath, Cs, nesheath);
            output.write(" -> %e, %e, %e\n", Ne(r.ind, mesh->ystart, jz),
                         Ne(r.ind, mesh->ystart - 1, jz), g);
            exit(1);
          }

          // Constant gradient on other cells
          for (int jy = mesh->ystart - 2; jy >= 0; jy--) {
            Vi(r.ind, jy, jz) =
                2. * Vi(r.ind, jy + 1, jz) - Vi(r.ind, jy + 2, jz);
            Ve(r.ind, jy, jz) =
                2. * Ve(r.ind, jy + 1, jz) - Ve(r.ind, jy + 2, jz);
            NVi(r.ind, jy, jz) =
                2. * NVi(r.ind, jy + 1, jz) - NVi(r.ind, jy + 2, jz);

            Ne(r.ind, jy, jz) =
                2. * Ne(r.ind, jy + 1, jz) - Ne(r.ind, jy + 2, jz);
            Te(r.ind, jy, jz) =
                2. * Te(r.ind, jy + 1, jz) - Te(r.ind, jy + 2, jz);
            Pe(r.ind, jy, jz) =
                2. * Pe(r.ind, jy + 1, jz) - Pe(r.ind, jy + 2, jz);

            phi(r.ind, jy, jz) =
                2. * phi(r.ind, jy + 1, jz) - phi(r.ind, jy + 2, jz);
            Vort(r.ind, jy, jz) =
                2. * Vort(r.ind, jy + 1, jz) - Vort(r.ind, jy + 2, jz);
            Jpar(r.ind, jy, jz) =
                2. * Jpar(r.ind, jy + 1, jz) - Jpar(r.ind, jy + 2, jz);
          }
        }
      }
      break;
    }
    case 2: { // Bohm sheath with free density
      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          // Zero-gradient density
          BoutReal nesheath = 0.5 * (3. * Ne(r.ind, mesh->ystart, jz) -
                                     Ne(r.ind, mesh->ystart + 1, jz));
          if (nesheath < 0.0)
            nesheath = 0.0;

          // Temperature at the sheath entrance
          BoutReal tesheath = floor(Te(r.ind, mesh->ystart, jz), 0.0);
          BoutReal tisheath = floor(Ti(r.ind, mesh->ystart, jz), 0.0);

          // Zero-gradient potential
          BoutReal phisheath = phi(r.ind, mesh->ystart, jz);

          // Ion velocity goes to the sound speed
          BoutReal visheath =
              -sqrt(tesheath + tisheath); // Sound speed outwards

          if (Vi(r.ind, mesh->ystart, jz) < visheath) {
            // If plasma is faster, go to plasma velocity
            visheath = Vi(r.ind, mesh->ystart, jz);
          }

          // Sheath current
          BoutReal phi_te =
              floor(phisheath / Telim(r.ind, mesh->ystart, jz), 0.0);
          BoutReal vesheath =
              -sqrt(tesheath) * (sqrt(mi_me) / (2. * sqrt(PI))) * exp(-phi_te);
          // J = n*(Vi - Ve)
          BoutReal jsheath = nesheath * (visheath - vesheath);
          if (nesheath < 1e-10) {
            vesheath = visheath;
            jsheath = 0.0;
          }

          // Apply boundary condition half-way between cells
          for (int jy = mesh->ystart - 1; jy >= 0; jy--) {
            // Neumann conditions
            phi(r.ind, jy, jz) = phisheath;
            Vort(r.ind, jy, jz) = Vort(r.ind, mesh->ystart, jz);

            // Here zero-gradient Te, heat flux applied later
            Te(r.ind, jy, jz) = Te(r.ind, mesh->ystart, jz);
            Ti(r.ind, jy, jz) = Ti(r.ind, mesh->ystart, jz);

            // Dirichlet conditions
            Ne(r.ind, jy, jz) = 2. * nesheath - Ne(r.ind, mesh->ystart, jz);
            Pe(r.ind, jy, jz) =
                2. * nesheath * tesheath - Pe(r.ind, mesh->ystart, jz);
            Pi(r.ind, jy, jz) =
                2. * nesheath * tisheath - Pi(r.ind, mesh->ystart, jz);

            Vi(r.ind, jy, jz) = 2. * visheath - Vi(r.ind, mesh->ystart, jz);
            Ve(r.ind, jy, jz) = 2. * vesheath - Ve(r.ind, mesh->ystart, jz);
            Jpar(r.ind, jy, jz) = 2. * jsheath - Jpar(r.ind, mesh->ystart, jz);
            NVi(r.ind, jy, jz) =
                2. * nesheath * visheath - NVi(r.ind, mesh->ystart, jz);
          }
        }
      }
      break;
    }
    case 3: { // Insulating Bohm sheath with free density
      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          // Free density
          BoutReal nesheath = 0.5 * (3. * Ne(r.ind, mesh->ystart, jz) -
                                     Ne(r.ind, mesh->ystart + 1, jz));
          if (nesheath < 0.0)
            nesheath = 0.0;

          // Temperature at the sheath entrance
          BoutReal tesheath = floor(Te(r.ind, mesh->ystart, jz), 0.0);
          BoutReal tisheath = floor(Ti(r.ind, mesh->ystart, jz), 0.0);

          // Zero-gradient potential
          // NOTE: This should probably not be zero-gradient,
          // since insulator could have electric field across it
          BoutReal phisheath = phi(r.ind, mesh->ystart, jz);

          // Ion velocity goes to the sound speed
          BoutReal visheath =
              -sqrt(tesheath + tisheath); // Sound speed outwards

          if (Vi(r.ind, mesh->ystart, jz) < visheath) {
            // If plasma is faster, go to plasma velocity
            visheath = Vi(r.ind, mesh->ystart, jz);
          }

          // Sheath current set to zero, as insulating boundary
          BoutReal vesheath = visheath;

          // Apply boundary condition half-way between cells
          for (int jy = mesh->ystart - 1; jy >= 0; jy--) {
            Ne(r.ind, jy, jz) = 2. * nesheath - Ne(r.ind, mesh->ystart, jz);
            phi(r.ind, jy, jz) = 2. * phisheath - phi(r.ind, mesh->ystart, jz);
            Vort(r.ind, jy, jz) = Vort(r.ind, mesh->ystart, jz);

            // Here zero-gradient Te, heat flux applied later
            Te(r.ind, jy, jz) = Te(r.ind, mesh->ystart, jz);
            Ti(r.ind, jy, jz) = Ti(r.ind, mesh->ystart, jz);

            Pe(r.ind, jy, jz) = Pe(r.ind, mesh->ystart, jz);
            Pi(r.ind, jy, jz) = Pi(r.ind, mesh->ystart, jz);

            // Dirichlet conditions
            Vi(r.ind, jy, jz) = 2. * visheath - Vi(r.ind, mesh->ystart, jz);
            Ve(r.ind, jy, jz) = 2. * vesheath - Ve(r.ind, mesh->ystart, jz);
            Jpar(r.ind, jy, jz) = -Jpar(r.ind, mesh->ystart, jz);
            NVi(r.ind, jy, jz) =
                2. * nesheath * visheath - NVi(r.ind, mesh->ystart, jz);
          }
        }
      }
      break;
    }
    }
  } else {
    // sheath_ydown == false
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        for (int jy = mesh->ystart - 1; jy >= 0; jy--) {
          // Zero-gradient Te
          Te(r.ind, jy, jz) = Te(r.ind, mesh->ystart, jz);
          Telim(r.ind, jy, jz) = Telim(r.ind, mesh->ystart, jz);

          Ti(r.ind, jy, jz) = Ti(r.ind, mesh->ystart, jz);
          Tilim(r.ind, jy, jz) = Tilim(r.ind, mesh->ystart, jz);

          Pe(r.ind, jy, jz) = Ne(r.ind, jy, jz) * Te(r.ind, jy, jz);
          Pelim(r.ind, jy, jz) = Nelim(r.ind, jy, jz) * Telim(r.ind, jy, jz);

          Pi(r.ind, jy, jz) = Ne(r.ind, jy, jz) * Ti(r.ind, jy, jz);
          Pilim(r.ind, jy, jz) = Nelim(r.ind, jy, jz) * Tilim(r.ind, jy, jz);

          // No flows
          Vi(r.ind, jy, jz) = -Vi(r.ind, mesh->ystart, jz);
          NVi(r.ind, jy, jz) = -NVi(r.ind, mesh->ystart, jz);
          Ve(r.ind, jy, jz) = -Ve(r.ind, mesh->ystart, jz);
          Jpar(r.ind, jy, jz) = -Jpar(r.ind, mesh->ystart, jz);
        }
      }
    }
  }

  if (sheath_yup) {
    switch (sheath_model) {
    case 0: { // Normal Bohm sheath
      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          // Zero-gradient density
          BoutReal nesheath = floor(Ne(r.ind, mesh->yend, jz), 0.0);

          // Temperature at the sheath entrance
          BoutReal tesheath = floor(Te(r.ind, mesh->yend, jz), 0.0);
          BoutReal tisheath = floor(Ti(r.ind, mesh->yend, jz), 0.0);

          // Zero-gradient potential
          BoutReal phisheath = phi(r.ind, mesh->yend, jz);

          // Ion velocity goes to the sound speed
          BoutReal visheath = sqrt(tesheath + tisheath); // Sound speed outwards

          if (Vi(r.ind, mesh->yend, jz) > visheath) {
            // If plasma is faster, go to plasma velocity
            visheath = Vi(r.ind, mesh->yend, jz);
          }

          // Sheath current
          // Note that phi/Te >= 0.0 since for phi < 0
          // vesheath is the electron saturation current
          BoutReal phi_te =
              floor(phisheath / Telim(r.ind, mesh->yend, jz), 0.0);
          BoutReal vesheath =
              sqrt(tesheath) * (sqrt(mi_me) / (2. * sqrt(PI))) * exp(-phi_te);
          // J = n*(Vi - Ve)
          BoutReal jsheath = nesheath * (visheath - vesheath);
          if (nesheath < 1e-10) {
            vesheath = visheath;
            jsheath = 0.0;
          }

          // Apply boundary condition half-way between cells
          for (int jy = mesh->yend + 1; jy < mesh->LocalNy; jy++) {
            // Neumann conditions
            Ne(r.ind, jy, jz) = nesheath;
            phi(r.ind, jy, jz) = phisheath;
            Vort(r.ind, jy, jz) = Vort(r.ind, mesh->yend, jz);

            // Here zero-gradient Te, heat flux applied later
            Te(r.ind, jy, jz) = Te(r.ind, mesh->yend, jz);
            Ti(r.ind, jy, jz) = Ti(r.ind, mesh->yend, jz);

            Pe(r.ind, jy, jz) = Pe(r.ind, mesh->yend, jz);
            Pelim(r.ind, jy, jz) = Pelim(r.ind, mesh->yend, jz);
            Pi(r.ind, jy, jz) = Pi(r.ind, mesh->yend, jz);
            Pilim(r.ind, jy, jz) = Pilim(r.ind, mesh->yend, jz);

            // Dirichlet conditions
            Vi(r.ind, jy, jz) = 2. * visheath - Vi(r.ind, mesh->yend, jz);
            Ve(r.ind, jy, jz) = 2. * vesheath - Ve(r.ind, mesh->yend, jz);
            Jpar(r.ind, jy, jz) = 2. * jsheath - Jpar(r.ind, mesh->yend, jz);
            NVi(r.ind, jy, jz) =
                2. * nesheath * visheath - NVi(r.ind, mesh->yend, jz);
          }
        }
      }
      break;
    }
    case 1: {
      /*
        Loizu boundary conditions

        Temperature

        Grad_par(Te) = 0

        Density equation

        Grad_par(n) = -(n/Cs) Grad_par(Vi)
        -> n_p - n_m = - (n_p + n_m)/(2Cs) (Vi_p - Vi_m)

        Pressure

        Grad_par(Pe) = Te Grad_par(n)

        Potential

        Grad_par(phi) = -Cs Grad_par(Vi)
        -> phi_p - phi_m = -Cs (Vi_p - Vi_m)
       */

      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          // Temperature at the sheath entrance
          BoutReal tesheath = floor(Te(r.ind, mesh->yend, jz), 0.0);

          // Zero gradient Te
          Te(r.ind, mesh->yend + 1, jz) = Te(r.ind, mesh->yend, jz);
          BoutReal Cs = sqrt(tesheath); // Sound speed

          // Ion velocity goes to the sound speed
          // Dirichlet boundary condition
          Vi(r.ind, mesh->yend + 1, jz) = 2. * Cs - Vi(r.ind, mesh->yend, jz);

          BoutReal g = 0.0;
          if (tesheath > 0.1 / Tnorm) {
            // Only divide by Cs if the temperature is greater than 0.1eV
            // to avoid divide-by-zero errors
            g = (Vi(r.ind, mesh->yend + 1, jz) - Vi(r.ind, mesh->yend, jz)) /
                (2. * Cs);
          }
	  
          // Mixed boundary condition on n
          Ne(r.ind, mesh->yend + 1, jz) =
              Ne(r.ind, mesh->yend, jz) * (1 - g) / (1 + g);

          // Make sure nesheath doesn't go negative
          Ne(r.ind, mesh->yend + 1, jz) =
              floor(Ne(r.ind, mesh->yend + 1, jz), -Ne(r.ind, mesh->yend, jz));
          // Density at the sheath
          BoutReal nesheath =
              0.5 * (Ne(r.ind, mesh->yend, jz) + Ne(r.ind, mesh->yend + 1, jz));

          // Momentum
          NVi(r.ind, mesh->yend + 1, jz) = NVi(r.ind, mesh->yend, jz);
          if (NVi(r.ind, mesh->yend + 1, jz) < 0.0) {
            // Limit flux to be >= 0
            NVi(r.ind, mesh->yend + 1, jz) = -NVi(r.ind, mesh->yend, jz);
          }

          // Pressure
          Pe(r.ind, mesh->yend + 1, jz) =
              Pe(r.ind, mesh->yend, jz) +
              tesheath *
                  (Ne(r.ind, mesh->yend + 1, jz) - Ne(r.ind, mesh->yend, jz));

          // Potential
          phi(r.ind, mesh->yend + 1, jz) =
              phi(r.ind, mesh->yend, jz) -
              Cs * (Vi(r.ind, mesh->yend + 1, jz) - Vi(r.ind, mesh->yend, jz));
          BoutReal phisheath = 0.5 * (phi(r.ind, mesh->yend, jz) +
                                      phi(r.ind, mesh->yend + 1, jz));

          // Sheath current
          BoutReal phi_te = floor(phisheath / tesheath, 0.0);
          BoutReal vesheath =
              Cs * (sqrt(mi_me) / (2. * sqrt(PI))) * exp(-phi_te);

          BoutReal jsheath = nesheath * (Cs - vesheath);

          // Electron velocity
          Ve(r.ind, mesh->yend + 1, jz) =
              2. * vesheath - Ve(r.ind, mesh->yend, jz);

          // Parallel velocity
          Jpar(r.ind, mesh->yend + 1, jz) =
              2. * jsheath - Jpar(r.ind, mesh->yend, jz);

          if (currents && !finite(Jpar(r.ind, mesh->yend + 1, jz))) {
            output.write("JPAR: %d, %d: %e, %e, %e, %e\n", r.ind, jz, jsheath,
                         vesheath, Cs, nesheath);
            output.write(" -> %e, %e, %e\n", Ne(r.ind, mesh->yend, jz),
                         Ne(r.ind, mesh->yend + 1, jz), g);
            exit(1);
          }

          // Constant gradient on other cells
          for (int jy = mesh->yend + 2; jy < mesh->LocalNy; jy++) {
            Vi(r.ind, jy, jz) =
                2. * Vi(r.ind, jy - 1, jz) - Vi(r.ind, jy - 2, jz);
            Ve(r.ind, jy, jz) =
                2. * Ve(r.ind, jy - 1, jz) - Ve(r.ind, jy - 2, jz);
            NVi(r.ind, jy, jz) =
                2. * NVi(r.ind, jy - 1, jz) - NVi(r.ind, jy - 2, jz);

            Ne(r.ind, jy, jz) =
                2. * Ne(r.ind, jy - 1, jz) - Ne(r.ind, jy - 2, jz);
            Te(r.ind, jy, jz) =
                2. * Te(r.ind, jy - 1, jz) - Te(r.ind, jy - 2, jz);
            Pe(r.ind, jy, jz) =
                2. * Pe(r.ind, jy - 1, jz) - Pe(r.ind, jy - 2, jz);

            phi(r.ind, jy, jz) =
                2. * phi(r.ind, jy - 1, jz) - phi(r.ind, jy - 2, jz);
            Vort(r.ind, jy, jz) =
                2. * Vort(r.ind, jy - 1, jz) - Vort(r.ind, jy - 2, jz);
            Jpar(r.ind, jy, jz) =
                2. * Jpar(r.ind, jy - 1, jz) - Jpar(r.ind, jy - 2, jz);
          }
        }
      }
      break;
    }
    case 2: { // Bohm sheath with free density
      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          // Zero-gradient density
          BoutReal nesheath = 0.5 * (3. * Ne(r.ind, mesh->yend, jz) -
                                     Ne(r.ind, mesh->yend - 1, jz));
          if (nesheath < 0.0)
            nesheath = 0.0;

          // Temperature at the sheath entrance
          BoutReal tesheath = floor(Te(r.ind, mesh->yend, jz), 0.0);
          BoutReal tisheath = floor(Ti(r.ind, mesh->yend, jz), 0.0);

          // Zero-gradient potential
          BoutReal phisheath = phi(r.ind, mesh->yend, jz);

          // Ion velocity goes to the sound speed
          BoutReal visheath = sqrt(tesheath + tisheath); // Sound speed outwards

          if (Vi(r.ind, mesh->yend, jz) > visheath) {
            // If plasma is faster, go to plasma velocity
            visheath = Vi(r.ind, mesh->yend, jz);
          }

          // Sheath current
          BoutReal phi_te =
              floor(phisheath / Telim(r.ind, mesh->yend, jz), 0.0);
          BoutReal vesheath =
              sqrt(tesheath) * (sqrt(mi_me) / (2. * sqrt(PI))) * exp(-phi_te);
          // J = n*(Vi - Ve)
          BoutReal jsheath = nesheath * (visheath - vesheath);
          if (nesheath < 1e-10) {
            vesheath = visheath;
            jsheath = 0.0;
          }

          // Apply boundary condition half-way between cells
          for (int jy = mesh->yend + 1; jy < mesh->LocalNy; jy++) {
            // Neumann conditions
            phi(r.ind, jy, jz) = phisheath;
            Vort(r.ind, jy, jz) = Vort(r.ind, mesh->yend, jz);

            // Here zero-gradient Te, heat flux applied later
            Te(r.ind, jy, jz) = tesheath;
            Ti(r.ind, jy, jz) = Ti(r.ind, mesh->yend, jz);

            // Dirichlet conditions
            Ne(r.ind, jy, jz) = 2. * nesheath - Ne(r.ind, mesh->yend, jz);
            Pe(r.ind, jy, jz) =
                2. * nesheath * tesheath - Pe(r.ind, mesh->yend, jz);
            Pi(r.ind, jy, jz) =
                2. * nesheath * tisheath - Pi(r.ind, mesh->yend, jz);

            Vi(r.ind, jy, jz) = 2. * visheath - Vi(r.ind, mesh->yend, jz);
            Ve(r.ind, jy, jz) = 2. * vesheath - Ve(r.ind, mesh->yend, jz);
            Jpar(r.ind, jy, jz) = 2. * jsheath - Jpar(r.ind, mesh->yend, jz);
            NVi(r.ind, jy, jz) =
                2. * nesheath * visheath - NVi(r.ind, mesh->yend, jz);
          }
        }
      }
      break;
    }
    case 3: { // Insulating Bohm sheath with free density
      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          // Zero-gradient density
          BoutReal nesheath = 0.5 * (3. * Ne(r.ind, mesh->yend, jz) -
                                     Ne(r.ind, mesh->yend - 1, jz));
          if (nesheath < 0.0) {
            nesheath = 0.0;
          }

          // Temperature at the sheath entrance
          BoutReal tesheath = floor(Te(r.ind, mesh->yend, jz), 0.0);
          BoutReal tisheath = floor(Ti(r.ind, mesh->yend, jz), 0.0);

          // Zero-gradient potential
          BoutReal phisheath = phi(r.ind, mesh->yend, jz);

          // Ion velocity goes to the sound speed
          BoutReal visheath = sqrt(tesheath + tisheath); // Sound speed outwards

          if (Vi(r.ind, mesh->yend, jz) > visheath) {
            // If plasma is faster, go to plasma velocity
            visheath = Vi(r.ind, mesh->yend, jz);
          }

          // Zero sheath current
          BoutReal vesheath = visheath;

          // Apply boundary condition half-way between cells
          for (int jy = mesh->yend + 1; jy < mesh->LocalNy; jy++) {
            Ne(r.ind, jy, jz) = 2. * nesheath - Ne(r.ind, mesh->yend, jz);
            phi(r.ind, jy, jz) = 2. * phisheath - phi(r.ind, mesh->yend, jz);
            Vort(r.ind, jy, jz) = Vort(r.ind, mesh->yend, jz);

            // Here zero-gradient Te, heat flux applied later
            Te(r.ind, jy, jz) = Te(r.ind, mesh->yend, jz);
            Ti(r.ind, jy, jz) = Ti(r.ind, mesh->yend, jz);

            Pe(r.ind, jy, jz) = Pe(r.ind, mesh->yend, jz);
            Pi(r.ind, jy, jz) = Pi(r.ind, mesh->yend, jz);

            // Dirichlet conditions
            Vi(r.ind, jy, jz) = 2. * visheath - Vi(r.ind, mesh->yend, jz);
            Ve(r.ind, jy, jz) = 2. * vesheath - Ve(r.ind, mesh->yend, jz);
            Jpar(r.ind, jy, jz) = -Jpar(r.ind, mesh->yend, jz);
            NVi(r.ind, jy, jz) =
                2. * nesheath * visheath - NVi(r.ind, mesh->yend, jz);
          }
        }
      }
      break;
    }
    }
  }

  //////////////////////////////////////////////////////////////
  // Test boundary conditions for UEDGE comparison
  // This applies boundary conditions to Vi, Te and Ti at the targets
  // then updates NVi, Pe and Pi boundaries
  if (test_boundaries) {
    // Lower Y target
    int jy = 1;
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        // Apply Vi = -3e4 m/s (into target)
        BoutReal vi_bndry = -3e4 / Cs0;
        Vi(r.ind, jy, jz) = 2. * vi_bndry - Vi(r.ind, jy + 1, jz);

        // Apply Te = Ti = 10eV
        BoutReal te_bndry = 10. / Tnorm;
        BoutReal ti_bndry = 10. / Tnorm;
        Te(r.ind, jy, jz) = 2. * te_bndry - Te(r.ind, jy + 1, jz);
        Ti(r.ind, jy, jz) = 2. * ti_bndry - Ti(r.ind, jy + 1, jz);

        // Get density at the boundary
        BoutReal ne_bndry = 0.5 * (Ne(r.ind, jy, jz) + Ne(r.ind, jy + 1, jz));
        if (ne_bndry < 1e-5)
          ne_bndry = 1e-5;

        NVi(r.ind, jy, jz) = 2 * ne_bndry * vi_bndry - NVi(r.ind, jy + 1, jz);
        Pe(r.ind, jy, jz) = 2 * ne_bndry * te_bndry - Pe(r.ind, jy + 1, jz);
        Pi(r.ind, jy, jz) = 2 * ne_bndry * ti_bndry - Pi(r.ind, jy + 1, jz);
      }
    }
    // Upper Y target
    jy = mesh->yend + 1;
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        // Apply Vi = 3e4 m/s
        BoutReal vi_bndry = 3e4 / Cs0;
        Vi(r.ind, jy, jz) = 2. * vi_bndry - Vi(r.ind, jy - 1, jz);

        // Apply Te = Ti = 10eV
        BoutReal te_bndry = 10. / Tnorm;
        BoutReal ti_bndry = 10. / Tnorm;
        Te(r.ind, jy, jz) = 2. * te_bndry - Te(r.ind, jy - 1, jz);
        Ti(r.ind, jy, jz) = 2. * ti_bndry - Ti(r.ind, jy - 1, jz);

        // Get density at the boundary
        BoutReal ne_bndry = 0.5 * (Ne(r.ind, jy, jz) + Ne(r.ind, jy - 1, jz));
        if (ne_bndry < 1e-5)
          ne_bndry = 1e-5;

        NVi(r.ind, jy, jz) = 2 * ne_bndry * vi_bndry - NVi(r.ind, jy - 1, jz);
        Pe(r.ind, jy, jz) = 2 * ne_bndry * te_bndry - Pe(r.ind, jy - 1, jz);
        Pi(r.ind, jy, jz) = 2 * ne_bndry * ti_bndry - Pi(r.ind, jy - 1, jz);
      }
    }
  }

  // Shift from field aligned

  Ne = fromFieldAligned(Ne);
  
  Te = fromFieldAligned(Te);
  Ti = fromFieldAligned(Ti);
  Telim = fromFieldAligned(Telim);
  Tilim = fromFieldAligned(Tilim);
  
  Pe = fromFieldAligned(Pe);
  Pi = fromFieldAligned(Pi);
  Pelim = fromFieldAligned(Pelim);
  Pilim = fromFieldAligned(Pilim);
  
  Ve = fromFieldAligned(Ve);
  Vi = fromFieldAligned(Vi);
  NVi = fromFieldAligned(NVi);

  phi = fromFieldAligned(phi);
  Jpar = fromFieldAligned(Jpar);
  Vort = fromFieldAligned(Vort);


  if (!currents) {
    // No currents, so reset Ve to be equal to Vi
    // VePsi also reset, so saved in restart file correctly
    Ve = Vi;
    VePsi = Ve;
  }

  // Ensure that Nelim, Telim and Pelim are calculated in guard cells
  Nelim = floor(Ne, 1e-5);
  Telim = floor(Te, 0.1 / Tnorm);
  Pelim = Telim * Nelim;
  Tilim = floor(Ti, 0.1 / Tnorm);
  Pilim = Tilim * Nelim;
  
  //////////////////////////////////////////////////////////////
  // Fix profiles on lower Y in SOL region by applying
  // a Dirichlet boundary condition.
  // This is to remain consistent with no-flow boundary conditions
  // on velocity fields, and to avoid spurious fluxes of energy
  // through the boundaries.
  //
  // A generator is used (sol_ne, sol_te), and where sol_ne gives a negative
  // value, no boundary condition is applied. This is to allow
  // parts of the domain to be Dirichlet, and parts (e.g. PF) to be Neumann

  if (sol_fix_profiles) {
    TRACE("Fix profiles");
    for (RangeIterator idwn = mesh->iterateBndryLowerY(); !idwn.isDone();
         idwn.next()) {

      BoutReal xnorm = mesh->GlobalX(idwn.ind);
      BoutReal ynorm =
          0.5 * (mesh->GlobalY(mesh->ystart) + mesh->GlobalY(mesh->ystart - 1));

      BoutReal neval = sol_ne->generate(xnorm, TWOPI * ynorm, 0.0, t);
      BoutReal teval = sol_te->generate(xnorm, TWOPI * ynorm, 0.0, t);

      if ((neval < 0.0) || (teval < 0.0))
        continue; // Skip, leave as previous boundary

      for (int jy = mesh->ystart - 1; jy >= 0; jy--) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          Ne(idwn.ind, jy, jz) = 2. * neval - Ne(idwn.ind, mesh->ystart, jz);
          Te(idwn.ind, jy, jz) = 2. * teval - Te(idwn.ind, mesh->ystart, jz);

          Pe(idwn.ind, jy, jz) = Ne(idwn.ind, jy, jz) * Te(idwn.ind, jy, jz);

          Telim(idwn.ind, jy, jz) = floor(Te(idwn.ind, jy, jz), 0.1 / Tnorm);

          // Zero gradient on Vi to allow flows through boundary
          Vi(idwn.ind, jy, jz) = Vi(idwn.ind, mesh->ystart, jz);

          NVi(idwn.ind, jy, jz) = Vi(idwn.ind, jy, jz) * Ne(idwn.ind, jy, jz);

          // At the boundary, Ve = Vi so no currents
          Ve(idwn.ind, jy, jz) =
              2. * Vi(idwn.ind, jy, jz) - Ve(idwn.ind, mesh->ystart, jz);
        }
      }
    }
  }

  //////////////////////////////////////////////////////////////
  // Plasma quantities calculated.
  // At this point we have calculated all boundary conditions,
  // and auxilliary variables like jpar, phi, psi
  
  // Now can update the neutral gas model, so that we can use neutral
  // gas densities and interaction frequencies in the collision processes.

  // Neutral gas
  if (neutrals) {
    TRACE("Neutral gas model update");
    
    // Update neutral gas model
    neutrals->update(Ne, Te, Ti, Vi);
  }

  //////////////////////////////////////////////////////////////
  // Collisions and stress tensor
  TRACE("Collisions");
  
  // Normalised electron collision time
  tau_e = (Cs0 / rho_s0) * tau_e0 * pow(Telim, 1.5) / Nelim;

  // Normalised ion-ion collision time
  tau_i = (Cs0 / rho_s0) * tau_i0 * pow(Tilim, 1.5) / Nelim;

  if (ion_neutral && ( neutrals || (ion_neutral_rate > 0.0))) {
    // Include ion-neutral collisions in collision time
    
    Field3D neutral_rate = neutrals ? neutrals->Fperp : ion_neutral_rate;
    
    // Add collision frequencies (1/tau_i + neutral rate)
    tau_i /= 1.0 + tau_i * neutral_rate;
  }

  // Collisional damping (normalised)
  if (resistivity || (!electromagnetic && !FiniteElMass)) {
    // Need to calculate nu if electrostatic and zero electron mass
    nu = resistivity_multiply / (1.96 * tau_e * mi_me);

    if (electron_neutral && neutrals) {
      /*
       * Include electron-neutral collisions. These can dominate
       * the resistivity at low temperatures (~1eV)
       *
       * This assumes a fixed cross-section, independent of energy
       * 
       */
      BoutReal a0 = PI*SQ(5.29e-11); // Cross-section [m^2]
      
      // Electron thermal speed
      Field3D vth_e = sqrt(mi_me*Telim);
      
      // Electron-neutral collision rate
      Field3D nu_ne = vth_e * Nnorm * neutrals->getDensity() * a0 * rho_s0;
      
      // Add collision rate to the electron-ion rate
      nu += nu_ne;
    }

    if (resistivity_boundary_width > 0 &&
        ((mesh->getGlobalXIndex(mesh->xstart) - mesh->xstart)
         < resistivity_boundary_width)) {
      // A region near the boundary with different resistivity
      
      int imax = mesh->xstart + resistivity_boundary_width - 1 -
        (mesh->getGlobalXIndex(mesh->xstart) - mesh->xstart);
      if (imax > mesh->xend) {
        imax = mesh->xend;
      }
      int imin = mesh->xstart;
      if (!mesh->firstX()) {
        --imin; // Calculate in guard cells, for radial fluxes
      }
      for (int i = imin; i <= imax; ++i) {
        for (int j = mesh->ystart; j <= mesh->yend; ++j) {
          for (int k = 0; k < mesh->LocalNz; ++k) {
            nu(i, j, k) = resistivity_boundary;
          }
        }
      }
    }
    
    // Number of points in outer guard cells
    int nguard = mesh->LocalNx - mesh->xend - 1;
    
    if (resistivity_boundary_width > 0 &&
        (mesh->GlobalNx - nguard - mesh->getGlobalXIndex(mesh->xend) <=
         resistivity_boundary_width)) {
      // Outer boundary
      int imin =
        mesh->GlobalNx - nguard - resistivity_boundary_width - mesh->getGlobalXIndex(0);
      if (imin < mesh->xstart) {
        imin = mesh->xstart;
      }
      for (int i = imin; i <= mesh->xend; ++i) {
        for (int j = mesh->ystart; j <= mesh->yend; ++j) {
          for (int k = 0; k < mesh->LocalNz; ++k) {
            nu(i, j, k) = resistivity_boundary;
          }
        }
      }
    }
  }

  if (thermal_conduction || sinks) {
    // Braginskii expression for electron parallel conduction
    // kappa ~ n * v_th^2 * tau
    kappa_epar = 3.16 * mi_me * Telim * Nelim * tau_e;

    if (kappa_limit_alpha > 0.0) {
      /*
       * Flux limiter, as used in SOLPS.
       *
       * Calculate the heat flux from Spitzer-Harm and flux limit
       *
       * Typical value of alpha ~ 0.2 for electrons
       *
       * R.Schneider et al. Contrib. Plasma Phys. 46, No. 1-2, 3 – 191 (2006)
       * DOI 10.1002/ctpp.200610001
       */

      // Spitzer-Harm heat flux
      Field3D q_SH = kappa_epar * Grad_par(Te);
      Field3D q_fl = kappa_limit_alpha * Nelim * Telim * sqrt(mi_me * Telim);

      kappa_epar = kappa_epar / (1. + abs(q_SH / q_fl));

      // Values of kappa on cell boundaries are needed for fluxes
      mesh->communicate(kappa_epar);
      kappa_epar.applyBoundary("dirichlet");
    }

    // Ion parallel heat conduction
    kappa_ipar = 3.9 * Tilim * Nelim * tau_i;

    // Boundary conditions on heat conduction coefficients
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
	kappa_epar(r.ind, mesh->ystart - 1, jz) = kappa_epar(r.ind, mesh->ystart, jz);
	kappa_ipar(r.ind, mesh->ystart - 1, jz) = kappa_ipar(r.ind, mesh->ystart, jz);
      }
    }
    
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
	kappa_epar(r.ind, mesh->yend + 1, jz) = kappa_epar(r.ind, mesh->yend, jz);
	kappa_ipar(r.ind, mesh->yend + 1, jz) = kappa_ipar(r.ind, mesh->yend, jz);
      }
    }  
  }

  nu.applyBoundary(t);

  if (ion_viscosity) {
    ///////////////////////////////////////////////////////////
    // Ion stress tensor. Split into
    // Pi_ci = Pi_ciperp + Pi_cipar
    //
    // In the parallel ion momentum equation the Pi_cipar term
    // is solved as a parallel diffusion, so is treated separately
    // All other terms are added to Pi_ciperp, even if they are
    // not really parallel parts

    Field3D Tifree = copy(Ti);
    Tifree.applyBoundary("free_o2");

    // For nonlinear terms, need to evaluate qipar and qi squared
    Field3D qipar = -kappa_ipar * Grad_par(Tifree);

    // Limit the maximum value of tau_i
    tau_i = ceil(tau_i, 1e4);

    // Square of total heat flux, parallel and perpendicular
    // The first Pi term cancels the parallel part of the second term
    // Doesn't include perpendicular collisional transport
    Field3D qisq =
        (SQ(kappa_ipar) - SQ((5. / 2) * Pilim)) * SQ(Grad_par(Tifree)) +
        SQ((5. / 2) * Pilim) * SQ(Grad(Tifree)); // This term includes a
                                                 // parallel component which is
                                                 // cancelled in first term
    // Perpendicular part from curvature
    Pi_ciperp =
        -0.5 * 0.96 * Pi * tau_i *
            (Curlb_B * Grad(phi + 1.61 * Ti) -
             Curlb_B * Grad(Pi) / Nelim) // q perpendicular
        +
        0.96 * tau_i * (1.42 / B32) *
            FV::Div_par_K_Grad_par(B32 * kappa_ipar, Tifree) // q parallel
        -
        0.49 * (qipar / Pilim) *
            (2.27 * Grad_par(log(Tilim)) - Grad_par(log(Pilim))) +
        0.75 * (0.2 * SQ(qipar) - 0.085 * qisq) / (Pilim * Tilim);
    
    // Parallel part
    Pi_cipar = -0.96 * Pi * tau_i *
               (2. * Grad_par(Vi) + Vi * Grad_par(log(coord->Bxy)));
    // Could also be written as:
    // Pi_cipar = -
    // 0.96*Pi*tau_i*2.*Grad_par(sqrt(coord->Bxy)*Vi)/sqrt(coord->Bxy);
    
    if (mesh->firstX()) {
      // First cells in X subject to boundary effects.
      for(int i=mesh->xstart; (i <= mesh->xend) && (i < 4); i++) {
        for (int j = mesh->ystart; j <= mesh->yend; j++) {
          for (int k = 0; k < mesh->LocalNz; k++) {
            Pi_ciperp(i, j, k) *= float(i - mesh->xstart)/4.0;
            Pi_cipar(i, j, k) *= float(i - mesh->xstart)/4.0;
          }
        }
      }
    }

    // Limit size of stress tensor components
    // If the off-diagonal components of the pressure tensor are large compared
    // to the scalar pressure, then the model is probably breaking down.
    // This can happen in very low density regions
    BOUT_FOR(i, Pi_ciperp.getRegion("RGN_ALL")) {
      if (fabs(Pi_ciperp[i]) > Pi[i]) {
        Pi_ciperp[i] = SIGN(Pi_ciperp[i]) * Pi[i];
      }
      if (fabs(Pi_cipar[i]) > Pi[i]) {
        Pi_cipar[i] = SIGN(Pi_cipar[i]) * Pi[i];
      }
    }
    
    mesh->communicateXZ(Pi_ciperp, Pi_cipar);

    Pi_ciperp.clearParallelSlices();
    Pi_cipar.clearParallelSlices();
    
    Pi_ciperp = toFieldAligned(Pi_ciperp);
    Pi_cipar = toFieldAligned(Pi_cipar);
    {
      int jy = 1;
      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
	for (int jz = 0; jz < mesh->LocalNz; jz++) {
	  Pi_ciperp(r.ind, jy, jz) = Pi_ciperp(r.ind, jy + 1, jz);
	  Pi_cipar(r.ind, jy, jz) = Pi_cipar(r.ind, jy + 1, jz);
	}
      }
      jy = mesh->yend + 1;
      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
	for (int jz = 0; jz < mesh->LocalNz; jz++) {
	  Pi_ciperp(r.ind, jy, jz) = Pi_ciperp(r.ind, jy - 1, jz);
	  Pi_cipar(r.ind, jy, jz) = Pi_cipar(r.ind, jy - 1, jz);
	}
      }
    }
    
    // Smooth Pi_ciperp along Y
    //Pi_ciperp.applyBoundary("neumann_o2");
    Field3D Pi_ciperp_orig = copy(Pi_ciperp);
    for (auto &i : Pi_ciperp.getRegion("RGN_NOBNDRY")) {
      Pi_ciperp[i] = 0.5*Pi_ciperp_orig[i] + 0.25*(Pi_ciperp_orig[i.ym()] + Pi_ciperp_orig[i.yp()]);
    }
    Pi_ciperp = fromFieldAligned(Pi_ciperp);
    Pi_cipar = fromFieldAligned(Pi_cipar);
    
    mesh->communicateXZ(Pi_ciperp);
    
    // Apply boundary conditions
    Pi_ciperp.applyBoundary("neumann_o2");
    Pi_cipar.applyBoundary("neumann_o2");

    Pi_ci = Pi_cipar + Pi_ciperp;
  }

  ///////////////////////////////////////////////////////////
  // Density
  // This is the electron density equation
  TRACE("density");

  if (currents) {
    // ExB drift, only if electric field is evolved
    ddt(Ne) = -Div_n_bxGrad_f_B_XPPM(Ne, phi, ne_bndry_flux, poloidal_flows,
                                     true); // ExB drift
  } else {
    ddt(Ne) = 0.0;
  }

  // Parallel flow
  if (parallel_flow) {
    if (currents) {
      // Parallel wave speed increased to electron sound speed
      // since electrostatic & electromagnetic waves are supported
      ddt(Ne) -= FV::Div_par(Ne, Ve, sqrt(mi_me) * sound_speed);
    } else {
      // Parallel wave speed is ion sound speed
      ddt(Ne) -= FV::Div_par(Ne, Ve, sound_speed);
    }
  }

  if (j_diamag) {
    // Diamagnetic drift, formulated as a magnetic drift
    // i.e Grad-B + curvature drift
    ddt(Ne) -= FV::Div_f_v(Ne, -Telim * Curlb_B, ne_bndry_flux);
  }

  if (ramp_mesh && (t < ramp_timescale)) {
    ddt(Ne) += NeTarget / ramp_timescale;
  }

  if (classical_diffusion) {
    // Classical perpendicular diffusion
    // The only term here comes from the resistive drift

    Dn = (Telim + Tilim) / (tau_e * mi_me * SQ(coord->Bxy));
    ddt(Ne) += FV::Div_a_Laplace_perp(Dn, Ne);
    ddt(Ne) += FV::Div_a_Laplace_perp(Ne / (tau_e * mi_me * SQ(coord->Bxy)),
                                      Ti - 0.5 * Te);
  }
  if (anomalous_D > 0.0) {
    ddt(Ne) += FV::Div_a_Laplace_perp(anomalous_D, DC(Ne));
  }

  // Source
  if (adapt_source) {
    // Add source. Ensure that sink will go to zero as Ne -> 0
    Field2D NeErr = averageY(DC(Ne) - NeTarget);

    if (core_sources) {
      // Sources only in core (periodic Y) domain
      // Try to keep NeTarget

      ddt(Sn) = 0.0;
      for (int x = mesh->xstart; x <= mesh->xend; x++) {
        if (!mesh->periodicY(x))
          continue; // Not periodic, so skip

        for (int y = mesh->ystart; y <= mesh->yend; y++) {
          Sn(x, y) -= source_p * NeErr(x, y);
          ddt(Sn)(x, y) = -source_i * NeErr(x, y);

          if (Sn(x, y) < 0.0) {
            Sn(x, y) = 0.0;
            if (ddt(Sn)(x, y) < 0.0)
              ddt(Sn)(x, y) = 0.0;
          }
        }
      }

      NeSource = Sn;
    } else {
      // core_sources = false
      NeSource = Sn * where(Sn, NeTarget, Ne);
      NeSource -= source_p * NeErr / NeTarget;

      ddt(Sn) = -source_i * NeErr;
    }

    if (source_vary_g11) {
      NeSource *= g11norm;
    }
  }
  
  ddt(Ne) += NeSource;
  
  if (low_n_diffuse) {
    // Diffusion which kicks in at very low density, in order to
    // help prevent negative density regions
    ddt(Ne) += FV::Div_par_K_Grad_par(SQ(coord->dy) * coord->g_22 * 1e-4 / Nelim, Ne);
  }
  if (low_n_diffuse_perp) {
    ddt(Ne) += Div_Perp_Lap_FV_Index(1e-4 / Nelim, Ne, ne_bndry_flux);
  }

  if (ne_hyper_z > 0.) {
    ddt(Ne) -= ne_hyper_z * SQ(SQ(coord->dz)) * D4DZ4(Ne);
  }

  ///////////////////////////////////////////////////////////
  // Vorticity
  // This is the current continuity equation

  TRACE("vorticity");

  ddt(Vort) = 0.0;

  if (currents) {
    // Only evolve vorticity if any diamagnetic or parallel currents
    // are included.

    if (j_par) {
      TRACE("Vort:j_par");
      // Parallel current
      
      // This term is central differencing so that it balances the parallel gradient
      // of the potential in Ohm's law
      ddt(Vort) += Div_par(Jpar);
    }

    if (j_diamag) {
      // Electron diamagnetic current

      // Note: This term is central differencing so that it balances
      // the corresponding compression term in the pressure equation
      ddt(Vort) += Div((Pi + Pe) * Curlb_B);
    }

    // Advection of vorticity by ExB
    if (boussinesq) {
      TRACE("Vort:boussinesq");
      // Using the Boussinesq approximation

      ddt(Vort) -= Div_n_bxGrad_f_B_XPPM(0.5 * Vort, phi, vort_bndry_flux,
                                         poloidal_flows);

      // V_ExB dot Grad(Pi)
      Field3D vEdotGradPi = bracket(phi, Pi, BRACKET_ARAKAWA);
      vEdotGradPi.applyBoundary("free_o2");
      // delp2(phi) term
      Field3D DelpPhi_2B2 = 0.5 * Delp2(phi) / SQ(coord->Bxy);
      DelpPhi_2B2.applyBoundary("free_o2");

      mesh->communicate(vEdotGradPi, DelpPhi_2B2);

      ddt(Vort) -= FV::Div_a_Laplace_perp(0.5 / SQ(coord->Bxy), vEdotGradPi);

    } else {
      // When the Boussinesq approximation is not made,
      // then the changing ion density introduces a number
      // of other terms.

      throw BoutException("Hot ion non-Boussinesq not implemented yet\n");
    }

    if (classical_diffusion) {
      TRACE("Vort:classical_diffusion");
      // Perpendicular viscosity
      Field3D mu = 0.3 * Tilim / (tau_i * SQ(coord->Bxy));
      ddt(Vort) += FV::Div_a_Laplace_perp(mu, Vort);
    }

    if (ion_viscosity) {
      TRACE("Vort:ion_viscosity");
      // Ion collisional viscosity.
      // Contains poloidal viscosity

      ddt(Vort) += Div(0.5 * Pi_ci * Curlb_B) -
                   Div_n_bxGrad_f_B_XPPM(1. / 3, Pi_ci, vort_bndry_flux);
    }

    if (anomalous_nu > 0.0) {
      TRACE("Vort:anomalous_nu");
      // Perpendicular anomalous momentum diffusion
      ddt(Vort) += FV::Div_a_Laplace_perp(anomalous_nu, DC(Vort));
    }
    
    if (ion_neutral_rate > 0.0) {
      // Sink of vorticity due to ion-neutral friction
      ddt(Vort) -= ion_neutral_rate * Vort;
    }
    
    if (x_hyper_viscos > 0) {
      // Form of hyper-viscosity to suppress zig-zags in X
      ddt(Vort) -= x_hyper_viscos * D4DX4_FV_Index(Vort);
    }
    if (y_hyper_viscos > 0) {
      // Form of hyper-viscosity to suppress zig-zags in Y
      ddt(Vort) -= y_hyper_viscos * FV::D4DY4_Index(Vort, false);
    }
    if (z_hyper_viscos > 0) {
      // Form of hyper-viscosity to suppress zig-zags in Z
      ddt(Vort) -= z_hyper_viscos * SQ(SQ(coord->dz)) * D4DZ4(Vort);
    }

    if (vort_dissipation) {
      // Adds dissipation term like in other equations
      // Maximum speed either electron sound speed or Alfven speed
      Field3D max_speed = Bnorm * coord->Bxy /
                          sqrt(SI::mu0 * AA * SI::Mp * Nnorm * Nelim) /
                          Cs0; // Alfven speed (normalised by Cs0)
      Field3D elec_sound = sqrt(mi_me) * sound_speed; // Electron sound speed
      for (auto& i : max_speed.getRegion(RGN_ALL)) {
	if (elec_sound[i] > max_speed[i]) {
	  max_speed[i] = elec_sound[i];
	}

        // Limit to 100x reference sound speed or light speed
        BoutReal lim = BOUTMIN(100., 3e8/Cs0);
        if (max_speed[i] > lim) {
          max_speed[i] = lim;
        }
      }

      ddt(Vort) -= FV::Div_par(Vort, 0.0, max_speed);
    }
  }

  ///////////////////////////////////////////////////////////
  // Ohm's law
  // VePsi = Ve - Vi + 0.5*mi_me*beta_e*psi
  TRACE("Ohm's law");

  ddt(VePsi) = 0.0;

  Field3D NelimVe = Nelim;

  if (currents && (electromagnetic || FiniteElMass)) {
    // Evolve VePsi except for electrostatic and zero electron mass case

    if (resistivity) {
      ddt(VePsi) -= mi_me * nu * (Ve - Vi);
      // External electric field 
      // ddt(VePsi) += mi_me*nu*(Jpar - Jpar0)/NelimVe;
    }

    // Parallel electric field
    if (j_par) {
      ddt(VePsi) += mi_me * Grad_parP(phi);
    }

    // Parallel electron pressure
    if (pe_par) {
      ddt(VePsi) -= mi_me * Grad_parP(Pelim) / NelimVe;
    }
    
    if (thermal_force) {
      ddt(VePsi) -= mi_me * 0.71 * Grad_parP(Te);
    }

    if (electron_viscosity) {
      // Electron parallel viscosity (Braginskii)
      Field3D ve_eta = 0.973 * mi_me * tau_e * Telim;
      
      if (eta_limit_alpha > 0.) {
        // SOLPS-style flux limiter
        // Values of alpha ~ 0.5 typically
        Field3D q_cl = ve_eta * Grad_par(Ve);           // Collisional value
        Field3D q_fl = eta_limit_alpha * Pelim * mi_me; // Flux limit

        ve_eta = ve_eta / (1. + abs(q_cl / q_fl));

        mesh->communicate(ve_eta);
        ve_eta.applyBoundary("neumann");
      }
      ddt(VePsi) += FV::Div_par_K_Grad_par(ve_eta, Ve);
    }
    
    if (FiniteElMass) {
      // Finite Electron Mass. Small correction needed to conserve energy
      ddt(VePsi) -= Vi * Grad_par(Ve - Vi); // Parallel advection
      ddt(VePsi) -= bracket(phi, Ve - Vi, BRACKET_ARAKAWA);  // ExB advection
      // Should also have ion polarisation advection here
    }
    
    if (numdiff > 0.0) {
      ddt(VePsi) += sqrt(mi_me) * numdiff * Div_par_diffusion_index(Ve);
    }

    if (hyper > 0.0) {
      ddt(VePsi) -= hyper * mi_me * nu * Delp2(Jpar) / Nelim;
    }

    if (hyperpar > 0.0) {
      ddt(VePsi) -= hyperpar * FV::D4DY4_Index(Ve - Vi);
    }
    
    if (vepsi_dissipation) {
      // Adds dissipation term like in other equations
      // Maximum speed either electron sound speed or Alfven speed
      Field3D max_speed = Bnorm * coord->Bxy /
                          sqrt(SI::mu0 * AA * SI::Mp * Nnorm * Nelim) /
                          Cs0; // Alfven speed (normalised by Cs0)
      Field3D elec_sound = sqrt(mi_me) * sound_speed; // Electron sound speed
      for (auto& i : max_speed.getRegion(RGN_ALL)) {
        // Maximum of Alfven or thermal electron speed
	if (elec_sound[i] > max_speed[i]) {
	  max_speed[i] = elec_sound[i];
	}

        // Limit to 100x reference sound speed or light speed
        BoutReal lim = BOUTMIN(100., 3e8/Cs0);
        if (max_speed[i] > lim) {
          max_speed[i] = lim;
        }
      }

      ddt(VePsi) -= FV::Div_par(Ve - Vi, 0.0, max_speed);
    }
  }

  ///////////////////////////////////////////////////////////
  // Ion velocity
  if (ion_velocity) {
    TRACE("Ion velocity");

    if (currents) {
      // ExB drift, only if electric field calculated
      ddt(NVi) = -Div_n_bxGrad_f_B_XPPM(NVi, phi, ne_bndry_flux,
                                        poloidal_flows); // ExB drift
    } else {
      ddt(NVi) = 0.0;
    }

    if (j_diamag) {
      // Magnetic drift
      ddt(NVi) -= FV::Div_f_v(NVi, Tilim * Curlb_B,
                              ne_bndry_flux); // Grad-B, curvature drift
    }

    ddt(NVi) -= FV::Div_par_fvv(Ne, Vi, sound_speed, false); //FV::Div_par(NVi, Vi, sound_speed, false);
    
    // Ignoring polarisation drift for now
    if (pe_par) {
      ddt(NVi) -= Grad_parP(Pe + Pi);
    }
    
    if (ion_viscosity) {
      TRACE("NVi:ion viscosity");
      // Poloidal flow damping

      // The parallel part is solved as a diffusion term
      ddt(NVi) += 1.28 * sqrtB *
                  FV::Div_par_K_Grad_par(Pi * tau_i / (coord->Bxy), sqrtB * Vi);

      if (currents) {
        // Perpendicular part. B32 = B^{3/2}
        // This is only included if ExB flow is included
        ddt(NVi) -= (2. / 3) * B32 * Grad_par(Pi_ciperp / B32);
      }
    }

    // Ion-neutral friction
    if (ion_neutral_rate > 0.0) {
      ddt(NVi) -= ion_neutral_rate * NVi;
    }

    if (numdiff > 0.0) {
      ddt(NVi) += numdiff * Div_par_diffusion_index(Vi);
      // ddt(NVi) += FV::Div_par_K_Grad_par(SQ(mesh->dy)*mesh->g_22*numdiff, Vi);
    }

    if (density_inflow) {
      // Particles arrive in cell at rate NeSource
      // This should come from a flow through the cell edge
      // with a flow velocity, and hence momentum

      ddt(NVi) += NeSource * (NeSource / Ne) * coord->dy * sqrt(coord->g_22);
    }

    if (classical_diffusion) {
      // Using same cross-field drift as in density equation
      ddt(NVi) += FV::Div_a_Laplace_perp(Vi * Dn, Ne);

      ddt(NVi) += FV::Div_a_Laplace_perp(NVi / (tau_e * mi_me * SQ(coord->Bxy)),
                                         Ti - 0.5 * Te);
    }

    if ((anomalous_D > 0.0) && anomalous_D_nvi) {
      ddt(NVi) += FV::Div_a_Laplace_perp(DC(Vi) * anomalous_D, DC(Ne));
    }
    
    if (anomalous_nu > 0.0) {
      ddt(NVi) += FV::Div_a_Laplace_perp(DC(Ne) * anomalous_nu, DC(Vi));
    }
    
    if (hyperpar > 0.0) {
      ddt(NVi) -= hyperpar * FV::D4DY4_Index(Vi) / mi_me;
    }

    if (low_n_diffuse) {
      // Diffusion which kicks in at very low density, in order to
      // help prevent negative density regions
      ddt(NVi) += FV::Div_par_K_Grad_par(SQ(coord->dy) * coord->g_22 * 1e-4 / Nelim, NVi);
    }
    if (low_n_diffuse_perp) {
      ddt(NVi) += Div_Perp_Lap_FV_Index(1e-4 / Nelim, NVi, ne_bndry_flux);
    }
  }

  ///////////////////////////////////////////////////////////
  // Pressure equation
  TRACE("Electron pressure");

  if (currents) {
    // Divergence of heat flux due to ExB advection
    ddt(Pe) =
        -Div_n_bxGrad_f_B_XPPM(Pe, phi, pe_bndry_flux, poloidal_flows, true);
  } else {
    ddt(Pe) = 0.0;
  }

  if (parallel_flow_p_term) {
    // Parallel flow
    if (currents) {
      // Like Ne term, parallel wave speed increased
      ddt(Pe) -= FV::Div_par(Pe, Ve, sqrt(mi_me) * sound_speed);
    } else {
      ddt(Pe) -= FV::Div_par(Pe, Ve, sound_speed);
    }
  }

  if (j_diamag) { // Diamagnetic flow
    // Magnetic drift (curvature) divergence.
    ddt(Pe) -= (5. / 3) * FV::Div_f_v(Pe, -Telim * Curlb_B, pe_bndry_flux);

    // This term energetically balances diamagnetic term
    // in the vorticity equation
    ddt(Pe) -= (2. / 3) * Pe * (Curlb_B * Grad(phi));
  }

  // Parallel heat conduction
  if (thermal_conduction) {
    ddt(Pe) += (2. / 3) * FV::Div_par_K_Grad_par(kappa_epar, Te);
  }

  if (thermal_flux) {
    // Parallel heat convection
    ddt(Pe) += (2. / 3) * 0.71 * Div_parP(Te * Jpar);
  }

  if (currents && resistivity) {
    // Ohmic heating
    ddt(Pe) += nu * Jpar * (Jpar - Jpar0) / Nelim;
  }

  if (pe_hyper_z > 0.0) {
    ddt(Pe) -= pe_hyper_z * SQ(SQ(coord->dz)) * D4DZ4(Pe);
  }

  ///////////////////////////////////
  // Heat transmission through sheath
  // Note: Have to calculate in field-aligned coordinates

  Field3D Ne_FA = toFieldAligned(Ne);
  Field3D Te_FA = toFieldAligned(Te);
  Field3D Ti_FA = toFieldAligned(Ti);

  wall_power = 0.0; // Diagnostic output
  if (sheath_yup) {
    TRACE("electron sheath yup heat transmission");

    Field3D sheath_dpe{zeroFrom(Ne_FA)}; // Field aligned

    switch (sheath_model) {
    case 0:
    case 2:
    case 3: {
      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          // Temperature and density at the sheath entrance
          BoutReal tesheath = floor(
              0.5 * (Te_FA(r.ind, mesh->yend, jz) + Te_FA(r.ind, mesh->yend + 1, jz)),
              0.0);
          BoutReal tisheath = floor(
              0.5 * (Ti_FA(r.ind, mesh->yend, jz) + Ti_FA(r.ind, mesh->yend + 1, jz)),
              0.0);
          BoutReal nesheath = floor(
              0.5 * (Ne_FA(r.ind, mesh->yend, jz) + Ne_FA(r.ind, mesh->yend + 1, jz)),
              0.0);

          // Sound speed (normalised units)
          BoutReal Cs = sqrt(tesheath + tisheath);

          // Heat flux
          BoutReal q = (sheath_gamma_e - 1.5) * tesheath * nesheath * Cs;

          // Multiply by cell area to get power
          BoutReal flux = q * (coord->J(r.ind, mesh->yend) +
                               coord->J(r.ind, mesh->yend + 1)) /
                          (sqrt(coord->g_22(r.ind, mesh->yend)) +
                           sqrt(coord->g_22(r.ind, mesh->yend + 1)));

          // Divide by volume of cell, and 2/3 to get pressure
          BoutReal power =
              flux / (coord->dy(r.ind, mesh->yend) * coord->J(r.ind, mesh->yend));
          sheath_dpe(r.ind, mesh->yend, jz) = -(2. / 3) * power;
          wall_power(r.ind, mesh->yend) += power;
        }
      }
      break;
    }
    }
    ddt(Pe) += fromFieldAligned(sheath_dpe);
  }
  if (sheath_ydown) {
    TRACE("electron sheath ydown heat transmission");

    Field3D sheath_dpe{zeroFrom(Te_FA)};

    switch (sheath_model) {
    case 0:
    case 2:
    case 3: {
      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          // Temperature and density at the sheath entrance
          BoutReal tesheath = floor(0.5 * (Te_FA(r.ind, mesh->ystart, jz) +
                                           Te_FA(r.ind, mesh->ystart - 1, jz)),
                                    0.0);
          BoutReal tisheath = floor(0.5 * (Ti_FA(r.ind, mesh->ystart, jz) +
                                           Ti_FA(r.ind, mesh->ystart - 1, jz)),
                                    0.0);
          BoutReal nesheath = floor(0.5 * (Ne_FA(r.ind, mesh->ystart, jz) +
                                           Ne_FA(r.ind, mesh->ystart - 1, jz)),
                                    0.0);

          // Sound speed (normalised units)
          BoutReal Cs = sqrt(tesheath + tisheath);

          // Heat flux
          BoutReal q =
              (sheath_gamma_e - 1.5) * tesheath * nesheath * Cs; // NB: positive

          // Multiply by cell area to get power
          BoutReal flux = q * (coord->J(r.ind, mesh->ystart) +
                               coord->J(r.ind, mesh->ystart - 1)) /
                          (sqrt(coord->g_22(r.ind, mesh->ystart)) +
                           sqrt(coord->g_22(r.ind, mesh->ystart - 1)));

          // Divide by volume of cell, and 2/3 to get pressure
          BoutReal power = flux / (coord->dy(r.ind, mesh->ystart) *
                                   coord->J(r.ind, mesh->ystart));
          sheath_dpe(r.ind, mesh->ystart, jz) = -(2. / 3) * power;
          wall_power(r.ind, mesh->ystart) += power;
        }
      }
      break;
    }
    }

    ddt(Pe) += fromFieldAligned(sheath_dpe);
  }
  // Transfer and source terms
  if (thermal_force) {
    ddt(Pe) -= (2. / 3) * 0.71 * Jpar * Grad_parP(Te);
  }

  if (pe_par_p_term) {
    // This term balances energetically the pressure term
    // in Ohm's law
    ddt(Pe) -= (2. / 3) * Pelim * Div_par(Ve);
  }
  if (ramp_mesh && (t < ramp_timescale)) {
    ddt(Pe) += PeTarget / ramp_timescale;
  }

  //////////////////////
  // Classical diffusion

  if (classical_diffusion) {

    // Combined resistive drift and cross-field heat diffusion
    // nu_rho2 = nu_ei * rho_e^2 in normalised units
    Field3D nu_rho2 = Telim / (tau_e * mi_me * SQ(coord->Bxy));

    ddt(Pe) +=
        (2. / 3) * (FV::Div_a_Laplace_perp(nu_rho2, Pe + Pi) +
                    (11. / 12) * FV::Div_a_Laplace_perp(nu_rho2 * Ne, Te));
  }

  //////////////////////
  // Anomalous diffusion

  if ((anomalous_D > 0.0) && anomalous_D_pepi) {
    ddt(Pe) += FV::Div_a_Laplace_perp(anomalous_D * DC(Te), DC(Ne));
  }
  if (anomalous_chi > 0.0) {
    ddt(Pe) += (2. / 3) * FV::Div_a_Laplace_perp(anomalous_chi * DC(Ne), DC(Te));
  }

  //////////////////////
  // Sources

  if (adapt_source) {
    // Add source. Ensure that sink will go to zero as Pe -> 0
    Field2D PeErr = averageY(DC(Pe) - PeTarget);

    if (core_sources) {
      // Sources only in core

      ddt(Spe) = 0.0;
      for (int x = mesh->xstart; x <= mesh->xend; x++) {
        if (!mesh->periodicY(x))
          continue; // Not periodic, so skip

        for (int y = mesh->ystart; y <= mesh->yend; y++) {
          Spe(x, y) -= source_p * PeErr(x, y);
          ddt(Spe)(x, y) = -source_i * PeErr(x, y);

          if (Spe(x, y) < 0.0) {
            Spe(x, y) = 0.0;
            if (ddt(Spe)(x, y) < 0.0)
              ddt(Spe)(x, y) = 0.0;
          }
        }
      }

      if (energy_source) {
        // Add the same amount of energy to each particle
        PeSource = Spe * Nelim / DC(Nelim);
      } else {
        PeSource = Spe;
      }
    } else {

      Spe -= source_p * PeErr / PeTarget;
      ddt(Spe) = -source_i * PeErr;

      if (energy_source) {
        // Add the same amount of energy to each particle
        PeSource = Spe * Nelim / DC(Nelim);
      } else {
        PeSource = Spe * where(Spe, PeTarget, Pe);
      }
    }
    
    if (source_vary_g11) {
      PeSource *= g11norm;
    }
    
  } else {
    // Not adapting sources

    if (energy_source) {
      // Add the same amount of energy to each particle
      PeSource = Spe * Nelim / DC(Nelim);
      
      if (source_vary_g11) {
        PeSource *= g11norm;
      }
    } else {
      // Add the same amount of energy per volume
      // If no particle source added, then this can lead to
      // a small number of particles with a lot of energy!
    }
  }

  ddt(Pe) += PeSource;
  
  ///////////////////////////////////////////////////////////
  // Ion pressure equation
  // Similar to electron pressure equation
  TRACE("Ion pressure");

  if (currents) {
    // ExB advection
    ddt(Pi) =
        -Div_n_bxGrad_f_B_XPPM(Pi, phi, pe_bndry_flux, poloidal_flows, true);
  } else {
    ddt(Pi) = 0.0;
  }

  // Parallel flow
  if (parallel_flow_p_term) {
    ddt(Pi) -= FV::Div_par(Pi, Vi, sound_speed);
  }

  if (j_diamag) { // Diamagnetic flow
    // Magnetic drift (curvature) divergence
    ddt(Pi) -= (5. / 3) * FV::Div_f_v(Pi, Tilim * Curlb_B, pe_bndry_flux);

    // Compression of ExB flow
    // These terms energetically balances diamagnetic term
    // in the vorticity equation
    ddt(Pi) -= (2. / 3) * Pi * (Curlb_B * Grad(phi));

    ddt(Pi) += Pi * Div((Pe + Pi) * Curlb_B);
  }

  if (j_par) {
    if (boussinesq) {
      ddt(Pi) -= (2. / 3) * Jpar * Grad_parP(Pi);
    } else {
      ddt(Pi) -= (2. / 3) * Jpar * Grad_parP(Pi) / Nelim;
    }
  }

  // Parallel heat conduction
  if (thermal_conduction) {
    ddt(Pi) += (2. / 3) * FV::Div_par_K_Grad_par(kappa_ipar, Ti);
  }

  // Parallel pressure gradients (sound waves)
  if (pe_par_p_term) {
    // This term balances energetically the pressure term
    // in the parallel momentum equation
    ddt(Pi) -= (2. / 3) * Pilim * Div_par(Vi);
  }
  
  if (electron_ion_transfer) {
    // Electron-ion heat transfer
    Wi = (3. / mi_me) * Nelim * (Te - Ti) / tau_e;
    ddt(Pi) += (2. / 3) * Wi;
    ddt(Pe) -= (2. / 3) * Wi;
  }
  
  //////////////////////
  // Classical diffusion

  if (classical_diffusion) {
    // Cross-field heat conduction
    // kappa_perp = 2 * n * nu_ii * rho_i^2

    ddt(Pi) += (2. / 3) * FV::Div_a_Laplace_perp(
                              2. * Pilim / (SQ(coord->Bxy) * tau_i), Ti);

    // Resistive drift terms

    // nu_rho2 = (Ti/Te) * nu_ei * rho_e^2 in normalised units
    Field3D nu_rho2 = Tilim / (tau_e * mi_me * SQ(coord->Bxy));

    ddt(Pi) += (5. / 3) *
               (FV::Div_a_Laplace_perp(nu_rho2, Pe + Pi) -
                (3. / 2) * FV::Div_a_Laplace_perp(nu_rho2 * Ne, Te));

    // Collisional heating from perpendicular viscosity
    // in the vorticity equation

    if (currents) {
      Vector3D Grad_perp_vort = Grad(Vort);
      Grad_perp_vort.y = 0.0; // Zero parallel component
      ddt(Pi) -= (2. / 3) * (3. / 10) * Tilim / (SQ(coord->Bxy) * tau_i) *
                 (Grad_perp_vort * Grad(phi + Pi));
    }
  }

  if (ion_viscosity) {
    // Collisional heating due to parallel viscosity
    
    ddt(Pi) += (2. / 3) * 1.28 * (Pi * tau_i / sqrtB) * Grad_par(sqrtB * Vi) * Div_par(Vi);
    
    if (currents) {
      ddt(Pi) -= (4. / 9) * Pi_ciperp * Div_par(Vi);
      //(4. / 9) * Vi * B32 * Grad_par(Pi_ciperp / B32);

      ddt(Pi) -= (2. / 6) * Pi_ci * Curlb_B * Grad(phi + Pi);
      ddt(Pi) += (2. / 9) * bracket(Pi_ci, phi + Pi, BRACKET_ARAKAWA);
    }
  }
  
  //////////////////////
  // Anomalous diffusion

  if ((anomalous_D > 0.0) && anomalous_D_pepi) {
    ddt(Pi) += FV::Div_a_Laplace_perp(anomalous_D * DC(Ti), DC(Ne));
  }

  if (anomalous_chi > 0.0) {
    ddt(Pi) +=
        (2. / 3) * FV::Div_a_Laplace_perp(anomalous_chi * DC(Ne), DC(Ti));
  }

  ///////////////////////////////////
  // Heat transmission through sheath

  if (sheath_yup) {
    TRACE("ion sheath yup heat transmission");
    
    Field3D sheath_dpi{zeroFrom(Te_FA)}; // Field aligned

    switch (sheath_model) {
    case 0:
    case 2:
    case 3: {
      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          // Temperature and density at the sheath entrance
          BoutReal tesheath = floor(
              0.5 * (Te_FA(r.ind, mesh->yend, jz) + Te_FA(r.ind, mesh->yend + 1, jz)),
              0.0);
          BoutReal tisheath = floor(
              0.5 * (Ti_FA(r.ind, mesh->yend, jz) + Ti_FA(r.ind, mesh->yend + 1, jz)),
              0.0);
          BoutReal nesheath = floor(
              0.5 * (Ne_FA(r.ind, mesh->yend, jz) + Ne_FA(r.ind, mesh->yend + 1, jz)),
              0.0);

          // Sound speed (normalised units)
          BoutReal Cs = sqrt(tesheath + tisheath);

          // Heat flux
          BoutReal q = (sheath_gamma_i - 1.5) * tisheath * nesheath * Cs;

          // Multiply by cell area to get power
          BoutReal flux = q * (coord->J(r.ind, mesh->yend) +
                               coord->J(r.ind, mesh->yend + 1)) /
                          (sqrt(coord->g_22(r.ind, mesh->yend)) +
                           sqrt(coord->g_22(r.ind, mesh->yend + 1)));

          // Divide by volume of cell, and 2/3 to get pressure
          BoutReal power =
              flux / (coord->dy(r.ind, mesh->yend) * coord->J(r.ind, mesh->yend));
          sheath_dpi(r.ind, mesh->yend, jz) = -(2. / 3) * power;
          wall_power(r.ind, mesh->yend) += power;
        }
      }
      break;
    }
    }
    ddt(Pi) += fromFieldAligned(sheath_dpi);
  }
  if (sheath_ydown) {
    TRACE("ion sheath ydown heat transmission");

    Field3D sheath_dpi{zeroFrom(Ti_FA)};

    switch (sheath_model) {
    case 0:
    case 2:
    case 3: {
      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          // Temperature and density at the sheath entrance
          BoutReal tesheath = floor(0.5 * (Te_FA(r.ind, mesh->ystart, jz) +
                                           Te_FA(r.ind, mesh->ystart - 1, jz)),
                                    0.0);
          BoutReal tisheath = floor(0.5 * (Ti_FA(r.ind, mesh->ystart, jz) +
                                           Ti_FA(r.ind, mesh->ystart - 1, jz)),
                                    0.0);
          BoutReal nesheath = floor(0.5 * (Ne_FA(r.ind, mesh->ystart, jz) +
                                           Ne_FA(r.ind, mesh->ystart - 1, jz)),
                                    0.0);

          // Sound speed (normalised units)
          BoutReal Cs = sqrt(tesheath + tisheath);

          // Heat flux (positive)
          BoutReal q = (sheath_gamma_i - 1.5) * tisheath * nesheath * Cs;

          // Multiply by cell area to get power
          BoutReal flux = q * (coord->J(r.ind, mesh->ystart) +
                               coord->J(r.ind, mesh->ystart - 1)) /
                          (sqrt(coord->g_22(r.ind, mesh->ystart)) +
                           sqrt(coord->g_22(r.ind, mesh->ystart - 1)));

          // Divide by volume of cell, and 2/3 to get pressure
          BoutReal power = flux / (coord->dy(r.ind, mesh->ystart) *
                                   coord->J(r.ind, mesh->ystart));
          sheath_dpi(r.ind, mesh->ystart, jz) = -(2. / 3) * power;
          wall_power(r.ind, mesh->ystart) += power;
        }
      }
      break;
    }
    }

    ddt(Pi) += fromFieldAligned(sheath_dpi);
  }

  //////////////////////
  // Sources

  if (adapt_source) {
    // Add source. Ensure that sink will go to zero as Pe -> 0
    Field2D PiErr = averageY(DC(Pi) - PiTarget);

    if (core_sources) {
      // Sources only in core

      ddt(Spi) = 0.0;
      for (int x = mesh->xstart; x <= mesh->xend; x++) {
        if (!mesh->periodicY(x))
          continue; // Not periodic, so skip

        for (int y = mesh->ystart; y <= mesh->yend; y++) {
          Spi(x, y) -= source_p * PiErr(x, y);
          ddt(Spi)(x, y) = -source_i * PiErr(x, y);

          if (Spi(x, y) < 0.0) {
            Spi(x, y) = 0.0;
            if (ddt(Spi)(x, y) < 0.0)
              ddt(Spi)(x, y) = 0.0;
          }
        }
      }

      if (energy_source) {
        // Add the same amount of energy to each particle
        PiSource = Spi * Nelim / DC(Nelim);
      } else {
        PiSource = Spi;
      }
    } else {

      Spi -= source_p * PiErr / PiTarget;
      ddt(Spi) = -source_i * PiErr;

      if (energy_source) {
        // Add the same amount of energy to each particle
        PiSource = Spi * Nelim / DC(Nelim);
      } else {
        PiSource = Spi * where(Spi, PiTarget, Pi);
      }
    }
    
    if (source_vary_g11) {
      PiSource *= g11norm;
    }
    
  } else {
    // Not adapting sources

    if (energy_source) {
      // Add the same amount of energy to each particle
      PiSource = Spi * Nelim / DC(Nelim);

      if (source_vary_g11) {
        PiSource *= g11norm;
      }
      
    } else {
      // Add the same amount of energy per volume
      // If no particle source added, then this can lead to
      // a small number of particles with a lot of energy!
    }
  }
  
  ddt(Pi) += PiSource;
  
  ///////////////////////////////////////////////////////////
  // Radial buffer regions for turbulence simulations

  if (radial_buffers) {
    /// Radial buffer regions

    // Calculate flux sZ averages
    Field2D PeDC = averageY(DC(Pe));
    Field2D PiDC = averageY(DC(Pi));
    Field2D NeDC = averageY(DC(Ne));
    Field2D VortDC = averageY(DC(Vort));

    if ((mesh->getGlobalXIndex(mesh->xstart) - mesh->xstart) < radial_inner_width) {
      // This processor contains points inside the inner radial boundary

      int imax = mesh->xstart + radial_inner_width - 1 -
                 (mesh->getGlobalXIndex(mesh->xstart) - mesh->xstart);
      if (imax > mesh->xend) {
        imax = mesh->xend;
      }

      int imin = mesh->xstart;
      if (!mesh->firstX()) {
        --imin; // Calculate in guard cells, for radial fluxes
      }
      int ncz = mesh->LocalNz;

      for (int i = imin; i <= imax; ++i) {
        // position inside the boundary (0 = on boundary, 0.5 = first cell)
        BoutReal pos =
            static_cast<BoutReal>(mesh->getGlobalXIndex(i) - mesh->xstart) + 0.5;

        // Diffusion coefficient which increases towards the boundary
        BoutReal D = radial_buffer_D * (1. - pos / radial_inner_width);

        for (int j = mesh->ystart; j <= mesh->yend; ++j) {
          BoutReal dx = coord->dx(i, j);
          BoutReal dx_xp = coord->dx(i + 1, j);
          BoutReal J = coord->J(i, j);
          BoutReal J_xp = coord->J(i + 1, j);

          // Calculate metric factors for radial fluxes
          BoutReal rad_flux_factor = 0.25 * (J + J_xp) * (dx + dx_xp);
          BoutReal x_factor = rad_flux_factor / (J * dx);
          BoutReal xp_factor = rad_flux_factor / (J_xp * dx_xp);

          for (int k = 0; k < ncz; ++k) {
            // Relax towards constant value on flux surface
            ddt(Pe)(i, j, k) -= D * (Pe(i, j, k) - PeDC(i, j));
            ddt(Pi)(i, j, k) -= D * (Pi(i, j, k) - PiDC(i, j));
            ddt(Ne)(i, j, k) -= D * (Ne(i, j, k) - NeDC(i, j));
            ddt(Vort)(i, j, k) -= D * (Vort(i, j, k) - VortDC(i, j));
            ddt(NVi)(i, j, k) -= D * NVi(i, j, k);
            
            // Radial fluxes
            BoutReal f = D * (Ne(i + 1, j, k) - Ne(i, j, k));
            ddt(Ne)(i, j, k) += f * x_factor;
            ddt(Ne)(i + 1, j, k) -= f * xp_factor;

            f = D * (Pe(i + 1, j, k) - Pe(i, j, k));
            ddt(Pe)(i, j, k) += f * x_factor;
            ddt(Pe)(i + 1, j, k) -= f * xp_factor;

            f = D * (Pi(i + 1, j, k) - Pi(i, j, k));
            ddt(Pi)(i, j, k) += f * x_factor;
            ddt(Pi)(i + 1, j, k) -= f * xp_factor;

            f = D * (Vort(i + 1, j, k) - Vort(i, j, k));
            ddt(Vort)(i, j, k) += f * x_factor;
            ddt(Vort)(i + 1, j, k) -= f * xp_factor;
          }
        }
      }
    }
    // Number of points in outer guard cells
    int nguard = mesh->LocalNx - mesh->xend - 1;

    if (mesh->GlobalNx - nguard - mesh->getGlobalXIndex(mesh->xend) <=
        radial_outer_width) {

      // Outer boundary
      int imin =
          mesh->GlobalNx - nguard - radial_outer_width - mesh->getGlobalXIndex(0);
      if (imin < mesh->xstart) {
        imin = mesh->xstart;
      }
      int ncz = mesh->LocalNz;
      for (int i = imin; i <= mesh->xend; ++i) {

        // position inside the boundary
        BoutReal pos =
            static_cast<BoutReal>(mesh->GlobalNx - nguard - mesh->getGlobalXIndex(i)) -
            0.5;

        // Diffusion coefficient which increases towards the boundary
        BoutReal D = radial_buffer_D * (1. - pos / radial_outer_width);

        for (int j = mesh->ystart; j <= mesh->yend; ++j) {
          BoutReal dx = coord->dx(i, j);
          BoutReal dx_xp = coord->dx(i + 1, j);
          BoutReal J = coord->J(i, j);
          BoutReal J_xp = coord->J(i + 1, j);

          // Calculate metric factors for radial fluxes
          BoutReal rad_flux_factor = 0.25 * (J + J_xp) * (dx + dx_xp);
          BoutReal x_factor = rad_flux_factor / (J * dx);
          BoutReal xp_factor = rad_flux_factor / (J_xp * dx_xp);

          for (int k = 0; k < ncz; ++k) {
            ddt(Pe)(i, j, k) -= D * (Pe(i, j, k) - PeDC(i, j));
            ddt(Pi)(i, j, k) -= D * (Pi(i, j, k) - PiDC(i, j));
            ddt(Ne)(i, j, k) -= D * (Ne(i, j, k) - NeDC(i, j));
            ddt(Vort)(i, j, k) -= D * (Vort(i, j, k) - VortDC(i, j));
            // ddt(Vort)(i,j,k) -= D*Vort(i,j,k);

            BoutReal f = D * (Vort(i + 1, j, k) - Vort(i, j, k));
            ddt(Vort)(i, j, k) += f * x_factor;
            ddt(Vort)(i + 1, j, k) -= f * xp_factor;
          }
        }
      }
    }
  }

  ///////////////////////////////////////////////////////////
  // Neutral gas
  if (neutrals) {
    TRACE("Neutral gas model add sources");
    
    // Add sources/sinks to plasma equations
    ddt(Ne) -= neutrals->S;             // Sink of plasma density
    ddt(NVi) -= neutrals->F;            // Plasma momentum
    ddt(Pe) -= (2. / 3) * neutrals->Rp; // Plasma radiated energy
    ddt(Pi) -= (2. / 3) * neutrals->Qi; // Ion energy transferred to neutrals

    // Calculate atomic rates

    if (neutral_friction) {
      // Vorticity
      if (boussinesq) {
        ddt(Vort) -= FV::Div_a_Laplace_perp(
            neutrals->Fperp / (Nelim * SQ(coord->Bxy)), phi);
      } else {
        ddt(Vort) -=
            FV::Div_a_Laplace_perp(neutrals->Fperp / SQ(coord->Bxy), phi);
      }
    }

    ///////////////////////////////////////////////////////////////////
    // Recycling at the boundary
    TRACE("Neutral recycling fluxes");
    wall_flux = 0.0;

    if (sheath_ydown) {
      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        // Calculate flux of ions into target from Ne and Vi boundary
        // This calculation is supposed to be consistent with the flow
        // of plasma from Div_par_FV(Ne, Ve)

        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          BoutReal flux_ion =
              -0.5 *
              (Ne(r.ind, mesh->ystart, jz) + Ne(r.ind, mesh->ystart - 1, jz)) *
              0.5 * (Ve(r.ind, mesh->ystart, jz) +
                     Ve(r.ind, mesh->ystart - 1, jz)); // Flux through surface
                                                       // [m^-2 s^-1], should be
                                                       // positive since Ve <
                                                       // 0.0

          // Flow of neutrals inwards
          BoutReal flow = frecycle * flux_ion *
                          (coord->J(r.ind, mesh->ystart) +
                           coord->J(r.ind, mesh->ystart - 1)) /
                          (sqrt(coord->g_22(r.ind, mesh->ystart)) +
                           sqrt(coord->g_22(r.ind, mesh->ystart - 1)));

          // Rate of change of neutrals in final cell
          BoutReal dndt = flow / (coord->J(r.ind, mesh->ystart) *
                                  coord->dy(r.ind, mesh->ystart));

          // Add mass, momentum and energy to the neutrals

          neutrals->addDensity(r.ind, mesh->ystart, jz, dndt);
          neutrals->addPressure(r.ind, mesh->ystart, jz,
                                dndt * (3.5 / Tnorm)); // Franck-Condon energy
          neutrals->addMomentum(r.ind, mesh->ystart, jz,
                                dndt * neutral_vwall * sqrt(3.5 / Tnorm));

          // Power deposited onto the wall due to surface recombination
          wall_power(r.ind, mesh->ystart) += (13.6 / Tnorm) * dndt;
        }
      }
    }

    if (sheath_yup) {
      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        // Calculate flux of ions into target from Ne and Vi boundary
        // This calculation is supposed to be consistent with the flow
        // of plasma from FV::Div_par(Ne, Ve)

        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          // Flux through surface [m^-2 s^-1], should be positive
          BoutReal flux_ion =
              frecycle * 0.5 *
              (Ne(r.ind, mesh->yend, jz) + Ne(r.ind, mesh->yend + 1, jz)) *
              0.5 * (Ve(r.ind, mesh->yend, jz) + Ve(r.ind, mesh->yend + 1, jz));

          // Flow of neutrals inwards
          BoutReal flow = flux_ion * (coord->J(r.ind, mesh->yend) +
                                      coord->J(r.ind, mesh->yend + 1)) /
                          (sqrt(coord->g_22(r.ind, mesh->yend)) +
                           sqrt(coord->g_22(r.ind, mesh->yend + 1)));

          // Rate of change of neutrals in final cell
          BoutReal dndt =
              flow / (coord->J(r.ind, mesh->yend) * coord->dy(r.ind, mesh->yend));

          // Add mass, momentum and energy to the neutrals

          neutrals->addDensity(r.ind, mesh->yend, jz, dndt);
          neutrals->addPressure(r.ind, mesh->yend, jz,
                                dndt * (3.5 / Tnorm)); // Franck-Condon energy
          neutrals->addMomentum(r.ind, mesh->yend, jz,
                                -dndt * neutral_vwall * sqrt(3.5 / Tnorm));

          // Power deposited onto the wall due to surface recombination
          wall_power(r.ind, mesh->yend) += (13.6 / Tnorm) * dndt;
        }
      }
    }
  }

  //////////////////////////////////////////////////////////////
  // Impurities

  if ((carbon_fraction > 0.0) || impurity_adas) {
    
    if (carbon_fraction > 0.0) {
      TRACE("Carbon impurity radiation");
      Rzrad = carbon_rad->power(Te * Tnorm, Ne * Nnorm,
                                Ne * (Nnorm * carbon_fraction)); // J / m^3 / s
    } else {
      Rzrad = 0.0;
    }

    if (impurity_adas) {
      for (auto &i : Rzrad.getRegion("RGN_NOY")) {
        // Calculate cell centre (C), left (L) and right (R) values
        
        BoutReal Te_C = Te[i],
          Te_L = 0.5 * (Te[i.ym()] + Te[i]),
          Te_R = 0.5 * (Te[i] + Te[i.yp()]);
        BoutReal Ne_C = Ne[i],
          Ne_L = 0.5 * (Ne[i.ym()] + Ne[i]),
          Ne_R = 0.5 * (Ne[i] + Ne[i.yp()]);
        
        BoutReal Rz_L = computeRadiatedPower(*impurity,
                                             Te_L * Tnorm,        // electron temperature [eV]
                                             Ne_L * Nnorm,        // electron density [m^-3]
                                             fimp * Ne_L * Nnorm, // impurity density [m^-3]
                                             0.0);       // Neutral density [m^-3]
        
        BoutReal Rz_C = computeRadiatedPower(*impurity,
                                             Te_C * Tnorm,        // electron temperature [eV]
                                             Ne_C * Nnorm,        // electron density [m^-3]
                                             fimp * Ne_C * Nnorm, // impurity density [m^-3]
                                             0.0);       // Neutral density [m^-3]
        
        BoutReal Rz_R = computeRadiatedPower(*impurity,
                                             Te_R * Tnorm,        // electron temperature [eV]
                                             Ne_R * Nnorm,        // electron density [m^-3]
                                             fimp * Ne_R * Nnorm, // impurity density [m^-3]
                                             0.0);       // Neutral density [m^-3]
        
        // Jacobian (Cross-sectional area)
        BoutReal J_C = coord->J[i],
          J_L = 0.5 * (coord->J[i.ym()] + coord->J[i]),
          J_R = 0.5 * (coord->J[i] + coord->J[i.yp()]);
        
        // Simpson's rule, calculate average over cell
        Rzrad[i] += (J_L * Rz_L +
                     4. * J_C * Rz_C +
                     J_R * Rz_R) / (6. * J_C);
      }
    }
    
    Rzrad /= SI::qe * Tnorm * Nnorm * Omega_ci;                // Normalise
    ddt(Pe) -= (2. / 3) * Rzrad;
  }

  //////////////////////////////////////////////////////////////
  // Parallel closures for 2D simulations

  if (sinks) {
    // Sink terms for 2D simulations

    // Field3D nsink = 0.5*Ne*sqrt(Telim)*sink_invlpar;   // n C_s/ (2L)  //
    // Sound speed flow to targets
    Field3D nsink = 0.5 * sqrt(Ti) * Ne * sink_invlpar;
    nsink = floor(nsink, 0.0);

    ddt(Ne) -= nsink;

    Field3D conduct = (2. / 3) * kappa_epar * Te * SQ(sink_invlpar);
    conduct = floor(conduct, 0.0);
    ddt(Pe) -= conduct      // Heat conduction
               + Te * nsink // Advection
        ;

    if (sheath_closure) {
      ///////////////////////////
      // Sheath dissipation closure

      Field3D phi_te = floor(phi / Telim, 0.0);
      Field3D jsheath = Nelim * sqrt(Telim) *
                        (1 - (sqrt(mi_me) / (2. * sqrt(PI))) * exp(-phi_te));

      ddt(Vort) += jsheath * sink_invlpar;
    } else {
      ///////////////////////////
      // Vorticity closure
      ddt(Vort) -= FV::Div_a_Laplace_perp(nsink / SQ(coord->Bxy), phi);
      // ddt(Vort) -= (nsink / Ne)*Vort;
      // Note: If nsink = n * nu then this reduces to
      // ddt(Vort) -= nu * Vort
    }

    if (drift_wave) {
      // Include a drift-wave closure in the core region
      // as in Hasegawa-Wakatani and SOLT models

      Field2D Tedc = DC(Telim); // Zonal average
      Field3D Q = phi - Telim * log(Nelim);

      // Drift-wave operator
      Field3D Adw = alpha_dw * pow(Tedc, 1.5) * (Q - DC(Q));

      ddt(Ne) += Adw;
      ddt(Vort) += Adw;
    }

    // Electron and ion parallel dynamics not evolved
  }

  if (low_pass_z >= 0) {
    // Low pass Z filtering, keeping up to and including low_pass_z
    ddt(Ne) = lowPass(ddt(Ne), low_pass_z);
    ddt(Pe) = lowPass(ddt(Pe), low_pass_z);

    if (currents) {
      ddt(Vort) = lowPass(ddt(Vort), low_pass_z);
      if (electromagnetic || FiniteElMass) {
        ddt(VePsi) = lowPass(ddt(VePsi), low_pass_z);
      }
    }
    if (ion_velocity) {
      ddt(NVi) = lowPass(ddt(NVi), low_pass_z);
    }
  }

  if (!evolve_plasma) {
    ddt(Ne) = 0.0;
    ddt(Pe) = 0.0;
    ddt(Vort) = 0.0;
    ddt(VePsi) = 0.0;
    ddt(NVi) = 0.0;
  }

  return 0;
}

/*!
 * Preconditioner. Solves the heat conduction
 *
 * @param[in] t  The simulation time
 * @param[in] gamma   Factor in front of the Jacobian in (I - gamma*J). Related
 * to timestep
 * @param[in] delta   Not used here
 */
int Hermes::precon(BoutReal t, BoutReal gamma, BoutReal delta) {
  static InvertPar *inv = NULL;
  if (!inv) {
    // Initialise parallel inversion class
    inv = InvertPar::Create();
    inv->setCoefA(1.0);
  }
  if (thermal_conduction) {
    // Set the coefficient in front of Grad2_par2
    inv->setCoefB(-(2. / 3) * gamma * kappa_epar);
    Field3D dT = ddt(Pe);
    dT.applyBoundary("neumann");
    ddt(Pe) = inv->solve(dT);
  }

  // Neutral gas preconditioning
  if (neutrals)
    neutrals->precon(t, gamma, delta);

  return 0;
}

const Field3D Hermes::Grad_parP(const Field3D &f) {
  return Grad_par(f); //+ 0.5*beta_e*bracket(psi, f, BRACKET_ARAKAWA);
}

const Field3D Hermes::Div_parP(const Field3D &f) {
  return Div_par(f);
  //+ 0.5*beta_e*coord->Bxy*bracket(psi, f/coord->Bxy, BRACKET_ARAKAWA);
}

// Standard main() function
BOUTMAIN(Hermes);
