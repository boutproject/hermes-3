
#include "../include/vorticity.hxx"
#include "../include/div_ops.hxx"

#include <bout/fv_ops.hxx>
#include <invert_laplace.hxx>
#include <bout/invert/laplacexy.hxx>
#include <bout/constants.hxx>
#include <difops.hxx>
#include <derivs.hxx>

using bout::globals::mesh;

Vorticity::Vorticity(std::string name, Options &alloptions, Solver *solver) {
  AUTO_TRACE();
  
  solver->add(Vort, "Vort");

  SAVE_REPEAT(phi);
  
  auto& options = alloptions[name];

  exb_advection = options["exb_advection"]
                      .doc("Include ExB advection (nonlinear term)?")
                      .withDefault<bool>(true);

  diamagnetic = options["diamagnetic"]
                    .doc("Include diamagnetic current?")
                    .withDefault<bool>(true);

  sheath_boundary =
      options["sheath_boundary"]
          .doc("Set potential to j=0 sheath at boundaries? (default = 0)")
          .withDefault<bool>(false);

  diamagnetic_polarisation =
      options["diamagnetic_polarisation"]
          .doc("Include diamagnetic drift in polarisation current?")
          .withDefault<bool>(true);

  collisional_friction =
    options["collisional_friction"]
    .doc("Damp vorticity based on mass-weighted collision frequency")
    .withDefault<bool>(false);

  average_atomic_mass =
      options["average_atomic_mass"]
          .doc("Weighted average atomic mass, for polarisaion current "
               "(Boussinesq approximation)")
          .withDefault<BoutReal>(2.0); // Deuterium

  bndry_flux = options["bndry_flux"]
                   .doc("Allow flows through radial boundaries")
                   .withDefault<bool>(true);

  poloidal_flows = options["poloidal_flows"]
                       .doc("Include poloidal ExB flow")
                       .withDefault<bool>(true);

  split_n0 = options["split_n0"]
                 .doc("Split phi into n=0 and n!=0 components")
                 .withDefault<bool>(false);

  hyper_z = options["hyper_z"]
    .doc("Hyper-viscosity in Z. < 0 -> off")
    .withDefault(-1.0);

  // Numerical dissipation terms
  // These are required to suppress parallel zig-zags in
  // cell centred formulations. Essentially adds (hopefully small)
  // parallel currents

  vort_dissipation = options["vort_dissipation"]
    .doc("Parallel dissipation of vorticity")
    .withDefault<bool>(false);

  phi_dissipation = options["phi_dissipation"]
    .doc("Parallel dissipation of potential [Recommended]")
    .withDefault<bool>(true);

  phi_boundary_relax = options["phi_boundary_relax"]
    .doc("Relax x boundaries of phi towards Neumann?")
    .withDefault<bool>(false);

  // Add phi to restart files so that the value in the boundaries
  // is restored on restart. This is done even when phi is not evolving,
  // so that phi can be saved and re-loaded
  get_restart_datafile()->addOnce(phi, "phi");

  // Set initial value. Will be overwritten if restarting
  phi = 0.0;

  auto coord = mesh->getCoordinates();

  if (split_n0) {
    // Create an XY solver for n=0 component
    laplacexy = new LaplaceXY(mesh);
    // Set coefficients for Boussinesq solve
    laplacexy->setCoefs(average_atomic_mass / SQ(coord->Bxy), 0.0);
  }
  phiSolver = Laplacian::create(&options["laplacian"]);
  // Set coefficients for Boussinesq solve
  phiSolver->setCoefC(average_atomic_mass / SQ(coord->Bxy));
  
  if (phi_boundary_relax) {
    // Set the last update time to -1, so it will reset
    // the first time RHS function is called
    phi_boundary_last_update = -1.;

    phi_boundary_timescale = options["phi_boundary_timescale"]
      .doc("Timescale for phi boundary relaxation [seconds]")
      .withDefault(1e-4)
      / get<BoutReal>(alloptions["units"]["seconds"]);
    // Normalise to internal time units

    phiSolver->setInnerBoundaryFlags(INVERT_SET);
    phiSolver->setOuterBoundaryFlags(INVERT_SET);
  }

  // Read curvature vector
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
      if (diamagnetic) {
        // Need curvature
        throw;
      } else {
        output_warn.write("No curvature vector in input grid");
        Curlb_B = 0.0;
      }
    }
  }

  if (Options::root()["mesh"]["paralleltransform"]["type"].as<std::string>() == "shifted") {
    Field2D I;
    mesh->get(I, "sinty");
    Curlb_B.z += I * Curlb_B.x;
  }
  
  Options& units = alloptions["units"];
  BoutReal Bnorm = units["Tesla"];
  BoutReal Lnorm = units["meters"];
  
  Curlb_B.x /= Bnorm;
  Curlb_B.y *= SQ(Lnorm);
  Curlb_B.z *= SQ(Lnorm);

  Curlb_B *= 2. / coord->Bxy;
  
  Bsq = SQ(coord->Bxy);
}

void Vorticity::transform(Options &state) {
  AUTO_TRACE();

  auto& fields = state["fields"];

  // Set the boundary of phi. Both 2D and 3D fields are kept, though the 3D field
  // is constant in Z. This is for efficiency, to reduce the number of conversions.
  // Note: For now the boundary values are all at the midpoint,
  //       and only phi is considered, not phi + Pi which is handled in Boussinesq solves
  Field3D Pi_sum = 0.0; ///< Contribution from ion pressure, weighted by atomic mass
  if (diamagnetic_polarisation) {
    // Diamagnetic term in vorticity. Note this is weighted by the mass
    // This includes all species, including electrons
    Options& allspecies = state["species"];
    for (auto& kv : allspecies.getChildren()) {
      Options& species = allspecies[kv.first]; // Note: need non-const

      if (!(IS_SET_NOBOUNDARY(species["pressure"]) and
            species.isSet("charge") and species.isSet("AA"))) {
        continue; // No pressure, charge or mass -> no polarisation current
      }

      // Don't need sheath boundary
      auto P = GET_NOBOUNDARY(Field3D, species["pressure"]);
      auto AA = get<BoutReal>(species["AA"]);

      Pi_sum += P * (AA / average_atomic_mass);
    }
  }

  Pi_sum.applyBoundary("neumann");

  if (phi_boundary_relax) {
    // Update the boundary regions by relaxing towards zero gradient
    // on a given timescale.

    BoutReal time = get<BoutReal>(state["time"]);

    if (phi_boundary_last_update < 0.0) {
      // First time this has been called.
      phi_boundary_last_update = time;
      
    } else if (time > phi_boundary_last_update) {
      // Only update if time has advanced
      // Uses an exponential decay of the weighting of the value in the boundary
      // so that the solution is well behaved for arbitrary steps
      BoutReal weight = exp(-(time - phi_boundary_last_update) / phi_boundary_timescale);
      phi_boundary_last_update = time;

      if (mesh->firstX()) {
        for (int j = mesh->ystart; j <= mesh->yend; j++) {
          BoutReal phivalue = 0.0;
          for (int k = 0; k < mesh->LocalNz; k++) {
            phivalue += phi(mesh->xstart, j, k);
          }
          phivalue /= mesh->LocalNz; // Average in Z of point next to boundary

          // Old value of phi at boundary
          BoutReal oldvalue =
            0.5 * (phi(mesh->xstart - 1, j, 0) + phi(mesh->xstart, j, 0));

          // New value of phi at boundary, relaxing towards phivalue
          BoutReal newvalue =
            weight * oldvalue + (1. - weight) * phivalue;

          // Set phi at the boundary to this value
          for (int k = 0; k < mesh->LocalNz; k++) {
            phi(mesh->xstart - 1, j, k) = 2.*newvalue - phi(mesh->xstart, j, k);

            // Note: This seems to make a difference, but don't know why.
            // Without this, get convergence failures with no apparent instability
            // (all fields apparently smooth, well behaved)
            phi(mesh->xstart - 2, j, k) = phi(mesh->xstart - 1, j, k);
          }
        }
      }

      if (mesh->lastX()) {
        for (int j = mesh->ystart; j <= mesh->yend; j++) {
          BoutReal phivalue = 0.0;
          for (int k = 0; k < mesh->LocalNz; k++) {
            phivalue += phi(mesh->xend, j, k);
          }
          phivalue /= mesh->LocalNz; // Average in Z of point next to boundary

          // Old value of phi at boundary
          BoutReal oldvalue = 0.5 * (phi(mesh->xend + 1, j, 0) + phi(mesh->xend, j, 0));

          // New value of phi at boundary, relaxing towards phivalue
          BoutReal newvalue = weight * oldvalue + (1. - weight) * phivalue;

          // Set phi at the boundary to this value
          for (int k = 0; k < mesh->LocalNz; k++) {
            phi(mesh->xend + 1, j, k) = 2.*newvalue - phi(mesh->xend, j, k);

            // Note: This seems to make a difference, but don't know why.
            // Without this, get convergence failures with no apparent instability
            // (all fields apparently smooth, well behaved)
            phi(mesh->xend + 2, j, k) = phi(mesh->xend + 1, j, k);
          }
        }
      }
    }
  } else {
    // phi_boundary_relax = false
    //
    // Set boundary from temperature, to be consistent with j=0 at sheath

    // Sheath multiplier Te -> phi (2.84522 for Deuterium)
    BoutReal sheathmult = 0.0;
    if (sheath_boundary) {
      BoutReal Me_Mp = get<BoutReal>(state["species"]["e"]["AA"]);
      sheathmult = log(0.5 * sqrt(1. / (Me_Mp * PI)));
    }

    Field3D Te; // Electron temperature, use for outer boundary conditions
    if (state["species"]["e"].isSet("temperature")) {
      // Electron temperature set
      Te = GET_NOBOUNDARY(Field3D, state["species"]["e"]["temperature"]);
    } else {
      Te = 0.0;
    }

    // Sheath multiplier Te -> phi (2.84522 for Deuterium if Ti = 0)
    if (mesh->firstX()) {
      for (int j = mesh->ystart; j <= mesh->yend; j++) {
        BoutReal teavg = 0.0; // Average Te in Z
        for (int k = 0; k < mesh->LocalNz; k++) {
          teavg += Te(mesh->xstart, j, k);
        }
        teavg /= mesh->LocalNz;
        BoutReal phivalue = sheathmult * teavg;
        // Set midpoint (boundary) value
        for (int k = 0; k < mesh->LocalNz; k++) {
          phi(mesh->xstart - 1, j, k) = 2. * phivalue - phi(mesh->xstart, j, k);

          // Note: This seems to make a difference, but don't know why.
          // Without this, get convergence failures with no apparent instability
          // (all fields apparently smooth, well behaved)
          phi(mesh->xstart - 2, j, k) = phi(mesh->xstart - 1, j, k);
        }
      }
    }

    if (mesh->lastX()) {
      for (int j = mesh->ystart; j <= mesh->yend; j++) {
        BoutReal teavg = 0.0; // Average Te in Z
        for (int k = 0; k < mesh->LocalNz; k++) {
          teavg += Te(mesh->xend, j, k);
        }
        teavg /= mesh->LocalNz;
        BoutReal phivalue = sheathmult * teavg;
        // Set midpoint (boundary) value
        for (int k = 0; k < mesh->LocalNz; k++) {
          phi(mesh->xend + 1, j, k) = 2. * phivalue - phi(mesh->xend, j, k);

          // Note: This seems to make a difference, but don't know why.
          // Without this, get convergence failures with no apparent instability
          // (all fields apparently smooth, well behaved)
          phi(mesh->xend + 2, j, k) = phi(mesh->xend + 1, j, k);
        }
      }
    }
  }

  // Update boundary conditions. Two issues:
  // 1) Solving here for phi + Pi, and then subtracting Pi from the result
  //    The boundary values should therefore include Pi
  // 2) The INVERT_SET flag takes the value in the guard (boundary) cell
  //    and sets the boundary between cells to this value.
  //    This shift by 1/2 grid cell is important.

  Field3D phi_plus_pi = phi + Pi_sum;

  if (mesh->firstX()) {
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      for (int k = 0; k < mesh->LocalNz; k++) {
        // Average phi + Pi at the boundary, and set the boundary cell
        // to this value. The phi solver will then put the value back
        // onto the cell mid-point
        phi_plus_pi(mesh->xstart - 1, j, k) =
          0.5
          * (phi_plus_pi(mesh->xstart - 1, j, k) +
             phi_plus_pi(mesh->xstart, j, k));
      }
    }
  }

  if (mesh->lastX()) {
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      for (int k = 0; k < mesh->LocalNz; k++) {
        phi_plus_pi(mesh->xend + 1, j, k) =
          0.5
          * (phi_plus_pi(mesh->xend + 1, j, k) +
             phi_plus_pi(mesh->xend, j, k));
      }
    }
  }

  // Calculate potential
  if (split_n0) {
    ////////////////////////////////////////////
    // Split into axisymmetric and non-axisymmetric components
    Field2D Vort2D = DC(Vort); // n=0 component
    Field2D phi_plus_pi_2d = DC(phi_plus_pi);
    phi_plus_pi -= phi_plus_pi_2d;

    phi_plus_pi_2d = laplacexy->solve(Vort2D, phi_plus_pi_2d);

    // Solve non-axisymmetric part using X-Z solver
    phi = phi_plus_pi_2d
      + phiSolver->solve((Vort - Vort2D) * (Bsq / average_atomic_mass),
                         phi_plus_pi)
      - Pi_sum;

  } else {
    phi = phiSolver->solve(Vort * (Bsq / average_atomic_mass),
                           phi_plus_pi)
      - Pi_sum;
  }

  // Ensure that potential is set in the communication guard cells
  mesh->communicate(phi);

  ddt(Vort) = 0.0;

  if (diamagnetic) {
    // Diamagnetic current. This is calculated here so that the energy sources/sinks
    // can be calculated for the evolving species.

    Vector3D Jdia; Jdia.x = 0.0; Jdia.y = 0.0; Jdia.z = 0.0;
    Jdia.covariant = Curlb_B.covariant;

    Options& allspecies = state["species"];
    
    for (auto& kv : allspecies.getChildren()) {
      Options& species = allspecies[kv.first]; // Note: need non-const

      if (!(IS_SET_NOBOUNDARY(species["pressure"]) and IS_SET(species["charge"]))) {
        continue; // No pressure or charge -> no diamagnetic current
      }
      // Note that the species must have a charge, but charge is not used,
      // because it cancels out in the expression for current
      
      auto P = GET_NOBOUNDARY(Field3D, species["pressure"]);

      Vector3D Jdia_species = P * Curlb_B; // Diamagnetic current for this species
      
      // This term energetically balances diamagnetic term
      // in the vorticity equation
      subtract(species["energy_source"],
               Jdia_species * Grad(phi));

      Jdia += Jdia_species; // Collect total diamagnetic current
    }
    
    // Note: This term is central differencing so that it balances
    // the corresponding compression term in the species pressure equations
    Field3D DivJdia = Div(Jdia);
    ddt(Vort) += DivJdia;

    if (diamagnetic_polarisation) {
      // Calculate energy exchange term nonlinear in pressure
      // ddt(Pi) += Pi * Div((Pe + Pi) * Curlb_B);
      for (auto& kv : allspecies.getChildren()) {
        Options& species = allspecies[kv.first]; // Note: need non-const

        if (!(IS_SET_NOBOUNDARY(species["pressure"]) and species.isSet("charge") and species.isSet("AA"))) {
          continue; // No pressure, charge or mass -> no polarisation current due to diamagnetic flow
        }
        auto P = GET_NOBOUNDARY(Field3D, species["pressure"]);
        auto AA = get<BoutReal>(species["AA"]);
        
        add(species["energy_source"],
            (3./2) * P * (AA / average_atomic_mass) * DivJdia);
      }
    }

    set(fields["DivJdia"], DivJdia);
  }

  set(fields["vorticity"], Vort);
  set(fields["phi"], phi);
}
  
void Vorticity::finally(const Options &state) {
  AUTO_TRACE();

  phi = get<Field3D>(state["fields"]["phi"]);

  if (exb_advection) {
    ddt(Vort) -= Div_n_bxGrad_f_B_XPPM(Vort, phi, bndry_flux, poloidal_flows);

    /*
      ddt(Vort) -=
        Div_n_bxGrad_f_B_XPPM(0.5 * Vort, phi, bndry_flux, poloidal_flows);

    // V_ExB dot Grad(Pi)
    Field3D vEdotGradPi = bracket(phi, Pi, BRACKET_ARAKAWA);
    vEdotGradPi.applyBoundary("free_o2");
    // delp2(phi) term
    Field3D DelpPhi_2B2 = 0.5 * Delp2(phi) / Bsq;
    DelpPhi_2B2.applyBoundary("free_o2");
    
    mesh->communicate(vEdotGradPi, DelpPhi_2B2);
    
    ddt(Vort) -= FV::Div_a_Laplace_perp(0.5 / Bsq, vEdotGradPi);
    */
  }
  
  if (state.isSection("fields") and state["fields"].isSet("DivJextra")) {
    auto DivJextra = get<Field3D>(state["fields"]["DivJextra"]);

    // Parallel current is handled here, to allow different 2D or 3D closures
    // to be used
    ddt(Vort) += DivJextra;
  }

  // Parallel current due to species parallel flow
  for (auto& kv : state["species"].getChildren()) {
    const Options& species = kv.second;

    if (!species.isSet("charge") or !species.isSet("momentum")) {
      continue; // Not charged, or no parallel flow
    }
    const BoutReal Z = get<BoutReal>(species["charge"]);
    if (fabs(Z) < 1e-5) {
      continue; // Not charged
    }

    const Field3D N = get<Field3D>(species["density"]);
    const Field3D NV = get<Field3D>(species["momentum"]);
    const BoutReal A = get<BoutReal>(species["AA"]);

    // Note: Using NV rather than N*V so that the cell boundary flux is correct
    ddt(Vort) += Div_par((Z / A) * NV);
  }

  if (collisional_friction) {
    // Damping of vorticity due to collisions

    // Calculate a mass-weighted collision frequency
    Field3D sum_A_nu_n = zeroFrom(Vort); // Sum of atomic mass * collision frequency * density
    Field3D sum_A_n = zeroFrom(Vort);    // Sum of atomic mass * density

    const Options& allspecies = state["species"];
    for (const auto& kv : allspecies.getChildren()) {
      const Options& species = kv.second;

      if (!(species.isSet("charge") and species.isSet("AA"))) {
        continue; // No charge or mass -> no current
      }
      if (fabs(get<BoutReal>(species["charge"])) < 1e-5) {
	continue; // Zero charge
      }

      const BoutReal A = get<BoutReal>(species["AA"]);
      const Field3D N = get<Field3D>(species["density"]);
      const Field3D AN = A * N;
      sum_A_n += AN;
      if (species.isSet("collision_frequency")) {
	sum_A_nu_n += AN * get<Field3D>(species["collision_frequency"]);
      }
    }

    const Field3D weighted_collision_frequency = sum_A_nu_n / sum_A_n;

    ddt(Vort) -= FV::Div_a_Laplace_perp(
		    weighted_collision_frequency * average_atomic_mass
		    / Bsq, phi);
  }

  if (state.isSet("sound_speed")) {
    Field3D sound_speed = get<Field3D>(state["sound_speed"]);
    if (vort_dissipation) {
      // Adds dissipation term like in other equations
      ddt(Vort) -= FV::Div_par(Vort, 0.0, sound_speed);
    }

    if (phi_dissipation) {
      // Adds dissipation term like in other equations, but depending on gradient of potential
      ddt(Vort) -= FV::Div_par(-phi, 0.0, sound_speed);
    }
  }

  if (hyper_z > 0) {
    // Form of hyper-viscosity to suppress zig-zags in Z
    auto* coord = Vort.getCoordinates();
    ddt(Vort) -= hyper_z * SQ(SQ(coord->dz)) * D4DZ4(Vort);
  }
}


