
#include "../include/vorticity.hxx"
#include "../include/div_ops.hxx"

#include <bout/constants.hxx>
#include <bout/fv_ops.hxx>
#include <bout/invert/laplacexy.hxx>
#include <bout/derivs.hxx>
#include <bout/difops.hxx>
#include <bout/invert_laplace.hxx>

using bout::globals::mesh;

namespace {
BoutReal floor(BoutReal value, BoutReal min) {
  if (value < min)
    return min;
  return value;
}

Ind3D indexAt(const Field3D& f, int x, int y, int z) {
  int ny = f.getNy();
  int nz = f.getNz();
  return Ind3D{(x * ny + y) * nz + z, ny, nz};
}

/// Limited free gradient of log of a quantity
/// This ensures that the guard cell values remain positive
/// while also ensuring that the quantity never increases
///
///  fm  fc | fp
///         ^ boundary
///
/// exp( 2*log(fc) - log(fm) )
///
BoutReal limitFree(BoutReal fm, BoutReal fc) {
  if (fm < fc) {
    return fc; // Neumann rather than increasing into boundary
  }
  if (fm < 1e-10) {
    return fc; // Low / no density condition
  }
  BoutReal fp = SQ(fc) / fm;
#if CHECKLEVEL >= 2
  if (!std::isfinite(fp)) {
    throw BoutException("SheathBoundary limitFree: {}, {} -> {}", fm, fc, fp);
  }
#endif

  return fp;
}

void applyDirichletBoundary(Field2D &f) {
  for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
    f(r.ind, mesh->ystart - 1) = -f(r.ind, mesh->ystart);
  }
  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
    f(r.ind, mesh->yend + 1) = -f(r.ind, mesh->yend);
  }

  if (mesh->lastX()) {
    for (int y = mesh->ystart; y <= mesh->yend; y++) {
      f(mesh->xend + 1, y) = -f(mesh->xend, y);
    }
  }
  // Not applying to core boundary
  if (mesh->firstX() && !mesh->periodicY(mesh->xstart)) {
    for (int y = mesh->ystart; y <= mesh->yend; y++) {
      f(mesh->xstart - 1, y) = -f(mesh->xstart, y);
    }
  }
}
}

Vorticity::Vorticity(std::string name, Options& alloptions, Solver* solver) {
  AUTO_TRACE();

  solver->add(Vort, "Vort");

  auto& options = alloptions[name];
  // Normalisations
  const Options& units = alloptions["units"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();
  const BoutReal Bnorm = units["Tesla"];
  const BoutReal Lnorm = units["meters"];

  exb_advection = options["exb_advection"]
                      .doc("Include ExB advection (nonlinear term)?")
                      .withDefault<bool>(true);

  exb_advection_simplified = options["exb_advection_simplified"]
                      .doc("Simplify nonlinear ExB advection form?")
                      .withDefault<bool>(true);

  diamagnetic =
      options["diamagnetic"].doc("Include diamagnetic current?").withDefault<bool>(true);

  sheath_boundary = options["sheath_boundary"]
                        .doc("Set potential to j=0 sheath at radial boundaries? (default = 0)")
                        .withDefault<bool>(false);

  diamagnetic_polarisation =
      options["diamagnetic_polarisation"]
          .doc("Include diamagnetic drift in polarisation current?")
          .withDefault<bool>(true);

  collisional_friction =
      options["collisional_friction"]
          .doc("Damp vorticity based on mass-weighted collision frequency")
          .withDefault<bool>(false);

  average_atomic_mass = options["average_atomic_mass"]
                            .doc("Weighted average atomic mass, for polarisation current "
                                 "(Boussinesq approximation)")
                            .withDefault<BoutReal>(2.0); // Deuterium

  bndry_flux = options["bndry_flux"]
                   .doc("Allow flows through radial boundaries")
                   .withDefault<bool>(true);

  poloidal_flows =
      options["poloidal_flows"].doc("Include poloidal ExB flow").withDefault<bool>(true);

  split_n0 = options["split_n0"]
                 .doc("Split phi into n=0 and n!=0 components")
                 .withDefault<bool>(false);

  viscosity = options["viscosity"]
    .doc("Kinematic viscosity [m^2/s]")
    .withDefault<BoutReal>(0.0)
    / (Lnorm * Lnorm * Omega_ci);
  viscosity.applyBoundary("dirichlet");

  hyper_z = options["hyper_z"].doc("Hyper-viscosity in Z. < 0 -> off").withDefault(-1.0);

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

  phi_sheath_dissipation = options["phi_sheath_dissipation"]
    .doc("Add dissipation when phi < 0.0 at the sheath")
    .withDefault<bool>(false);

  damp_core_vorticity = options["damp_core_vorticity"]
	  .doc("Damp vorticity at the core boundary?")
	  .withDefault<bool>(false);

  // Add phi to restart files so that the value in the boundaries
  // is restored on restart. This is done even when phi is not evolving,
  // so that phi can be saved and re-loaded

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

  } catch (BoutException& e) {
    try {
      // May be 2D, reading as 3D
      Vector2D curv2d;
      curv2d.covariant = false;
      mesh->get(curv2d, "bxcv");
      Curlb_B = curv2d;
    } catch (BoutException& e) {
      if (diamagnetic) {
        // Need curvature
        throw;
      } else {
        output_warn.write("No curvature vector in input grid");
        Curlb_B = 0.0;
      }
    }
  }

  if (Options::root()["mesh"]["paralleltransform"]["type"].as<std::string>()
      == "shifted") {
    Field2D I;
    mesh->get(I, "sinty");
    Curlb_B.z += I * Curlb_B.x;
  }

  Curlb_B.x /= Bnorm;
  Curlb_B.y *= SQ(Lnorm);
  Curlb_B.z *= SQ(Lnorm);

  Curlb_B *= 2. / coord->Bxy;

  applyDirichletBoundary(Curlb_B.x);
  applyDirichletBoundary(Curlb_B.y);
  applyDirichletBoundary(Curlb_B.z);

  Bsq = SQ(coord->Bxy);

  diagnose = options["diagnose"]
    .doc("Output additional diagnostics?")
    .withDefault<bool>(false);
}

void Vorticity::transform(Options& state) {
  AUTO_TRACE();

  auto& fields = state["fields"];

  // Set the boundary of phi. Both 2D and 3D fields are kept, though the 3D field
  // is constant in Z. This is for efficiency, to reduce the number of conversions.
  // Note: For now the boundary values are all at the midpoint,
  //       and only phi is considered, not phi + Pi which is handled in Boussinesq solves
  Pi_hat = 0.0; // Contribution from ion pressure, weighted by atomic mass / charge
  if (diamagnetic_polarisation) {
    // Diamagnetic term in vorticity. Note this is weighted by the mass
    // This includes all species, including electrons
    Options& allspecies = state["species"];
    for (auto& kv : allspecies.getChildren()) {
      Options& species = allspecies[kv.first]; // Note: need non-const

      if (!(IS_SET_NOBOUNDARY(species["pressure"]) and species.isSet("charge")
            and species.isSet("AA"))) {
        continue; // No pressure, charge or mass -> no polarisation current
      }

      const auto charge = get<BoutReal>(species["charge"]);
      if (fabs(charge) < 1e-5) {
        // No charge
        continue;
      }

      // Don't need sheath boundary
      const auto P = GET_NOBOUNDARY(Field3D, species["pressure"]);
      const auto AA = get<BoutReal>(species["AA"]);

      Pi_hat += P * (AA / average_atomic_mass / charge);
    }
  }

  Pi_hat.applyBoundary("neumann");

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
          BoutReal newvalue = weight * oldvalue + (1. - weight) * phivalue;

          // Set phi at the boundary to this value
          for (int k = 0; k < mesh->LocalNz; k++) {
            phi(mesh->xstart - 1, j, k) = 2. * newvalue - phi(mesh->xstart, j, k);

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
            phi(mesh->xend + 1, j, k) = 2. * newvalue - phi(mesh->xend, j, k);

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

  Field3D phi_plus_pi = phi + Pi_hat;

  if (mesh->firstX()) {
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      for (int k = 0; k < mesh->LocalNz; k++) {
        // Average phi + Pi at the boundary, and set the boundary cell
        // to this value. The phi solver will then put the value back
        // onto the cell mid-point
        phi_plus_pi(mesh->xstart - 1, j, k) =
            0.5 * (phi_plus_pi(mesh->xstart - 1, j, k) + phi_plus_pi(mesh->xstart, j, k));
      }
    }
  }

  if (mesh->lastX()) {
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      for (int k = 0; k < mesh->LocalNz; k++) {
        phi_plus_pi(mesh->xend + 1, j, k) =
            0.5 * (phi_plus_pi(mesh->xend + 1, j, k) + phi_plus_pi(mesh->xend, j, k));
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
          + phiSolver->solve((Vort - Vort2D) * (Bsq / average_atomic_mass), phi_plus_pi)
          - Pi_hat;

  } else {
    phi = phiSolver->solve(Vort * (Bsq / average_atomic_mass), phi_plus_pi) - Pi_hat;
  }

  // Ensure that potential is set in the communication guard cells
  mesh->communicate(phi);

  // Outer boundary cells
  if (mesh->firstX()) {
    for (int i = mesh->xstart - 2; i >= 0; --i) {
      for (int j = mesh->ystart; j <= mesh->yend; ++j) {
        for (int k = 0; k < mesh->LocalNz; ++k) {
          phi(i, j, k) = phi(i + 1, j, k);
        }
      }
    }
  }
  if (mesh->lastX()) {
    for (int i = mesh->xend + 2; i < mesh->LocalNx; ++i) {
      for (int j = mesh->ystart; j <= mesh->yend; ++j) {
        for (int k = 0; k < mesh->LocalNz; ++k) {
          phi(i, j, k) = phi(i - 1, j, k);
        }
      }
    }
  }

  ddt(Vort) = 0.0;

  if (diamagnetic) {
    // Diamagnetic current. This is calculated here so that the energy sources/sinks
    // can be calculated for the evolving species.

    Vector3D Jdia;
    Jdia.x = 0.0;
    Jdia.y = 0.0;
    Jdia.z = 0.0;
    Jdia.covariant = Curlb_B.covariant;

    Options& allspecies = state["species"];

    for (auto& kv : allspecies.getChildren()) {
      Options& species = allspecies[kv.first]; // Note: need non-const

      if (!(IS_SET_NOBOUNDARY(species["pressure"]) and IS_SET(species["charge"]))) {
        continue; // No pressure or charge -> no diamagnetic current
      }
      if (fabs(get<BoutReal>(species["charge"])) < 1e-5) {
        // No charge
        continue;
      }

      // Note that the species must have a charge, but charge is not used,
      // because it cancels out in the expression for current

      auto P = GET_NOBOUNDARY(Field3D, species["pressure"]);

      P.clearParallelSlices();
      phi.clearParallelSlices();
      
      // Note: We need boundary conditions on P, so apply the same
      //       free boundary condition as sheath_boundary.
      if (P.hasParallelSlices()) {
        Field3D &P_ydown = P.ydown();
        Field3D &P_yup = P.yup();
        for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            P_ydown(r.ind, mesh->ystart - 1, jz) = P(r.ind, mesh->ystart, jz); //2 * P(r.ind, mesh->ystart, jz) - P_yup(r.ind, mesh->ystart + 1, jz);
          }
        }
        for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            P_yup(r.ind, mesh->yend + 1, jz) = P(r.ind, mesh->yend, jz); //2 * P(r.ind, mesh->yend, jz) - P_ydown(r.ind, mesh->yend - 1, jz);
          }
        }
      } else {
        Field3D P_fa = toFieldAligned(P);
        for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            auto i = indexAt(P_fa, r.ind, mesh->ystart, jz);
            P_fa[i.ym()] = P_fa[i]; //limitFree(P_fa[i.yp()], P_fa[i]);
          }
        }
        for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            auto i = indexAt(P_fa, r.ind, mesh->yend, jz);
            P_fa[i.yp()] = P_fa[i]; //limitFree(P_fa[i.ym()], P_fa[i]);
          }
        }
        P = fromFieldAligned(P_fa);
      }

      // Note: This calculation requires phi derivatives at the Y boundaries
      //       Setting to free boundaries
      if (phi.hasParallelSlices()) {
        Field3D &phi_ydown = phi.ydown();
        Field3D &phi_yup = phi.yup();
        for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            phi_ydown(r.ind, mesh->ystart - 1, jz) = 2 * phi(r.ind, mesh->ystart, jz) - phi_yup(r.ind, mesh->ystart + 1, jz);
          }
        }
        for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            phi_yup(r.ind, mesh->yend + 1, jz) = 2 * phi(r.ind, mesh->yend, jz) - phi_ydown(r.ind, mesh->yend - 1, jz);
          }
        }
      } else {
        Field3D phi_fa = toFieldAligned(phi);
        for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            phi_fa(r.ind, mesh->ystart - 1, jz) = 2 * phi_fa(r.ind, mesh->ystart, jz) - phi_fa(r.ind, mesh->ystart + 1, jz);
          }
        }
        for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            phi_fa(r.ind, mesh->yend + 1, jz) = 2 * phi_fa(r.ind, mesh->yend, jz) - phi_fa(r.ind, mesh->yend - 1, jz);
          }
        }
        phi = fromFieldAligned(phi_fa);
      }

      Vector3D Jdia_species = P * Curlb_B; // Diamagnetic current for this species

      // This term energetically balances diamagnetic term
      // in the vorticity equation
      subtract(species["energy_source"], Jdia_species * Grad(phi));

      Jdia += Jdia_species; // Collect total diamagnetic current
    }

    // Note: This term is central differencing so that it balances
    // the corresponding compression term in the species pressure equations
    DivJdia = Div(Jdia);
    ddt(Vort) += DivJdia;

    set(fields["DivJdia"], DivJdia);
  }

  if (collisional_friction) {
    // Damping of vorticity due to collisions

    // Calculate a mass-weighted collision frequency
    Field3D sum_A_nu_n =
        zeroFrom(Vort); // Sum of atomic mass * collision frequency * density
    Field3D sum_A_n = zeroFrom(Vort); // Sum of atomic mass * density

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
      const Field3D N = GET_NOBOUNDARY(Field3D, species["density"]);
      const Field3D AN = A * N;
      sum_A_n += AN;
      if (IS_SET(species["collision_frequency"])) {
        sum_A_nu_n += AN * GET_VALUE(Field3D, species["collision_frequency"]);
      }
    }

    Field3D weighted_collision_frequency = sum_A_nu_n / sum_A_n;
    weighted_collision_frequency.applyBoundary("neumann");

    DivJcol = -FV::Div_a_Grad_perp(
        weighted_collision_frequency * average_atomic_mass / Bsq, phi + Pi_hat);

    ddt(Vort) += DivJcol;
    set(fields["DivJcol"], DivJcol);
  }

  set(fields["vorticity"], Vort);
  set(fields["phi"], phi);
}

void Vorticity::finally(const Options& state) {
  AUTO_TRACE();

  phi = get<Field3D>(state["fields"]["phi"]);

  if (exb_advection) {
    // These terms come from divergence of polarisation current

    if (exb_advection_simplified) {
      // By default this is a simplified nonlinear term
      ddt(Vort) -= Div_n_bxGrad_f_B_XPPM(Vort, phi, bndry_flux, poloidal_flows);

    } else {
      // If diamagnetic_polarisation = false and B is constant, then
      // this term reduces to the simplified form above.
      //
      // Because this is implemented in terms of an operation on the result
      // of an operation, we need to communicate and the resulting stencil is
      // wider than the simple form.
      ddt(Vort) -=
        Div_n_bxGrad_f_B_XPPM(0.5 * Vort, phi, bndry_flux, poloidal_flows);

      // V_ExB dot Grad(Pi)
      Field3D vEdotGradPi = bracket(phi, Pi_hat, BRACKET_ARAKAWA);
      vEdotGradPi.applyBoundary("free_o2");

      // delp2(phi) term
      Field3D DelpPhi_2B2 = 0.5 * average_atomic_mass * Delp2(phi) / Bsq;
      DelpPhi_2B2.applyBoundary("free_o2");

      mesh->communicate(vEdotGradPi, DelpPhi_2B2);

      ddt(Vort) -= FV::Div_a_Grad_perp(0.5 * average_atomic_mass / Bsq, vEdotGradPi);
      ddt(Vort) -= Div_n_bxGrad_f_B_XPPM(DelpPhi_2B2, phi + Pi_hat, bndry_flux,
                                         poloidal_flows);
    }
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

  // Viscosity
  ddt(Vort) += FV::Div_a_Grad_perp(viscosity, Vort);

  if (vort_dissipation) {
    // Adds dissipation term like in other equations
    Field3D sound_speed = get<Field3D>(state["sound_speed"]);
    ddt(Vort) -= FV::Div_par(Vort, 0.0, sound_speed);
  }

  if (phi_dissipation) {
    // Adds dissipation term like in other equations, but depending on gradient of
    // potential
    Field3D sound_speed = get<Field3D>(state["sound_speed"]);
    ddt(Vort) -= FV::Div_par(-phi, 0.0, sound_speed);
  }

  if (hyper_z > 0) {
    // Form of hyper-viscosity to suppress zig-zags in Z
    auto* coord = Vort.getCoordinates();
    ddt(Vort) -= hyper_z * SQ(SQ(coord->dz)) * D4DZ4(Vort);
  }

  if (phi_sheath_dissipation) {
    // Dissipation when phi < 0.0 at the sheath

    auto phi_fa = toFieldAligned(phi);
    Field3D dissipation{zeroFrom(phi_fa)};
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        auto i = indexAt(phi_fa, r.ind, mesh->ystart, jz);
        BoutReal phisheath = 0.5*(phi_fa[i] + phi_fa[i.ym()]);
        dissipation[i] = -floor(-phisheath, 0.0);
      }
    }

    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        auto i = indexAt(phi_fa, r.ind, mesh->yend, jz);
        BoutReal phisheath = 0.5*(phi_fa[i] + phi_fa[i.yp()]);
        dissipation[i] = -floor(-phisheath, 0.0);
      }
    }
    ddt(Vort) += fromFieldAligned(dissipation);
  }

  if (damp_core_vorticity) {
    // Damp axisymmetric vorticity near core boundary
    if (mesh->firstX() and mesh->periodicY(mesh->xstart)) {
      for (int j = mesh->ystart; j <= mesh->yend; j++) {
        BoutReal vort_avg = 0.0; // Average Vort in Z
        for (int k = 0; k < mesh->LocalNz; k++) {
          vort_avg += Vort(mesh->xstart, j, k);
        }
        vort_avg /= mesh->LocalNz;
        for (int k = 0; k < mesh->LocalNz; k++) {
          ddt(Vort)(mesh->xstart, j, k) -= 0.01 * vort_avg;
        }
      }
    }
  }
}

void Vorticity::outputVars(Options& state) {
  AUTO_TRACE();
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);

  state["Vort"].setAttributes({{"time_dimension", "t"},
                               {"units", "C m^-3"},
                               {"conversion", SI::qe * Nnorm},
                               {"long_name", "vorticity"},
                               {"source", "vorticity"}});

  set_with_attrs(state["phi"], phi,
                 {{"time_dimension", "t"},
                  {"units", "V"},
                  {"conversion", Tnorm},
                  {"standard_name", "potential"},
                  {"long_name", "plasma potential"},
                  {"source", "vorticity"}});

  if (diagnose) {
    set_with_attrs(state["ddt(Vort)"], ddt(Vort),
                   {{"time_dimension", "t"},
                    {"units", "A m^-3"},
                    {"conversion", SI::qe * Nnorm * Omega_ci},
                    {"long_name", "Rate of change of vorticity"},
                    {"source", "vorticity"}});

    if (diamagnetic) {
      set_with_attrs(state["DivJdia"], DivJdia,
                     {{"time_dimension", "t"},
                      {"units", "A m^-3"},
                      {"conversion", SI::qe * Nnorm * Omega_ci},
                      {"long_name", "Divergence of diamagnetic current"},
                      {"source", "vorticity"}});
    }
    if (collisional_friction) {
      set_with_attrs(state["DivJcol"], DivJcol,
                     {{"time_dimension", "t"},
                      {"units", "A m^-3"},
                      {"conversion", SI::qe * Nnorm * Omega_ci},
                      {"long_name", "Divergence of collisional current"},
                      {"source", "vorticity"}});
    }
  }
}
