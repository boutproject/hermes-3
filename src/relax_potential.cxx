#include <bout/fv_ops.hxx>
#include <bout/solver.hxx>

using bout::globals::mesh;

#include "../include/div_ops.hxx"
#include "../include/relax_potential.hxx"
#include <bout/constants.hxx>
#include "../include/hermes_build_config.hxx"


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
}



RelaxPotential::RelaxPotential(std::string name, Options& alloptions, Solver* solver) {
  AUTO_TRACE();


  solver->add(Vort, "Vort"); // Vorticity evolving
  solver->add(phi1, "phi1"); // Evolving scaled potential ϕ_1 = λ_2 ϕ

  const Options& units = alloptions["units"];
  // Normalisations
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();
  const BoutReal Bnorm = units["Tesla"];
  const BoutReal Lnorm = units["meters"];
  const BoutReal Tnorm = units["eV"];
  const BoutReal Nnorm = units["inv_meters_cubed"];

  auto& options = alloptions[name];

  exb_advection = options["exb_advection"]
                      .doc("Include nonlinear ExB advection?")
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

  boussinesq = options["boussinesq"]
                   .doc("Use the Boussinesq approximation?")
                   .withDefault<bool>(true);


  average_atomic_mass = options["average_atomic_mass"]
                            .doc("Weighted average atomic mass, for polarisaion current "
                                 "(Boussinesq approximation)")
                            .withDefault<BoutReal>(2.0); // Deuterium

  bndry_flux = options["bndry_flux"]
                   .doc("Allow flows through radial boundaries")
                   .withDefault<bool>(true);

  poloidal_flows =
      options["poloidal_flows"].doc("Include poloidal ExB flow").withDefault<bool>(true);

  viscosity = options["viscosity"]
    .doc("Kinematic viscosity [m^2/s]")
    .withDefault<Field3D>(0.0)
    / (Lnorm * Lnorm * Omega_ci);
  mesh->communicate(viscosity);
  viscosity.splitParallelSlices(); // We need this otherwise applyParallelBoundary gives and assertion error. 
  viscosity.applyBoundary("dirichlet");
  viscosity.applyParallelBoundary("parallel_dirichlet_o2");

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



  lambda_1 = options["lambda_1"].doc("λ_1 > 1").withDefault(100);
  lambda_2 = options["lambda_2"].doc("λ_2 > λ_1").withDefault(1e5);

  // NOTE(malamast): How do we do that?
  // Add phi to restart files so that the value in the boundaries
  // is restored on restart. This is done even when phi is not evolving,
  // so that phi can be saved and re-loaded

  // Set initial value. Will be overwritten if restarting
  phi1 = 0.0;
  Vort = 0.0;

  if (phi_boundary_relax) {
    // Set the last update time to -1, so it will reset
    // the first time RHS function is called
    phi_boundary_last_update = -1.;

    phi_core_averagey = options["phi_core_averagey"]
      .doc("Average phi core boundary in Y?")
      .withDefault<bool>(false) and mesh->periodicY(mesh->xstart);

    phi_boundary_timescale = options["phi_boundary_timescale"]
                                 .doc("Timescale for phi boundary relaxation [seconds]")
                                 .withDefault(1e-4)
                             / get<BoutReal>(alloptions["units"]["seconds"]);
  }


  auto coord = mesh->getCoordinates();

  if (diamagnetic) {
    // Read curvature vector
    try {
      Curlb_B.covariant = false; // Contravariant
      mesh->get(Curlb_B, "bxcv");

    } catch (BoutException& e) {
      // May be 2D, reading as 3D
      Vector2D curv2d;
      curv2d.covariant = false;
      if (mesh->get(curv2d, "bxcv")) {
        throw BoutException("Curvature vector not found in input");
      }
      Curlb_B = curv2d;
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
  }

  Bsq = SQ(coord->Bxy);

  diagnose = options["diagnose"]
    .doc("Output additional diagnostics?")
    .withDefault<bool>(false);

}

void RelaxPotential::transform(Options& state) {
  AUTO_TRACE();

  phi.name = "phi";
  auto& fields = state["fields"];

  // Scale potential
  phi1.applyBoundary("neumann");

  phi = phi1 / lambda_2;
  Vort.applyBoundary("neumann");

  // phi.applyBoundary("neumann");
  // Set the boundary of phi. 
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
        BoutReal phivalue = 0.0;
        if (phi_core_averagey) {
          BoutReal philocal = 0.0;
          for (int j = mesh->ystart; j <= mesh->yend; j++) {
            for (int k = 0; k < mesh->LocalNz; k++) {
              philocal += phi(mesh->xstart, j, k);
            }
          }
          MPI_Comm comm_inner = mesh->getYcomm(0);
          int np;
          MPI_Comm_size(comm_inner, &np);
          MPI_Allreduce(&philocal,
                        &phivalue,
                        1, MPI_DOUBLE,
                        MPI_SUM, comm_inner);
          phivalue /= (np * mesh->LocalNz * mesh->LocalNy);
        }

        for (int j = mesh->ystart; j <= mesh->yend; j++) {
          if (!phi_core_averagey) {
            phivalue = 0.0; // Calculate phi boundary for each Y index separately
            for (int k = 0; k < mesh->LocalNz; k++) {
              phivalue += phi(mesh->xstart, j, k);
            }
            phivalue /= mesh->LocalNz; // Average in Z of point next to boundary
          }

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

  // if (mesh->firstX()) {
  //   for (int j = mesh->ystart; j <= mesh->yend; j++) {
  //     for (int k = 0; k < mesh->LocalNz; k++) {
  //       // Average phi + Pi at the boundary, and set the boundary cell
  //       // to this value. The phi solver will then put the value back
  //       // onto the cell mid-point
  //       phi(mesh->xstart - 1, j, k) =
  //           0.5 * (phi(mesh->xstart - 1, j, k) + phi(mesh->xstart, j, k));
  //     }
  //   }
  // }

  // if (mesh->lastX()) {
  //   for (int j = mesh->ystart; j <= mesh->yend; j++) {
  //     for (int k = 0; k < mesh->LocalNz; k++) {
  //       phi(mesh->xend + 1, j, k) =
  //           0.5 * (phi(mesh->xend + 1, j, k) + phi(mesh->xend, j, k));
  //     }
  //   }
  // }


  // Ensure that potential is set in the communication guard cells
  mesh->communicate(Vort, phi); //NOTE(malamast): Should we include phi1?

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

    // Pre-calculate this rather than calculate for each species
    Vector3D Grad_phi = Grad(phi);

    for (auto& kv : allspecies.getChildren()) {
      Options& species = allspecies[kv.first]; // Note: need non-const

      if (!(IS_SET_NOBOUNDARY(species["pressure"]) and IS_SET(species["charge"])
            and (get<BoutReal>(species["charge"]) != 0.0))) {
        continue; // No pressure or charge -> no diamagnetic current
      }
      // Note that the species must have a charge, but charge is not used,
      // because it cancels out in the expression for current

      auto P = GET_NOBOUNDARY(Field3D, species["pressure"]);

      Vector3D Jdia_species = P * Curlb_B; // Diamagnetic current for this species

      // This term energetically balances diamagnetic term
      // in the vorticity equation
      subtract(species["energy_source"], Jdia_species * Grad_phi);

      Jdia += Jdia_species; // Collect total diamagnetic current
    }

    // Note: This term is central differencing so that it balances
    // the corresponding compression term in the species pressure equations
    // Field3D DivJdia = Div(Jdia);
    DivJdia = Div(Jdia);
    ddt(Vort) += DivJdia;

    if (diamagnetic_polarisation) {
      // Calculate energy exchange term nonlinear in pressure
      // ddt(Pi) += Pi * Div((Pe + Pi) * Curlb_B);
      for (auto& kv : allspecies.getChildren()) {
        Options& species = allspecies[kv.first]; // Note: need non-const

        if (!(IS_SET_NOBOUNDARY(species["pressure"]) and IS_SET(species["charge"])
              and IS_SET(species["AA"]))) {
          continue; // No pressure, charge or mass -> no polarisation current due to
                    // rate of change of diamagnetic flow
        }
        auto P = GET_NOBOUNDARY(Field3D, species["pressure"]);

        add(species["energy_source"], (3. / 2) * P * DivJdia);
      }
    }

    set(fields["DivJdia"], DivJdia);
  }

  set(fields["vorticity"], Vort);
  set(fields["phi"], phi);
}

void RelaxPotential::finally(const Options& state) {
  AUTO_TRACE();

  const Options& allspecies = state["species"];

  phi = get<Field3D>(state["fields"]["phi"]);
  Vort = get<Field3D>(state["fields"]["vorticity"]);

  if (exb_advection) {
    ddt(Vort) -= Div_n_bxGrad_f_B_XPPM(Vort, phi, bndry_flux, poloidal_flows);
  }

  if (state.isSection("fields") and state["fields"].isSet("DivJextra")) {
    auto DivJextra = get<Field3D>(state["fields"]["DivJextra"]);

    // Parallel current is handled here, to allow different 2D or 3D closures
    // to be used
    ddt(Vort) += DivJextra;
  }

  // Parallel current due to species parallel flow
  for (auto& kv : allspecies.getChildren()) {
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

    Field3D zero {0.0};
    zero.splitParallelSlices();
    zero.yup() = 0.0;
    zero.ydown() = 0.0;

    Field3D dummy;
    ddt(Vort) -= FV::Div_par_mod<hermes::Limiter>(-phi, zero, sound_speed, dummy);
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

  // Solve diffusion equation for potential

  if (boussinesq) {
    ddt(phi1) =
        lambda_1 * (FV::Div_a_Grad_perp(average_atomic_mass / Bsq, phi) - Vort);

    if (diamagnetic_polarisation) {
      for (auto& kv : allspecies.getChildren()) {
        // Note: includes electrons (should it?)

        const Options& species = kv.second;
        if (!species.isSet("charge")) {
          continue; // Not charged
        }
        const BoutReal Z = get<BoutReal>(species["charge"]);
        if (fabs(Z) < 1e-5) {
          continue; // Not charged
        }
        if (!species.isSet("pressure")) {
          continue; // No pressure
        }
        const BoutReal A = get<BoutReal>(species["AA"]);
        const Field3D P = get<Field3D>(species["pressure"]);
        ddt(phi1) += lambda_1 * FV::Div_a_Grad_perp(A / Bsq, P);
      }
    }
  } else {
    // Non-Boussinesq. Calculate mass density by summing over species

    // Calculate vorticity from potential phi
    Field3D phi_vort = 0.0;
    for (auto& kv : allspecies.getChildren()) {
      const Options& species = kv.second;

      if (!species.isSet("charge")) {
        continue; // Not charged
      }
      const BoutReal Zi = get<BoutReal>(species["charge"]);
      if (fabs(Zi) < 1e-5) {
        continue; // Not charged
      }

      const BoutReal Ai = get<BoutReal>(species["AA"]);
      const Field3D Ni = get<Field3D>(species["density"]);

      phi_vort += FV::Div_a_Grad_perp((Ai / Bsq) * Ni, phi);

      if (diamagnetic_polarisation and species.isSet("pressure")) {
        // Calculate the diamagnetic flow contribution
        const Field3D Pi = get<Field3D>(species["pressure"]);
        phi_vort += FV::Div_a_Grad_perp(Ai / Bsq, Pi);
      }
    }

    ddt(phi1) = lambda_1 * (phi_vort - Vort);
  }
}

void RelaxPotential::outputVars(Options& state) {
  AUTO_TRACE();
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);

  // NOTE(Malamast): I added the below from the vorticity component. I need to check the units for Vort in the present component.
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
                  {"source", "relax_potential"}});

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
  }

}
