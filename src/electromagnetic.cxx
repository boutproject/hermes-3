#include "../include/electromagnetic.hxx"

#include <bout/constants.hxx>
#include <bout/mesh.hxx>
#include <bout/invert_laplace.hxx>

Electromagnetic::Electromagnetic(std::string name, Options &alloptions, Solver*) {
  AUTO_TRACE();

  Options& units = alloptions["units"];
  BoutReal Bnorm = units["Tesla"];
  BoutReal Tnorm = units["eV"];
  BoutReal Nnorm = units["inv_meters_cubed"];

  // Note: This is NOT plasma beta as usually defined, but
  // is a factor of 2 too small (if Ti=0) or 4 (if Ti = Te)
  // This is the normalisation parameter in the Helmholtz equation
  beta_em = SI::mu0 * SI::qe * Tnorm * Nnorm / SQ(Bnorm);
  output_info.write("Electromagnetic: beta_em = {}", beta_em);

  auto& options = alloptions[name];

  // Use the "Naulin" solver because we need to include toroidal
  // variations of the density (A coefficient)
  if (!options["laplacian"].isSet("type")) {
    options["laplacian"]["type"] = "naulin";
  }
  aparSolver = Laplacian::create(&options["laplacian"]);

  const_gradient = options["const_gradient"]
    .doc("Extrapolate gradient of Apar into all radial boundaries?")
    .withDefault<bool>(false);

  // Give Apar an initial value because we solve Apar by iteration
  // starting from the previous solution
  Apar = 0.0;

  if (const_gradient) {
    // Set flags to take the gradient from the RHS
    aparSolver->setInnerBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD + INVERT_RHS);
    aparSolver->setOuterBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD + INVERT_RHS);
    last_time = 0.0;

    apar_boundary_timescale = options["apar_boundary_timescale"]
      .doc("Timescale for Apar boundary relaxation [seconds]")
      .withDefault(1e-8)
      / get<BoutReal>(alloptions["units"]["seconds"]);

  } else if (options["apar_boundary_neumann"]
      .doc("Neumann on all radial boundaries?")
      .withDefault<bool>(false)) {
    // Set zero-gradient (neumann) boundary condition DC on the core
    aparSolver->setInnerBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD);
    aparSolver->setOuterBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD);

  } else if (options["apar_core_neumann"]
      .doc("Neumann radial boundary in the core? False => Dirichlet")
        .withDefault<bool>(true)
             and bout::globals::mesh->periodicY(bout::globals::mesh->xstart)) {
    // Set zero-gradient (neumann) boundary condition DC on the core
    aparSolver->setInnerBoundaryFlags(INVERT_DC_GRAD);
  }

  diagnose = options["diagnose"]
    .doc("Output additional diagnostics?")
    .withDefault<bool>(false);
}

void Electromagnetic::transform(Options &state) {
  AUTO_TRACE();
  
  Options& allspecies = state["species"];

  // Sum coefficients over species
  //
  // Laplace(A||) - alpha_em * beta_em A|| = - beta_em * Ajpar
  alpha_em = 0.0;
  Ajpar = 0.0;
  for (auto& kv : allspecies.getChildren()) {
    const Options& species = kv.second;

    if (!IS_SET(species["charge"]) or !species.isSet("momentum")) {
      continue; // Not charged, or no parallel flow
    }
    const BoutReal Z = get<BoutReal>(species["charge"]);
    if (fabs(Z) < 1e-5) {
      continue; // Not charged
    }

    // Cannot rely on boundary conditions being set
    const Field3D N = GET_NOBOUNDARY(Field3D, species["density"]);
    // Non-final because we're going to change momentum
    const Field3D mom = getNonFinal<Field3D>(species["momentum"]);
    const BoutReal A = get<BoutReal>(species["AA"]);

    // Coefficient in front of A_||
    alpha_em += floor(N, 1e-5) * (SQ(Z) / A);

    // Right hand side
    Ajpar += mom * (Z / A);
  }

  // Invert Helmholtz equation for Apar
  aparSolver->setCoefA((-beta_em) * alpha_em);

  if (const_gradient) {
    // Set gradient boundary condition from gradient inside boundary
    Field3D rhs = (-beta_em) * Ajpar;

    const auto* mesh = Apar.getMesh();
    const auto* coords = Apar.getCoordinates();

    BoutReal time = get<BoutReal>(state["time"]);
    BoutReal weight = 1.0;
    if (time > last_time) {
      weight = exp((last_time - time) / apar_boundary_timescale);
    }
    last_time = time;

    if (mesh->firstX()) {
      const int x = mesh->xstart - 1;
      for (int y = mesh->ystart; y <= mesh->yend; y++) {
        for (int z = mesh->zstart; z <= mesh->zend; z++) {
          rhs(x, y, z) = (weight * (Apar(x + 1, y, z) - Apar(x, y, z)) +
                          (1 - weight) * (Apar(x + 2, y, z) - Apar(x + 1, y, z))) /
            (sqrt(coords->g_11(x, y)) * coords->dx(x, y));
        }
      }
    }
    if (mesh->lastX()) {
      const int x = mesh->xend + 1;
      for (int y = mesh->ystart; y <= mesh->yend; y++) {
        for (int z = mesh->zstart; z <= mesh->zend; z++) {
          rhs(x, y, z) =  (weight * (Apar(x, y, z) - Apar(x - 1, y, z)) +
                           (1 - weight) * (Apar(x - 1, y, z) - Apar(x - 2, y, z))) /
            sqrt(coords->g_11(x, y)) / coords->dx(x, y);
        }
      }
    }
    // Use previous value of Apar as initial guess
    Apar = aparSolver->solve(rhs, Apar);
  } else {
    Apar = aparSolver->solve((-beta_em) * Ajpar, Apar);
  }

  // Save in the state
  set(state["fields"]["Apar"], Apar);

  // Update momentum
  for (auto& kv : allspecies.getChildren()) {
    Options& species = allspecies[kv.first]; // Note: need non-const

    if (!species.isSet("charge") or !species.isSet("momentum")) {
      continue; // Not charged, or no parallel flow
    }
    const BoutReal Z = get<BoutReal>(species["charge"]);
    if (fabs(Z) < 1e-5) {
      continue; // Not charged
    }
    const BoutReal A = get<BoutReal>(species["AA"]);
    const Field3D N = GET_NOBOUNDARY(Field3D, species["density"]);

    Field3D nv = getNonFinal<Field3D>(species["momentum"]);
    nv -= Z * N * Apar;
    // Note: velocity is momentum / (A * N)
    Field3D v = getNonFinal<Field3D>(species["velocity"]);
    v -= (Z / A) * N * Apar / floor(N, 1e-5);
    // Need to update the guard cells
    bout::globals::mesh->communicate(nv, v);
    v.applyBoundary("dirichlet");
    nv.applyBoundary("dirichlet");

    set(species["momentum"], nv);
    set(species["velocity"], v);
  }
}

void Electromagnetic::outputVars(Options &state) {
  // Normalisations
  auto Bnorm = get<BoutReal>(state["Bnorm"]);
  auto rho_s0 = get<BoutReal>(state["rho_s0"]);

  set_with_attrs(state["beta_em"], beta_em, {
     {"long_name", "Helmholtz equation parameter"}
   });

  set_with_attrs(state["Apar"], Apar, {
      {"time_dimension", "t"},
      {"units", "T m"},
      {"conversion", Bnorm * rho_s0},
      {"standard_name", "b dot A"},
      {"long_name", "Parallel component of vector potential A"}
    });

  if (diagnose) {
    set_with_attrs(state["Ajpar"], Ajpar, {
      {"time_dimension", "t"},
    });

    set_with_attrs(state["alpha_em"], alpha_em, {
      {"time_dimension", "t"},
    });
  }
}
