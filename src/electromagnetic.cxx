#include "../include/electromagnetic.hxx"

#include <bout/constants.hxx>
#include <invert_laplace.hxx>

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

  aparSolver = Laplacian::create(&options["laplacian"]);
  // Set zero-gradient (neumann) boundary conditions
  aparSolver->setInnerBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD);
  aparSolver->setOuterBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD);

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
    alpha_em += N * (SQ(Z) / A);

    // Right hand side
    Ajpar += mom * (Z / A);
  }

  // Invert Helmholtz equation for Apar
  aparSolver->setCoefA((-beta_em) * alpha_em);
  Apar = aparSolver->solve((-beta_em) * Ajpar);

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

    subtract(species["momentum"], Z * N * Apar);
    // Note: velocity is momentum / (A * N)
    subtract(species["velocity"], (Z / A) * Apar);
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
