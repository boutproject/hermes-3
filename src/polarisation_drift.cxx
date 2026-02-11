#include "../include/polarisation_drift.hxx"

#include <bout/constants.hxx>
#include <bout/fv_ops.hxx>
#include <bout/output_bout_types.hxx>
#include <bout/invert_laplace.hxx>
#include <bout/difops.hxx>

using bout::globals::mesh;

PolarisationDrift::PolarisationDrift(std::string name, Options& alloptions,
                                     Solver* UNUSED(solver))
    // FIXME: There is a lot of complicated conditional logic which is not being captured
    // here. E.g., species without AA or momentum will be skipped.
    : Component({readIfSet("species:{all_species}:charge"),
                 readOnly("species:{charged}:{inputs}"),
                 readIfSet("species:{charged}:momentum"),
                 readIfSet("species:{charged}:pressure", Regions::Interior),
                 readIfSet("fields:{fields}"),
                 readWrite("species:{charged}:{outputs}")}) {

  // Get options for this component
  auto& options = alloptions[name];

  // Cache the B^2 value
  auto coord = mesh->getCoordinates();
  Bsq = SQ(coord->Bxy);

  phiSolver = Laplacian::create(&options["laplacian"]);

  // For zero polarisation drift fluxes through the boundary the
  // radial electric field at the boundary should be constant
  // (e.g. zero), so zero-gradient radial BC on phi_G.

  phiSolver->setInnerBoundaryFlags(0);
  phiSolver->setOuterBoundaryFlags(0);

  boussinesq = options["boussinesq"]
    .doc("Assume a uniform mass density in calculating the polarisation drift")
    .withDefault<bool>(true);

  if (boussinesq) {
    average_atomic_mass =
      options["average_atomic_mass"]
      .doc("Weighted average atomic mass, for polarisaion current "
           "(Boussinesq approximation)")
      .withDefault<BoutReal>(2.0); // Deuterium
  } else {
    average_atomic_mass = 1.0;
    // Use a density floor to prevent divide-by-zero errors
    density_floor = options["density_floor"].doc("Minimum density floor").withDefault(1e-5);
  }

  advection = options["advection"]
    .doc("Include advection by polarisation drift, using potential flow approximation?")
    .withDefault<bool>(true);

  diamagnetic_polarisation =
      options["diamagnetic_polarisation"]
          .doc("Include diamagnetic drift in polarisation current?")
          .withDefault<bool>(true);

  diagnose = options["diagnose"]
    .doc("Output additional diagnostics?")
    .withDefault<bool>(false);

  // Pressure interior only,
  std::vector<std::string> inputs = {"AA"}, fields = {"DivJdia"},
                           outputs = {"energy_source"};
  if (advection) {
    // Interior only
    inputs.push_back("density");
    fields.push_back("DivJextra");
    fields.push_back("DivJdia");
    outputs.push_back("density_source");
    outputs.push_back("momentum_source");
  }
  substitutePermissions("inputs", inputs);
  substitutePermissions("fields", fields);
  // FIXME: energy_source and momentum source are only set if pressure
  // and momentum were set, respectively
  substitutePermissions("outputs", outputs);
}

Field3D PolarisationDrift::calcDivJ(GuardedOptions& state) {
  // Iterate through all subsections
  GuardedOptions allspecies = state["species"];

  // Calculate divergence of all currents except the polarisation current
  Field3D DivJ;
  if (IS_SET(state["fields"]["DivJdia"])) {
    DivJ = get<Field3D>(state["fields"]["DivJdia"]);
  } else {
    DivJ = 0.0;
  }

  // Parallel current due to species parallel flow
  for (auto& kv : allspecies.getChildren()) {
    const GuardedOptions species = kv.second;

    if (!species.isSet("charge") or !species.isSet("momentum")) {
      continue; // Not charged, or no parallel flow
    }
    const BoutReal Z = get<BoutReal>(species["charge"]);
    if (fabs(Z) < 1e-5) {
      continue; // Not charged
    }

    const Field3D NV = GET_VALUE(Field3D, species["momentum"]);
    const BoutReal A = get<BoutReal>(species["AA"]);

    // Note: Using NV rather than N*V so that the cell boundary flux is correct
    DivJ += Div_par((Z / A) * NV);
  }

  if (IS_SET(state["fields"]["DivJextra"])) {
    DivJ += get<Field3D>(state["fields"]["DivJextra"]);
  }
  if (IS_SET(state["fields"]["DivJcol"])) {
    DivJ += get<Field3D>(state["fields"]["DivJcol"]);
  }

  return DivJ;
}

void PolarisationDrift::diamagneticCompression(GuardedOptions& state, Field3D DivJ) {
  // Compression of ion diamagnetic contribution to polarisation velocity
  // Note: For now this ONLY includes the parallel and diamagnetic current terms
  //       Other terms e.g. ion viscous current, are in their separate components
  if (!boussinesq) {
    throw BoutException("diamagnetic_polarisation not implemented for non-Boussinesq");
  }

  GuardedOptions allspecies = state["species"];

  // Calculate energy exchange term nonlinear in pressure
  // (3 / 2) ddt(Pi) += (Pi / n0) * Div((Pe + Pi) * Curlb_B + Jpar);
  for (auto& kv : allspecies.getChildren()) {
    GuardedOptions species = allspecies[kv.first]; // Note: need non-const

    if (!(IS_SET_NOBOUNDARY(species["pressure"]) and species.isSet("charge")
          and species.isSet("AA"))) {
      // No pressure, charge or mass -> no polarisation current due to
      // diamagnetic flow
      continue;
    }

    const auto charge = get<BoutReal>(species["charge"]);
    if (fabs(charge) < 1e-5) {
      // No charge
      continue;
    }

    const auto P = GET_NOBOUNDARY(Field3D, species["pressure"]);
    const auto AA = get<BoutReal>(species["AA"]);

    Field3D energy_source = P * (AA / average_atomic_mass / charge) * DivJ;

    if (diagnose) {
      set_with_attrs(diagnostics[fmt::format("E{}_pol", kv.first)], energy_source,
                     {{"time_dimension", "t"},
                      {"units", "W m^-3"},
                      //{"conversion", SI::qe * Tnorm * Nnorm * Omega_ci},
                      {"long_name", "Compression of polarisation drift"},
                      {"source", "polarisation_drift"}});
    }

    add(species["energy_source"], energy_source);
  }
}

Field3D PolarisationDrift::calcMassDensity(GuardedOptions& state) {
  // Calculate the total mass density of species
  // which contribute to polarisation current
  Field3D mass_density;
  if (boussinesq) {
    mass_density = average_atomic_mass;
  } else {
    mass_density = 0.0;
    GuardedOptions allspecies = state["species"];
    for (auto& kv : allspecies.getChildren()) {
      const GuardedOptions species = kv.second;

      if (!(species.isSet("charge") and species.isSet("AA") and species.isSet("density"))) {
        continue; // No charge or mass -> no current
      }
      if (fabs(get<BoutReal>(species["charge"])) < 1e-5) {
        continue; // Zero charge
      }

      const BoutReal A = get<BoutReal>(species["AA"]);
      const Field3D N = GET_NOBOUNDARY(Field3D, species["density"]);
      mass_density += A * N;
    }

    // Apply a floor to prevent divide-by-zero errors
    mass_density = floor(mass_density, density_floor);
  }
  return mass_density;
}

Field3D PolarisationDrift::calcPolFlowPotential(Field3D mass_density, Field3D DivJ) {
  phiSolver->setCoefC(mass_density / Bsq);

  // Calculate time derivative of generalised potential
  // The assumption is that the polarisation drift can be parameterised
  // as a potential flow
  Field3D phi_pol = phiSolver->solve(DivJ * Bsq / mass_density);

  // Ensure that potential is set in communication guard cells
  mesh->communicate(phi_pol);

  // Zero flux
  phi_pol.applyBoundary("neumann");

  return phi_pol;
}

void PolarisationDrift::polarisationAdvection(GuardedOptions& state, Field3D phi_pol) {
  GuardedOptions allspecies = state["species"];

  for (auto& kv : allspecies.getChildren()) {
    GuardedOptions species = allspecies[kv.first]; // Note: need non-const

    if (!(species.isSet("charge") and species.isSet("AA"))) {
      continue; // No charge or mass -> no current
    }
    const BoutReal Z = get<BoutReal>(species["charge"]);
    if (fabs(Z) < 1e-5) {
      continue; // Not charged
    }
    const BoutReal A = get<BoutReal>(species["AA"]);

    // Shared coefficient in polarisation velocity
    Field3D coef = (A / Z) / Bsq;

    if (IS_SET(species["density"])) {
      auto N = GET_VALUE(Field3D, species["density"]);
      add(species["density_source"], FV::Div_a_Grad_perp(N * coef, phi_pol));
    }

    if (IS_SET(species["pressure"])) {
      auto P = GET_VALUE(Field3D, species["pressure"]);
      add(species["energy_source"], (5. / 2) * FV::Div_a_Grad_perp(P * coef, phi_pol));
    }

    if (IS_SET(species["momentum"])) {
      auto NV = GET_VALUE(Field3D, species["momentum"]);
      add(species["momentum_source"], FV::Div_a_Grad_perp(NV * coef, phi_pol));
    }
  }
}

void PolarisationDrift::transform_impl(GuardedOptions& state) {

  // Calculate divergence of all current except polarisation
  DivJ = calcDivJ(state);

  // Iterate through all subsections
  GuardedOptions allspecies = state["species"];

  if (diamagnetic_polarisation) {
    // Calculate energy source due to compression of polarisation drift
    diamagneticCompression(state, DivJ);
  }

  if (!advection) {
    return;
  }

  // Calculate advection terms using a potential-flow approximation
  Field3D mass_density = calcMassDensity(state);
  phi_pol = calcPolFlowPotential(mass_density, DivJ);
  // Calculate the advection terms for each species
  polarisationAdvection(state, phi_pol);
}

void PolarisationDrift::outputVars(Options& state) {
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);

  if (diagnose) {
    set_with_attrs(state["DivJpol"], -DivJ,
                   {{"time_dimension", "t"},
                    {"units", "A m^-3"},
                    {"conversion", SI::qe * Nnorm * Omega_ci},
                    {"long_name", "Divergence of polarisation current"},
                    {"source", "polarisation_drift"}});

    if (advection) {
      set_with_attrs(state["phi_pol"], phi_pol,
                     {{"time_dimension", "t"},
                      {"units", "V / s"},
                      {"conversion", Tnorm * Omega_ci},
                      {"standard_name", "flow potential"},
                      {"long_name", "polarisation flow potential"},
                      {"source", "polarisation_drift"}});
    }

    // Copy diagnostics into output
    for (auto& kv : diagnostics.getChildren()) {
      state[kv.first] = kv.second.copy();
    }
  }
}
