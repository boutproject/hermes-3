#include "../include/polarisation_drift.hxx"
#include "../include/component.hxx"
#include "../include/guarded_options.hxx"
#include "../include/hermes_utils.hxx"
#include "../include/permissions.hxx"

#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/constants.hxx>
#include <bout/difops.hxx>
#include <bout/field.hxx>
#include <bout/field3d.hxx>
#include <bout/fv_ops.hxx>
#include <bout/globals.hxx>
#include <bout/invert_laplace.hxx>
#include <bout/mesh.hxx>
#include <bout/output_bout_types.hxx>
#include <bout/utils.hxx>

#include "../include/div_ops.hxx"
#include <bout/vecops.hxx>
#include <algorithm>
#include <limits>
#include <tuple>

#include <string>
#include <vector>

using bout::globals::mesh;

PolarisationDrift::PolarisationDrift(std::string name, Options& alloptions,
                                     Solver* UNUSED(solver))
    : NamedComponent(name, {}) {
  auto& options = alloptions[name];
  implementation = options["implementation"]
                       .doc("Polarisation implementation: legacy potential flow or "
                            "conservative ion energy transport")
                       .withDefault(bout::polarisation::Implementation::legacy);
  std::vector<std::string> ion_models;
  std::vector<std::string> enabled;
  const Options& scheduler =
      alloptions.isSet("components") ? alloptions : alloptions["hermes"];
  if (scheduler.isSet("components")) {
    for (const auto& item : strsplit(scheduler["components"].as<std::string>(), ',')) {
      enabled.push_back(trim(item, " \t\r\n()"));
    }
  }
  for (const auto& [component_name, section] : alloptions.getChildren()) {
    if (!enabled.empty()
        && std::find(enabled.begin(), enabled.end(), component_name) == enabled.end()) {
      continue;
    }
    const std::string component_type =
        section.isSet("type") ? section["type"].as<std::string>() : component_name;
    if (component_type == "evolve_vorticity_ion_energy") {
      ion_models.push_back(component_name);
    }
  }
  if (implementation == bout::polarisation::Implementation::conservative) {
    if (ion_models.size() != 1) {
      throw BoutException("Conservative polarisation_drift requires exactly one "
                          "evolve_vorticity_ion_energy component");
    }
    auto& ion = alloptions[ion_models.front()];
    for (const auto* key : {"boussinesq", "average_atomic_mass", "density_floor",
                            "advection", "diamagnetic_polarisation"}) {
      if (options.isSet(key)) {
        throw BoutException("polarisation_drift:{} is a legacy option; conservative "
                            "transport uses ion-energy recovery settings",
                            key);
      }
    }
    average_atomic_mass = ion["average_atomic_mass"].withDefault<BoutReal>(2.0);
    density_floor = ion["density_floor"].withDefault<BoutReal>(1e-7);
    polarisation_density_threshold =
        options["polarisation_density_threshold"]
            .doc("Density threshold in units of Nnorm for polarisation particle "
                 "transport; zero disables the limiter")
            .withDefault<BoutReal>(0.0);
    if (!std::isfinite(polarisation_density_threshold)
        || polarisation_density_threshold < 0.0) {
      throw BoutException("polarisation_density_threshold must be finite and >= 0");
    }
    poloidal_flows = ion["poloidal_flows"].withDefault<bool>(false);
    diagnose = options["diagnose"].withDefault(ion["diagnose"].withDefault<bool>(false));
    radial_model.walls = !mesh->periodicX;
    radial_model.chi_ac =
        options["polarisation_ac_boundary"]
            .doc("Auxiliary AC wall condition: neumann (default) or dirichlet; "
                 "DC uses the phi wall rate, neglecting the pressure-potential "
                 "wall rate")
            .withDefault(bout::ion_energy::PolarisationACBoundary::neumann);
    radial_model.chi_flags =
        radial_model.chi_ac == bout::ion_energy::PolarisationACBoundary::neumann
            ? INVERT_AC_GRAD
            : 0;
    polarisation_corrections =
        options["polarisation_corrections"]
            .doc("Fixed FV defect corrections after the initial spectral solve")
            .withDefault(10);
    if (polarisation_corrections < 0) {
      throw BoutException("polarisation_corrections must be >= 0");
    }
    auto& laplacian = options["laplacian"];
    if (!radial_model.walls && laplacian.getChildren().empty()) {
      laplacian = ion["laplacian"].copy();
    }
    for (const auto* key : {"inner_boundary_flags", "outer_boundary_flags"}) {
      if (laplacian.isSet(key)) {
        throw BoutException("polarisation_drift:laplacian:{} is controlled by "
                            "radial_boundary and polarisation_ac_boundary",
                            key);
      }
    }
    // Components are constructed in input order, before permission autosorting.
    // Resolve the factory default here without requiring ion recovery to have
    // constructed its inverse. An explicit drift type needs no ion-type read.
    if (!laplacian.isSet("type")) {
      laplacian["type"].withDefault(ion["laplacian"]["type"].withDefault(
          std::string{LaplaceFactory::default_type}));
    }
    if (radial_model.walls
        && (laplacian["global_flags"].withDefault(0)
            & (INVERT_ZERO_DC | INVERT_KX_ZERO))) {
      throw BoutException("polarisation_drift:laplacian must retain the DC solution");
    }
    phiSolver = Laplacian::create(&laplacian);
    if (radial_model.walls) {
      phiSolver->setInnerBoundaryFlags(radial_model.chi_flags);
      phiSolver->setOuterBoundaryFlags(radial_model.chi_flags);
    }
    wall_rate_inner.resize(static_cast<size_t>(mesh->LocalNy), 0.0);
    wall_rate_outer.resize(static_cast<size_t>(mesh->LocalNy), 0.0);
    setPermissions(readOnly("species:{all_species}:charge"));
    setPermissions(readOnly("species:{ions}:AA"));
    for (const auto* field : {"density", "pressure", "temperature"}) {
      setPermissions(readOnly(std::string{"species:{ions}:"} + field, Regions::Interior));
    }
    setPermissions(readIfSet("species:{ions}:velocity", Regions::Interior));
    // Boundary-final momentum writers (including the production sheath) must
    // run after recovery and before current continuity. Recovery reads interiors.
    setPermissions(readIfSet("species:{charged}:momentum", Regions::All));
    setPermissions(readIfSet("species:{charged}:AA"));
    setPermissions(readOnly("fields:phi", Regions::All));
    setPermissions(readIfSet("fields:DivJextra", Regions::Interior));
    setPermissions(readIfSet("fields:DivJcol", Regions::Interior));
    setPermissions(readOnly("fields:ion_boussinesq_coefficient"));
    for (const auto* field :
         {"vorticity", "ion_pressure_potential", "ion_mass_density", "ion_reynolds_x",
          "ion_reynolds_z", "DivJReynolds", "DivJdia", "ion_wall_rate"}) {
      setPermissions(readOnly(std::string{"fields:"} + field, Regions::Interior));
    }
    for (const auto* field :
         {"density_source", "momentum_source", "polarisation_density_source",
          "polarisation_momentum_source"}) {
      setPermissions(readWrite(std::string{"species:{ions}:"} + field));
    }
    for (const auto* field : {"ion_polarisation_energy_divergence", "DivJpol", "ion_chi",
                              "DivJpar", "ion_current_rhs", "ion_auxiliary_residual"}) {
      setPermissions(readWrite(std::string{"fields:"} + field));
    }
    return;
  }
  if (!ion_models.empty()) {
    throw BoutException(
        "Legacy polarisation_drift cannot pair with evolve_vorticity_ion_energy; select "
        "implementation=conservative");
  }
  setPermissions(readIfSet("species:{all_species}:charge"));
  setPermissions(readOnly("species:{charged}:{inputs}"));
  setPermissions(readIfSet("species:{charged}:momentum"));
  setPermissions(readIfSet("fields:{fields}"));
  setPermissions(readWrite("species:{charged}:{outputs}"));
  // Cache the B^2 value
  const auto* coord = mesh->getCoordinates();
  Bsq = SQ(coord->Bxy());

  phiSolver = Laplacian::create(&options["laplacian"]);

  // For zero polarisation drift fluxes through the boundary the
  // radial electric field at the boundary should be constant
  // (e.g. zero), so zero-gradient radial BC on phi_G.

  phiSolver->setInnerBoundaryFlags(0);
  phiSolver->setOuterBoundaryFlags(0);

  boussinesq =
      options["boussinesq"]
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
    density_floor =
        options["density_floor"].doc("Minimum density floor").withDefault(1e-5);
  }

  advection = options["advection"]
                  .doc("Include advection by polarisation drift, using potential flow "
                       "approximation?")
                  .withDefault<bool>(true);

  diamagnetic_polarisation =
      options["diamagnetic_polarisation"]
          .doc("Include diamagnetic drift in polarisation current?")
          .withDefault<bool>(true);

  diagnose =
      options["diagnose"].doc("Output additional diagnostics?").withDefault<bool>(false);

  std::vector<std::string> inputs = {"AA"};
  std::vector<std::string> fields = {"DivJdia", "DivJextra", "DivJcol"};
  std::vector<std::string> outputs = {};
  if (advection) {
    inputs.push_back("density");
    outputs.push_back("density_source");
    outputs.push_back("energy_source");
    outputs.push_back("momentum_source");
    setPermissions(readIfSet("species:{charged}:pressure", Regions::All));
  } else if (diamagnetic_polarisation) {
    outputs.push_back("energy_source");
    setPermissions(readIfSet("species:{charged}:pressure", Regions::Interior));
  }
  substitutePermissions("inputs", inputs);
  substitutePermissions("fields", fields);
  // FIXME: density_source, energy_source, and momentum_source are only set if density, pressure,
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

    const Field3D energy_source = P * (AA / average_atomic_mass / charge) * DivJ;

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

      if (!(species.isSet("charge") and species.isSet("AA")
            and species.isSet("density"))) {
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
    const Field3D coef = Field3D{(A / Z) / Bsq};

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
  if (implementation == bout::polarisation::Implementation::conservative) {
    transformConservative(state);
    return;
  }

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
  set(state["polarisation_drift_implementation"], toString(implementation));
  if (implementation == bout::polarisation::Implementation::conservative) {
    outputConservative(state);
    return;
  }
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

namespace {
template <typename Field, typename Value>
void fillRadialGuards(Field& field, Value value) {
  auto* localmesh = field.getMesh();
  if (localmesh->periodicX) {
    return;
  }
  const auto fill = [&](int interior, int step, int count) {
    for (int y = localmesh->ystart; y <= localmesh->yend; ++y) {
      for (int z = 0; z < localmesh->LocalNz; ++z) {
        for (int layer = 1; layer <= count; ++layer) {
          field(interior + step * layer, y, z) = value(interior, step, layer, y, z);
        }
      }
    }
  };
  if (localmesh->firstX()) {
    fill(localmesh->xstart, -1, localmesh->xstart);
  }
  if (localmesh->lastX()) {
    fill(localmesh->xend, 1, localmesh->LocalNx - localmesh->xend - 1);
  }
}

template <typename Function>
void forAllIons(GuardedOptions allspecies, Function func) {
  for (auto& kv : allspecies.getChildren()) {
    if (kv.first == "e") {
      continue; // Skip electrons
    }
    auto species = allspecies[kv.first];

    if (!IS_SET(species["charge"])) {
      continue;
    }
    if (fabs(get<BoutReal>(species["charge"])) < 1e-5) {
      // No charge
      continue;
    }
    func(species);
  }
}
} // namespace

namespace bout::ion_energy {
/// Reflect auxiliary-potential guards about a prescribed DC wall midpoint.
///
/// @param[in,out] chi Auxiliary potential whose physical X-wall guards are updated.
/// @param ac Neumann retains each mirrored interior layer's Z fluctuations;
///           Dirichlet negates them to impose a zero AC midpoint.
/// @param dc_inner Optional inner-wall DC midpoints in chi's units, indexed by local
///                 Y with LocalNy entries; nullptr imposes a zero DC midpoint.
/// @param dc_outer Optional outer-wall DC midpoints with the same convention.
///
/// Returns no value. Each radial guard reflects its corresponding interior layer's
/// mean about the DC value. All local Y rows, including guards, are processed;
/// interior cells and the supplied DC arrays are preserved. Periodic X is a no-op.
/// No communication occurs.
/// @throws BoutException If nonperiodic X does not have complete local periodic Z lines.
void applyPolarisationBoundary(Field3D& chi, PolarisationACBoundary ac,
                               const std::vector<BoutReal>* dc_inner,
                               const std::vector<BoutReal>* dc_outer) {
  auto* localmesh = chi.getMesh();
  if (localmesh->periodicX) {
    return;
  }
  if (localmesh->getNZPE() != 1 || localmesh->zstart != 0
      || localmesh->zend + 1 != localmesh->LocalNz) {
    throw BoutException(
        "Mixed polarisation boundaries require complete periodic Z lines");
  }
  const auto reflect = [&](int guard, int interior, const std::vector<BoutReal>* dc) {
    for (int y = 0; y < localmesh->LocalNy; ++y) {
      BoutReal mean = 0.0;
      for (int z = 0; z < localmesh->LocalNz; ++z) {
        mean += chi(interior, y, z);
      }
      mean /= localmesh->LocalNz;
      for (int z = 0; z < localmesh->LocalNz; ++z) {
        const BoutReal fluctuation = chi(interior, y, z) - mean;
        chi(guard, y, z) =
            -mean + (ac == PolarisationACBoundary::neumann ? fluctuation : -fluctuation);
        if (dc != nullptr) {
          chi(guard, y, z) += 2.0 * (*dc)[static_cast<size_t>(y)];
        }
      }
    }
  };
  if (localmesh->firstX()) {
    for (int x = 0; x < localmesh->xstart; ++x) {
      reflect(x, 2 * localmesh->xstart - 1 - x, dc_inner);
    }
  }
  if (localmesh->lastX()) {
    for (int x = localmesh->xend + 1; x < localmesh->LocalNx; ++x) {
      reflect(x, 2 * localmesh->xend + 1 - x, dc_outer);
    }
  }
}
} // namespace bout::ion_energy

/// Approximate the auxiliary potential using a fixed number of defect corrections.
///
/// @param coefficient Read-only Boussinesq coefficient c in div(c*grad_perp(chi)),
///                    nonzero on physical cells, with valid geometry and guards.
/// @param rhs Read-only selected current-continuity RHS; only interior values enter the
///            residual, excluding encoded wall-potential rates in radial guards.
/// @param[in,out] inverse Approximate spectral inverse configured for coefficient
///                       and homogeneous correction boundaries. Its solve method
///                       is invoked once initially and once per defect correction.
/// @param corrections Nonnegative number of finite-volume defect corrections.
/// @param ac Radial condition on the Z-dependent part of chi.
/// @param dc_inner Optional inner-wall DC midpoints in chi's units, indexed by local
///                 Y with LocalNy entries; nullptr imposes zero.
/// @param dc_outer Optional outer-wall DC midpoints with the same convention.
/// @returns A new chi field, starting from zero on every call, with imposed radial
///          guards and communicated inter-rank guards. A fixed correction count
///          can leave a finite forward-operator residual; no warm start is retained.
///
/// Uses the global mesh to communicate temporary chi and correction fields;
/// mesh face-region caches and the inverse's internal solver state may change.
/// Supplied fields and DC arrays are preserved. Wall values enter through each
/// finite-volume residual while spectral corrections use homogeneous boundaries.
/// @throws BoutException If corrections is negative or a nonperiodic X mesh lacks
///                       complete local periodic Z lines.
Field3D bout::ion_energy::fixedPolarisationPotential(
    const Field3D& coefficient, const Field3D& rhs, Laplacian& inverse, int corrections,
    PolarisationACBoundary ac, const std::vector<BoutReal>* dc_inner,
    const std::vector<BoutReal>* dc_outer) {
  if (corrections < 0) {
    throw BoutException("polarisation_corrections must be >= 0");
  }
  Field3D chi = zeroFrom(rhs);
  const auto apply_boundary = [&]() {
    mesh->communicate(chi);
    bout::ion_energy::applyPolarisationBoundary(chi, ac, dc_inner, dc_outer);
    mesh->communicate(chi);
  };

  // Start from zero every call, and use a fixed number of solves rather than a
  // tolerance, so chi is a smooth function of rhs (and so of the state) and
  // never depends on previous RHS evaluations. The first solve is the spectral
  // inverse; the rest are defect corrections to the finite-volume operator.
  for (int solve = 0; solve <= corrections; ++solve) {
    apply_boundary();
    Field3D residual{zeroFrom(rhs)};
    const Field3D forward = Div_a_Grad_perp_XZ(coefficient, chi, true);
    BOUT_FOR(i, residual.getRegion("RGN_NOBNDRY")) {
      residual[i] = (rhs[i] - forward[i]) / coefficient[i];
    }
    const Field3D approximate = inverse.solve(residual, zeroFrom(chi));
    // Only copy solved X-Z planes; unused Y guards need deterministic values.
    Field3D correction{zeroFrom(chi)};
    BOUT_FOR(i, correction.getRegion("RGN_NOY")) { correction[i] = approximate[i]; }
    mesh->communicate(correction);
    chi += correction;
  }
  apply_boundary();
  return chi;
}

void PolarisationDrift::solveConservativePotential(const Field3D& rhs) {
  if (mesh->periodicX) {
    // Each perpendicular plane has a constant null mode. A parallel-current
    // source with a nonzero volume-weighted mean cannot be inverted here.
    // Reduce X/Z over the corresponding X communicator, separately for each Y.
    const auto* coord = rhs.getCoordinates();
    std::vector<BoutReal> local(2 * mesh->LocalNy, 0.0), global(local.size(), 0.0);
    BOUT_FOR(i, rhs.getRegion("RGN_NOBNDRY")) {
      const BoutReal weighted =
          rhs[i] * coord->J()[i] * coord->dx()[i] * coord->dy()[i] * coord->dz()[i];
      local[2 * i.y()] += weighted;
      local[2 * i.y() + 1] += std::abs(weighted);
    }
    MPI_Allreduce(local.data(), global.data(), static_cast<int>(local.size()), MPI_DOUBLE,
                  MPI_SUM, mesh->getXcomm());
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      const BoutReal tolerance =
          128 * std::numeric_limits<BoutReal>::epsilon() * global[2 * y + 1];
      if (std::abs(global[2 * y]) > tolerance) {
        throw BoutException("Conservative polarisation has an unsolvable periodic "
                            "X-Z mean current source at local Y={}: integral={}, "
                            "roundoff bound={}",
                            y, global[2 * y], tolerance);
      }
    }
  }
  // Preserve the existing periodic inverse/gauge; physical X uses the
  // separate mixed-boundary inverse, matching applyPolarisationBoundary.
  Laplacian& inverse = *phiSolver;
  // chi approximates d(Phi)/dt, whose wall value is d(phi_w)/dt + d(psi_w(z))/dt.
  // The closure uses d(phi_w)/dt for DC and neglects the pressure-potential
  // wall rate; AC defaults to Neumann. Corrections retain homogeneous boundaries.
  // The inhomogeneous value enters through every FV residual, including the first.
  chi = bout::ion_energy::fixedPolarisationPotential(
      boussinesq_coeff, rhs, inverse, polarisation_corrections, radial_model.chi_ac,
      radial_model.walls ? &wall_rate_inner : nullptr,
      radial_model.walls ? &wall_rate_outer : nullptr);
  chi_iterations = 1 + polarisation_corrections;

  // Diagnostic only: residual against the finite-volume operator
  auxiliary_residual = rhs - bout::ion_energy::Div_a_Grad_perp_XZ(boussinesq_coeff, chi);
  chi_residual = max(abs(auxiliary_residual), true);
}

BoutReal PolarisationDrift::densityLimiter(BoutReal mass_density) const {
  if (polarisation_density_threshold == 0.0) {
    return 1.0; // softFloor requires a strictly positive threshold.
  }
  const BoutReal rho = std::max(mass_density, 0.0);
  return rho / softFloor(rho, average_atomic_mass * polarisation_density_threshold);
}

void PolarisationDrift::transformConservative(GuardedOptions& state) {
  balance.clear();
  auto allspecies = state["species"];
  auto fields = state["fields"];
  const Field3D phi = GET_VALUE(Field3D, fields["phi"]);
  const Field3D psi = GET_NOBOUNDARY(Field3D, fields["ion_pressure_potential"]);
  const Field3D vorticity = GET_NOBOUNDARY(Field3D, fields["vorticity"]);
  boussinesq_coeff = get<Coordinates::FieldMetric>(fields["ion_boussinesq_coefficient"]);
  phiSolver->setCoefC(boussinesq_coeff);
  Field3D mass_density = GET_NOBOUNDARY(Field3D, fields["ion_mass_density"]);
  const Field3D DivJReynolds = GET_NOBOUNDARY(Field3D, fields["DivJReynolds"]);
  const Field3D DivJdia = GET_NOBOUNDARY(Field3D, fields["DivJdia"]);
  const bout::ion_energy::XZFaceFluxes reynolds{
      GET_NOBOUNDARY(Field3D, fields["ion_reynolds_x"]),
      GET_NOBOUNDARY(Field3D, fields["ion_reynolds_z"])};
  const Field3D rates = GET_NOBOUNDARY(Field3D, fields["ion_wall_rate"]);
  if (radial_model.walls) {
    for (int y = 0; y < mesh->LocalNy; ++y) {
      wall_rate_inner[static_cast<size_t>(y)] = rates(mesh->xstart - 1, y, 0);
      wall_rate_outer[static_cast<size_t>(y)] = rates(mesh->xend + 1, y, 0);
    }
  }
  //////////////////////////////////////////////////////////
  // Polarisation drift

  // Auxiliary current-continuity potential, independent of phi relaxation.
  // DC Dirichlet fixes the constant mode; the AC enum controls the auxiliary
  // potential guards. Transport retains the currents used by the forward solve.
  // Excluding the viscous rate from chi means the transported polarisation
  // current contains the ion viscous drift current. Its work and enthalpy are
  // carried by DivEpol, so subtracting DivEviscous separately would duplicate
  // that work. This relation uses the actual DivJviscous for either stress-wall
  // setting: free_slip retains viscous X faces, fixed/no_slip close them, and
  // polarisation transport retains its own physical X-wall faces in both cases.
  Jpar = zeroFrom(phi);
  for (auto& [name, species] : allspecies.getChildren()) {
    if (!species.isSet("charge") || !species.isSet("momentum")) {
      continue;
    }
    const BoutReal charge = get<BoutReal>(species["charge"]);
    if (std::abs(charge) < 1e-5) {
      continue;
    }
    const BoutReal mass = get<BoutReal>(species["AA"]);
    Jpar += (charge / mass) * GET_VALUE(Field3D, species["momentum"]);
  }
  mesh->communicate(Jpar);
  DivJpar = Div_par(Jpar);
  current_rhs = DivJReynolds + DivJdia + DivJpar;
  for (const auto* optional : {"DivJextra", "DivJcol"}) {
    if (IS_SET_NOBOUNDARY(fields[optional])) {
      current_rhs += GET_NOBOUNDARY(Field3D, fields[optional]);
    }
  }
  solveConservativePotential(current_rhs);

  // Use the same geometric chi currents as the residual, including wall faces.
  const auto chi_faces = bout::ion_energy::perpendicularFaces(boussinesq_coeff, chi);
  polarisation_flux = {-reynolds.x - chi_faces.x, -reynolds.z - chi_faces.z};
  const auto divPol = [&](const Field3D& weight) -> Field3D {
    Field3D result = bout::ion_energy::divergence(polarisation_flux, &weight);
    if (poloidal_flows) {
      Field3D w = copy(weight);
      fillRadialGuards(w, [&](int x, int, int, int y, int z) { return w(x, y, z); });
      result -= Div_n_bxGrad_f_B_XPPM(vorticity, phi, false, true, false, &w)
                - Div_n_bxGrad_f_B_XPPM(vorticity, phi, false, false, false, &w);
    }
    return result;
  };
  // Energy, momentum and particles carried through a radial wall use the local
  // wall-face state. Neumann guards make these values equal to the adjacent cell.
  // Work uses the local Phi midpoint, including conducting-wall pressure work.
  const auto midpoint = [](const Field3D& value, int x, int y, int z) {
    return 0.5 * (value(x, y, z) + value(x + 1, y, z));
  };
  const auto* coord = phi.getCoordinates();
  Field3D wall_energy_weight = zeroFrom(phi);
  Field3D wall_work_weight = zeroFrom(phi);
  Field3D wall_thermal_weight = zeroFrom(phi);
  Field3D wall_pressure_mismatch = zeroFrom(phi);
  const auto wallWeighted = [&](const Field3D& weight,
                                const Field3D* prescribed = nullptr) -> Field3D {
    Field3D result = divPol(weight);
    for (int x : {mesh->xstart - 1, mesh->xend}) {
      const bool inner = x == mesh->xstart - 1;
      if (mesh->periodicX || (inner && !mesh->firstX()) || (!inner && !mesh->lastX())) {
        continue;
      }
      const int cell = inner ? x + 1 : x;
      for (int y = mesh->ystart; y <= mesh->yend; ++y) {
        for (int z = 0; z < mesh->LocalNz; ++z) {
          const BoutReal revised =
              prescribed ? (*prescribed)(x, y, z) : midpoint(weight, x, y, z);
          const BoutReal volume = coord->J()(cell, y, z) * coord->dx()(cell, y, z)
                                  * coord->dy()(cell, y, z) * coord->dz()(cell, y, z);
          result(cell, y, z) += (inner ? -1.0 : 1.0) * polarisation_flux.x(x, y, z)
                                * (revised - weight(cell, y, z)) / volume;
        }
      }
    }
    return result;
  };
  DivJpol = divPol(Field3D{1.0});
  bout::ion_energy::XZFaceFluxes wall_flux;
  if (diagnose) {
    wall_flux = {zeroFrom(phi), zeroFrom(phi)};
    for (int x : {mesh->xstart - 1, mesh->xend}) {
      if (!mesh->periodicX
          && ((x == mesh->xstart - 1 && mesh->firstX())
              || (x == mesh->xend && mesh->lastX()))) {
        for (int y = mesh->ystart; y <= mesh->yend; ++y) {
          for (int z = 0; z < mesh->LocalNz; ++z) {
            wall_flux.x(x, y, z) = polarisation_flux.x(x, y, z);
          }
        }
      }
    }
  }

  // Jps = (A_s n_s / rho) Jpol, Gamma_ps = Jps / Z_s. Use the actual
  // mass density, so the species currents sum to Jpol above the floor.
  mass_density = floor(mass_density, average_atomic_mass * density_floor);
  const auto wallFraction = [&](const Field3D& density, BoutReal mass, int x, int y,
                                int z) {
    return 0.5
           * (mass * density(x, y, z) / mass_density(x, y, z)
              + mass * density(x + 1, y, z) / mass_density(x + 1, y, z));
  };
  // The pressure potential uses the same reference mass normalisation as Phi
  // and K_perp. Transporting work as Phi Jpol keeps the local thermal/kinetic
  // partition independent of a constant potential gauge up to the auxiliary
  // current residual. Thermal and parallel kinetic weights below retain the
  // actual species mass-density partition.
  Field3D E_adv = phi + psi;
  diagnostics = Options{};
  forAllIons(allspecies, [&](GuardedOptions species) {
    const BoutReal AA = get<BoutReal>(species["AA"]);
    const BoutReal Z = get<BoutReal>(species["charge"]);
    const Field3D N = GET_NOBOUNDARY(Field3D, species["density"]);
    const Field3D P = GET_NOBOUNDARY(Field3D, species["pressure"]);
    const Field3D coefficient = (AA / Z) / mass_density;
    const Field3D density_source = -wallWeighted(coefficient * N);
    balance.set("particles_" + species.name(), "out_polarisation", -density_source);
    add(species["density_source"], density_source);
    set(species["polarisation_density_source"], density_source);
    if (diagnose) {
      const Field3D weight = coefficient * N;
      diagnostics[species.name()]["wall_loss"] =
          wallWeighted(weight) - divPol(weight)
          + bout::ion_energy::divergence(wall_flux, &weight);
    }
    E_adv += coefficient * (3. / 2) * P;
    if (species.isSet("velocity")) {
      const Field3D V = GET_NOBOUNDARY(Field3D, species["velocity"]);
      E_adv += coefficient * 0.5 * AA * N * SQ(V);
    }
    if (species.isSet("momentum")) {
      const Field3D NV = GET_NOBOUNDARY(Field3D, species["momentum"]);
      Field3D wall_momentum_weight = zeroFrom(phi);
      const Field3D V = GET_NOBOUNDARY(Field3D, species["velocity"]);
      for (int x : {mesh->xstart - 1, mesh->xend}) {
        const bool inner = x == mesh->xstart - 1;
        if (mesh->periodicX || (inner && !mesh->firstX()) || (!inner && !mesh->lastX())) {
          continue;
        }
        for (int y = mesh->ystart; y <= mesh->yend; ++y) {
          for (int z = 0; z < mesh->LocalNz; ++z) {
            wall_momentum_weight(x, y, z) =
                AA * wallFraction(N, AA, x, y, z) * midpoint(V, x, y, z) / Z;
          }
        }
      }
      // The local particle partition carries its local per-particle momentum.
      const Field3D momentum_source =
          -wallWeighted(coefficient * NV, &wall_momentum_weight);
      add(species["momentum_source"], momentum_source);
      set(species["polarisation_momentum_source"], momentum_source);
    }
  });
  // Phi + sum_s[(3/2 p_s + K_parallel_s) A_s/(Z_s rho)].
  DivEpol = divPol(E_adv);
  balance.add("ion_energy", "radial_polarisation_inner", 0.0);
  balance.add("ion_energy", "radial_polarisation_outer", 0.0);
  for (int x : {mesh->xstart - 1, mesh->xend}) {
    const bool inner = x == mesh->xstart - 1;
    if (mesh->periodicX || (inner && !mesh->firstX()) || (!inner && !mesh->lastX())) {
      continue;
    }
    const int cell = inner ? x + 1 : x;
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      const auto nz = static_cast<size_t>(mesh->LocalNz);
      std::vector<BoutReal> thermal(nz, 0.0), kinetic_parallel(nz, 0.0);
      std::vector<BoutReal> pressure_transport(nz, 0.0);
      forAllIons(allspecies, [&](GuardedOptions species) {
        const BoutReal AA = get<BoutReal>(species["AA"]);
        const BoutReal Z = get<BoutReal>(species["charge"]);
        const Field3D N = GET_NOBOUNDARY(Field3D, species["density"]);
        const Field3D T = GET_NOBOUNDARY(Field3D, species["temperature"]);
        const bool has_velocity = species.isSet("velocity");
        Field3D V;
        if (has_velocity) {
          V = GET_NOBOUNDARY(Field3D, species["velocity"]);
        }
        for (int z = 0; z < mesh->LocalNz; ++z) {
          const auto index = static_cast<size_t>(z);
          const BoutReal partition = wallFraction(N, AA, x, y, z);
          const BoutReal temp = midpoint(T, x, y, z);
          const BoutReal velocity = has_velocity ? midpoint(V, x, y, z) : 0.0;
          pressure_transport[index] += partition * temp / Z;
          thermal[index] += 1.5 * partition * temp / Z;
          kinetic_parallel[index] += 0.5 * AA * partition * velocity * velocity / Z;
        }
      });
      for (int z = 0; z < mesh->LocalNz; ++z) {
        const auto index = static_cast<size_t>(z);
        const BoutReal work =
            0.5 * (phi(x, y, z) + phi(x + 1, y, z) + psi(x, y, z) + psi(x + 1, y, z));
        const BoutReal pressure_mid = 0.5 * (psi(x, y, z) + psi(x + 1, y, z));
        const BoutReal transported = thermal[index] + kinetic_parallel[index];
        const BoutReal revised = work + transported;
        balance.add("ion_energy",
                    inner ? "radial_polarisation_inner" : "radial_polarisation_outer",
                    (inner ? -1.0 : 1.0) * polarisation_flux.x(x, y, z) * revised);
        const BoutReal volume = coord->J()(cell, y, z) * coord->dx()(cell, y, z)
                                * coord->dy()(cell, y, z) * coord->dz()(cell, y, z);
        DivEpol(cell, y, z) += (inner ? -1.0 : 1.0) * polarisation_flux.x(x, y, z)
                               * (revised - E_adv(cell, y, z)) / volume;
        wall_energy_weight(x, y, z) = revised;
        wall_work_weight(x, y, z) = work;
        wall_thermal_weight(x, y, z) = transported;
        wall_pressure_mismatch(x, y, z) = pressure_transport[index] - pressure_mid;
      }
    }
  }
  forAllIons(allspecies, [&](GuardedOptions species) {
    auto& sources = diagnostics[species.name()];
    sources["wall_energy_weight"] = wall_energy_weight;
    sources["wall_work_weight"] = wall_work_weight;
    sources["wall_thermal_weight"] = wall_thermal_weight;
    sources["wall_pressure_mismatch"] = wall_pressure_mismatch;
  });
  set(state["fields"]["ion_polarisation_energy_divergence"], DivEpol);
  set(state["fields"]["DivJpol"], DivJpol);
  set(state["fields"]["ion_chi"], chi);
  set(fields["DivJpar"], DivJpar);
  set(fields["ion_current_rhs"], current_rhs);
  set(fields["ion_auxiliary_residual"], auxiliary_residual);
}

void PolarisationDrift::outputConservative(Options& state) {
  balance.output(state, "polarisation_drift", diagnose);
  const auto Nnorm = get<BoutReal>(state["Nnorm"]);
  const auto Tnorm = get<BoutReal>(state["Tnorm"]);
  const auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  const auto rho_s0 = get<BoutReal>(state["rho_s0"]);
  const BoutReal Pnorm = SI::qe * Tnorm * Nnorm;
  set_with_attrs(state["polarisation_drift_version"], 1,
                 {{"description", "Conservative polarisation transport with local "
                                  "Phi wall work and local species transport weights"}});
  set_with_attrs(state["ion_chi"], chi,
                 {{"time_dimension", "t"},
                  {"units", "V s^-1"},
                  {"conversion", Tnorm * Omega_ci},
                  {"long_name", "Auxiliary polarisation current potential"},
                  {"description", "Approximates d(Phi)/dt; DC wall value is d(phi_w)/dt, "
                                  "the local pressure-potential wall rate is neglected, "
                                  "and AC defaults to Neumann"}});
  if (diagnose) {
    set_with_attrs(state["parallel_current"], Jpar,
                   {{"time_dimension", "t"},
                    {"units", "A m^-2"},
                    {"conversion", SI::qe * Nnorm * rho_s0 * Omega_ci},
                    {"source", "polarisation_drift"},
                    {"long_name", "Completed parallel current with sheath guards"}});
    for (const auto& [name, field] :
         {std::pair{"DivJpar", &DivJpar}, std::pair{"ion_current_rhs", &current_rhs},
          std::pair{"ion_auxiliary_residual", &auxiliary_residual}}) {
      set_with_attrs(state[name], *field,
                     {{"time_dimension", "t"},
                      {"units", "A m^-3"},
                      {"conversion", SI::qe * Nnorm * Omega_ci},
                      {"source", "polarisation_drift"}});
    }
    for (const auto& [species, sources] : diagnostics.getChildren()) {
      for (const auto& [name, description] :
           {std::pair{"wall_energy_weight",
                      "Physical X-wall face weight: local Phi plus local "
                      "thermal and parallel kinetic energy"},
            std::pair{"wall_work_weight", "Physical X-wall face work weight: local Phi "
                                          "midpoint, including pressure work"},
            std::pair{"wall_thermal_weight",
                      "Physical X-wall face weight: local thermal and "
                      "parallel kinetic energy"},
            std::pair{"wall_pressure_mismatch",
                      "Local species pressure transport weight minus the "
                      "local pressure-potential midpoint"}}) {
        set_with_attrs(state[std::string{"polarisation_"} + name],
                       sources[name].as<Field3D>(),
                       {{"time_dimension", "t"},
                        {"units", "V"},
                        {"conversion", Tnorm},
                        {"description", description}});
      }
      set_with_attrs(state["particle_pol_wall_loss_" + species],
                     sources["wall_loss"].as<Field3D>(),
                     {{"time_dimension", "t"},
                      {"units", "m^-3 s^-1"},
                      {"conversion", Nnorm * Omega_ci},
                      {"long_name", "Outward polarisation particle loss per volume"},
                      {"species", species}});
    }
    for (const auto& [name, faces, indexing] :
         {std::tuple{"polarisation_flux_x", &polarisation_flux.x,
                     "x(i,j,k): face (i+1/2,j,k)"},
          std::tuple{"polarisation_flux_z", &polarisation_flux.z,
                     "z(i,j,k): face (i,j,k+1/2)"}}) {
      set_with_attrs(
          state[name], *faces,
          {{"time_dimension", "t"},
           {"units", "A"},
           {"conversion", SI::qe * Nnorm * Omega_ci * pow(rho_s0, 3)},
           {"long_name", "Area-integrated X-Z polarisation current"},
           {"description", "Increasing-coordinate current: minus Reynolds and chi faces; "
                           "X-Z contribution, excluding optional poloidal X-Y transport"},
           {"face_index", indexing}});
    }
    set_with_attrs(state["DivJpol"], DivJpol,
                   {{"time_dimension", "t"},
                    {"units", "A m^-3"},
                    {"long_name", "Polarisation transport current divergence"},
                    {"conversion", SI::qe * Nnorm * Omega_ci}});
    set_with_attrs(
        state["ddt(Phi)"], chi,
        {{"time_dimension", "t"},
         {"units", "V s^-1"},
         {"conversion", Tnorm * Omega_ci},
         {"long_name", "Auxiliary polarisation current potential"},
         {"description", "Fixed-correction auxiliary solve; viscous stress "
                         "and numerical phi damping are excluded from its RHS; "
                         "nonperiodic X: midpoint Dirichlet DC, "
                         "radial AC selected by polarisation_ac_boundary; "
                         "approximate, not exact ddt(phi+psi)"}});
    set_with_attrs(state["ion_chi_iterations"], chi_iterations,
                   {{"time_dimension", "t"},
                    {"description", "Auxiliary chi solver calls: one initial solve plus "
                                    "polarisation_corrections defect corrections"}});
    set_with_attrs(state["ion_chi_residual"], chi_residual,
                   {{"time_dimension", "t"},
                    {"description", "Maximum absolute normalised finite-volume residual "
                                    "of L(chi)=DivJReynolds+DivJdia+DivJpar"
                                    "+DivJextra+DivJcol over physical cells; "
                                    "viscous stress excluded"}});
    set_with_attrs(
        state["DivEpol"], DivEpol,
        {{"time_dimension", "t"}, {"units", "W m^-3"}, {"conversion", Pnorm * Omega_ci}});
  }
  for (const auto* field :
       {"ion_chi", "ddt(Phi)", "ion_chi_iterations", "ion_chi_residual", "DivJpol",
        "DivEpol", "polarisation_flux_x", "polarisation_flux_z"}) {
    if (state.isSet(field)) {
      state[field].attributes["source"] = std::string{"polarisation_drift"};
    }
  }
}
