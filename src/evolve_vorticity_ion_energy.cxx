
#include <bout/assert.hxx>
#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/constants.hxx>
#include <bout/coordinates.hxx>
#include <bout/derivs.hxx>
#include <bout/difops.hxx>
#include <bout/field.hxx>
#include <bout/field3d.hxx>
#include <bout/fv_ops.hxx>
#include <bout/globals.hxx>
#include <bout/initialprofiles.hxx>
#include <bout/invert_pardiv.hxx>
#include <bout/options.hxx>
#include <bout/output.hxx>
#include <bout/output_bout_types.hxx>
#include <bout/region.hxx>
#include <bout/solver.hxx>
#include <bout/vecops.hxx>

#include "../include/component.hxx"
#include "../include/div_ops.hxx"
#include "../include/evolve_vorticity_ion_energy.hxx"
#include "../include/guarded_options.hxx"
#include "../include/hermes_build_config.hxx"
#include "../include/hermes_utils.hxx"
#include "../include/permissions.hxx"

#include <algorithm>
#include <memory>
#include <string>

using bout::globals::mesh;

namespace {
/// Call a function with each ion species
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
template <typename Function>
void forAllIons(const Options& allspecies, Function func) {
  for (const auto& kv : allspecies.getChildren()) {
    if (kv.first == "e") {
      continue; // Skip electrons
    }
    const auto& species = allspecies[kv.first];

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

Field3D Div_a_Grad_perp_XZ(const Field3D& a, const Field3D& f) {
  ASSERT0(a.getLocation() == f.getLocation());

  Field3D result{zeroFrom(f)};

  const Coordinates* coord = f.getCoordinates();

  const auto& J = coord->J();
  const auto& g11 = coord->g11();
  const auto& g33 = coord->g33();
  const auto& dx = coord->dx();
  const auto& dz = coord->dz();

  BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
    const BoutReal volume = J[i] * dx[i] * dz[i];

    // Flux in x
    const auto ixp = i.xp();
    const auto ixm = i.xm();

    // Area of cell edges
    const BoutReal area_xp = 0.5 * (J[i] * g11[i] * dz[i] + J[ixp] * g11[ixp] * dz[ixp]);
    const BoutReal area_xm = 0.5 * (J[i] * g11[i] * dz[i] + J[ixm] * g11[ixm] * dz[ixm]);

    // Coefficient `a` at cell edges
    const BoutReal coeff_xp = 0.5 * (a[i] + a[ixp]);
    const BoutReal coeff_xm = 0.5 * (a[i] + a[ixm]);

    // Gradient at cell edges
    const BoutReal grad_xp = 2. * (f[ixp] - f[i]) / (dx[i] + dx[ixp]);
    const BoutReal grad_xm = 2. * (f[i] - f[ixm]) / (dx[i] + dx[ixm]);

    // Flux in z
    const auto izp = i.zp();
    const auto izm = i.zm();

    const BoutReal area_zp = 0.5 * (J[i] * g33[i] * dx[i] + J[izp] * g33[izp] * dx[izp]);
    const BoutReal area_zm = 0.5 * (J[i] * g33[i] * dx[i] + J[izm] * g33[izm] * dx[izm]);

    const BoutReal coeff_zp = 0.5 * (a[i] + a[izp]);
    const BoutReal coeff_zm = 0.5 * (a[i] + a[izm]);

    const BoutReal grad_zp = 2. * (f[izp] - f[i]) / (dz[i] + dz[izp]);
    const BoutReal grad_zm = 2. * (f[i] - f[izm]) / (dz[i] + dz[izm]);

    result[i] = (+(area_xp * coeff_xp * grad_xp) - (area_xm * coeff_xm * grad_xm)
                 + (area_zp * coeff_zp * grad_zp) - (area_zm * coeff_zm * grad_zm))
                / volume;
  }

  return result;
}

}; // namespace

EvolveVorticityIonEnergy::EvolveVorticityIonEnergy(std::string name, Options& alloptions,
                                                   Solver* solver)
    : NamedComponent(name,
                     {readOnly("species:{all_species}:charge"),
                      readIfSet("species:{charged}:velocity"),
                      readIfSet("species:{charged}:momentum"),

                      // Electron inputs
                      readOnly("species:e:pressure"), readIfSet("species:e:AA"),
                      // Electron outputs
                      readWrite("species:e:energy_source"),
                      // Ion inputs
                      readOnly("species:{ions}:density"), readOnly("species:{ions}:AA"),
                      // Ion outputs
                      readWrite("species:{ions}:pressure"),
                      readWrite("species:{ions}:temperature"),
                      // Fields
                      readWrite("fields:vorticity"), readWrite("fields:phi")}) {

  auto& options = alloptions[name];

  // Evolve the vorticity and ion energy in time
  solver->add(Vort, "Vort");
  solver->add(Eion, "Eion");

  const Options& units = alloptions["units"];
  const BoutReal Nnorm = units["inv_meters_cubed"];
  const BoutReal Tnorm = units["eV"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();
  const BoutReal Bnorm = units["Tesla"];
  const BoutReal Lnorm = units["meters"];

  const auto* coord = mesh->getCoordinates();

  average_atomic_mass = options["average_atomic_mass"]
                            .doc("Weighted average atomic mass, for polarisation current "
                                 "(Boussinesq approximation)")
                            .withDefault<BoutReal>(2.0); // Deuterium

  poloidal_flows =
      options["poloidal_flows"].doc("Include poloidal ExB flow").withDefault<bool>(false);

  neumann_boundary_average_z = options["neumann_boundary_average_z"]
                                   .doc("Apply neumann boundary with Z average?")
                                   .withDefault<bool>(false);

  density_floor = options["density_floor"].doc("Minimum density floor").withDefault(1e-7);

  // Add phi to restart files so that the value in the boundaries
  // is restored on restart. This is done even when phi is not evolving,
  // so that phi can be saved and re-loaded

  // Set initial value. Will be overwritten if restarting
  Pi_sum = (2. / 3) * Eion;

  // Coefficients for Boussinesq solve
  boussinesq_coeff = average_atomic_mass / SQ(coord->Bxy());

  phiSolver = Laplacian::create(&options["laplacian"]);
  phiSolver->setCoefC(boussinesq_coeff);

  phi_boundary_relax = options["phi_boundary_relax"]
                           .doc("Relax x boundaries of phi towards Neumann?")
                           .withDefault<bool>(false);

  if (phi_boundary_relax) {
    // Set the last update time to -1, so it will reset
    // the first time RHS function is called
    phi_boundary_last_update = -1.;

    phi_core_averagey = options["phi_core_averagey"]
                            .doc("Average phi core boundary in Y?")
                            .withDefault<bool>(false)
                        and mesh->periodicY(mesh->xstart);

    phi_boundary_timescale = options["phi_boundary_timescale"]
                                 .doc("Timescale for phi boundary relaxation [seconds]")
                                 .withDefault(1e-4)
                             / get<BoutReal>(alloptions["units"]["seconds"]);
    // Normalise to internal time units

    phiSolver->setInnerBoundaryFlags(INVERT_SET);
    phiSolver->setOuterBoundaryFlags(INVERT_SET);
  }

  // Fixed point iteration settings
  phi_max_iterations = options["phi_max_iterations"]
                           .doc("Maximum number of phi <-> Pi iterations")
                           .withDefault(5);
  phi_atol =
      options["phi_atol"].doc("phi convergence absolute tolerance [eV]").withDefault(1.0)
      / Tnorm;
  phi_rtol =
      options["phi_rtol"].doc("phi convergence relative tolerance").withDefault(1e-3);

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
      output_warn.write("No curvature vector in input grid");
      Curlb_B = 0.0;
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

  Curlb_B *= 2. / coord->Bxy();

  diagnose =
      options["diagnose"].doc("Output additional diagnostics?").withDefault<bool>(false);
}

void EvolveVorticityIonEnergy::transform_impl(GuardedOptions& state) {
  mesh->communicate(Vort, Eion);

  auto allspecies = state["species"];
  auto electrons = allspecies["e"];

  // From vorticity and ion energy E we will recover
  // the potential phi (ϕ) and ion pressure p_i
  //
  // With Boussinesq approximation this can be solved once for Φ but
  // applying boundary conditions on pressure and phi separately
  // requires iteration.

  // Calculate K_||i = ( m_i n v_||i )^2 / (2 m_i n)

  ion_KE_parallel = 0.0;
  Field3D Ni_sum = 0.0;
  Field3D mass_charge_density = 0.0;
  forAllIons(allspecies, [&](GuardedOptions species) {
    const BoutReal AA = get<BoutReal>(species["AA"]);
    const BoutReal Z = get<BoutReal>(species["charge"]);
    auto N = GET_NOBOUNDARY(Field3D, species["density"]);
    Ni_sum += N;
    mass_charge_density += (AA / Z) * N;

    if (!species.isSet("velocity")) {
      return;
    }
    auto V = GET_NOBOUNDARY(Field3D, species["velocity"]);

    ion_KE_parallel += 0.5 * AA * N * SQ(V);
  });

  const Field3D Ni_sum_floor = softFloor(Ni_sum, density_floor);
  // With a common temperature, psi = Ti * sum_s(A_s n_s / Z_s) / Abar.
  // Use the same density floor as the temperature recovery.
  const Field3D pressure_weight =
      mass_charge_density / (average_atomic_mass * Ni_sum_floor);
  psi = pressure_weight * Pi_sum;

  int iteration = 0;
  BoutReal phi_err;
  do {
    // Both potentials must converge: phi can be unchanged on the first
    // iteration even though the recovered pressure potential has changed.
    const Field3D phi_last = copy(phi);
    const Field3D psi_last = copy(psi);

    // Solve
    // ∇·(c∇⊥Φ) = ω for Φ = ϕ + psi
    // with c = rho_0 / B^2 and psi = sum_s(A_s p_s / (Abar Z_s))

    phi = phiSolver->solve(
        (Vort - Div_a_Grad_perp_XZ(boussinesq_coeff, psi)) / boussinesq_coeff, phi);

    // With the Boussinesq approximation the perpendicular flow velocity is
    // v⊥0 = b × ∇Φ / B = v_E + v∗
    // The mass density is constant so kinetic energy K⊥ = 1/2 c |∇⊥Φ|^2,

    ion_KE_perp = calculateKEperp(phi, psi);

    // Calculate pressure
    // E = Cv * P + (1/2) m n v^2 + K_\perp
    // Ensure that P does not go negative.
    Pi_sum = floor((2. / 3) * (Eion - ion_KE_parallel - ion_KE_perp), 0.0);

    applyPressureBoundary(Pi_sum);
    psi = pressure_weight * Pi_sum;

    phi_err = max((abs(phi - phi_last)) / (abs(phi) + phi_atol), true);
    const BoutReal psi_err = max(abs(psi - psi_last) / (abs(psi) + phi_atol), true);
    phi_err = std::max(phi_err, psi_err);
    ++iteration;

    //output.write("Iteration {} phi_err {}\n", iteration, phi_err);

    if (iteration > phi_max_iterations) {
      output_warn.write("\nMaximum phi iterations {} reached. phi error {}\n", iteration,
                        phi_err);
      break;
    }
  } while (phi_err > phi_rtol);
  mesh->communicate(phi, Pi_sum, psi);

  auto fields = state["fields"];
  set(fields["vorticity"], Vort);
  set(fields["phi"], phi);

  //////////////////////////////////////////////////////////
  // Calculate common ion temperature
  // Set temperature and density of each ion species

  Ti = Pi_sum / Ni_sum_floor;
  forAllIons(allspecies, [&](GuardedOptions species) {
    set(species["temperature"], Ti);
    auto Ni = GET_NOBOUNDARY(Field3D, species["density"]);
    set(species["pressure"], Ni * Ti);
  });

  //////////////////////////////////////////////////////////
  // Reynolds stress

  // V_ExB dot Grad(psi)
  Field3D vEdotGradPsi = bracket(phi, psi, BRACKET_ARAKAWA);
  vEdotGradPsi.applyBoundary("neumann_o2");

  // delp2(phi) term
  Field3D DelpPhi_2B2 = 0.5 * boussinesq_coeff * Delp2(phi);
  DelpPhi_2B2.applyBoundary("neumann_o2");

  mesh->communicate(vEdotGradPsi, DelpPhi_2B2);

  // Divergence of Reynolds stress term DivJReynolds = -Div(G)
  DivJReynolds = -Div_n_bxGrad_f_B_XPPM(0.5 * Vort, phi, false, poloidal_flows)
                 - FV::Div_a_Grad_perp(Field3D{0.5 * boussinesq_coeff}, vEdotGradPsi)
                 - Div_n_bxGrad_f_B_XPPM(DelpPhi_2B2, phi + psi, false, poloidal_flows);

  // Divergence of Reynolds stress enegy flux
  DivEReynolds = FV::Div_a_Grad_perp(0.5 * boussinesq_coeff * vEdotGradPsi, phi + psi);

  //////////////////////////////////////////////////////////
  // Diamagnetic drift

  // Electron diamagnetic current
  const auto Pe = GET_NOBOUNDARY(Field3D, electrons["pressure"]);
  const Vector3D Jdia_e = Pe * Curlb_B;

  // Ion diamagnetic current
  const Vector3D Jdia_i = Pi_sum * Curlb_B;

  DivJdia = Div(Jdia_e + Jdia_i);
  const Field3D DivJdia_e = Div(Jdia_e);

  //////////////////////////////////////////////////////////
  // Parallel current

  Field3D Jpar_e; // Electron current
  if (electrons.isSet("momentum")) {
    // Note: Using NV rather than N*V so that the cell boundary flux is correct
    const Field3D NV = get<Field3D>(electrons["momentum"]);
    const BoutReal A = get<BoutReal>(electrons["AA"]);
    Jpar_e = (-1. / A) * NV;
  } else {
    Jpar_e = 0.0;
  }

  Field3D Jpar_i = 0.0; // Ion current
  forAllIons(allspecies, [&](GuardedOptions species) {
    if (!species.isSet("momentum")) {
      return;
    }
    const Field3D NV = get<Field3D>(species["momentum"]);
    const BoutReal A = get<BoutReal>(species["AA"]);
    const BoutReal Z = get<BoutReal>(species["charge"]);

    Jpar_i += (Z / A) * NV;
  });

  DivEelectric = Div(phi * Jdia_i) + Div_par(phi * Jpar_i);

  //////////////////////////////////////////////////////////
  // ExB drift

  // ExB drift of thermal and parallel ion kinetic energy
  DivEExB = Div_n_bxGrad_f_B_XPPM((3. / 2) * Pi_sum + ion_KE_parallel, phi, false,
                                  poloidal_flows, true);

  // Energy transfer to/from electrons: phi * Div(J_e)
  // Add to electrons, subtract from ions
  phi_Div_Je = phi * (DivJdia_e + Div_par(Jpar_e));

  add(electrons["energy_source"], phi_Div_Je);

  //////////////////////////////////////////////////////////
  // Polarisation drift

  // Note: This will be handled in the `polarisation_drift`
  // component, because it needs to include other currents
  // due to e.g. ion viscosity

  // Energy density advected by polarisation drift
  // pi + 1/2 m_i n_0 v^2_||i + en_0ϕ
  // Note the n_0 in the parallel kinetic energy and phi term
  Field3D E_adv = (5. / 2) * Pi_sum + phi;

  forAllIons(allspecies, [&](GuardedOptions species) {
    if (!species.isSet("velocity")) {
      return;
    }
    const BoutReal AA = get<BoutReal>(species["AA"]);
    auto V = GET_NOBOUNDARY(Field3D, species["velocity"]);

    // Density is n_0 due to Boussinesq approximation
    E_adv += 0.5 * AA * SQ(V);
  });

  // Calculate Div( E_adv v_pol )
  //
  // vpol = - c / (en_0) Grad_perp dPhi/dt - G / (en_0)

  // Solve for Div( c Grad_perp dPhi/dt ) = dVort/dt
  // Note: Check case where boundaries are relaxed
  const Field3D dVortdt = DivJReynolds + DivJdia;
  dPhidt = phiSolver->solve(dVortdt / boussinesq_coeff, dPhidt);

  // Set Neumann boundaries so energy fluxes through boundaries are disabled.
  // Solving for dPhidt with two Neumann boundaries is ill conditioned
  dPhidt.applyBoundary("neumann_o2");

  // Divergence of energy flux due carried by polarisation drift
  DivEpol =
      // Potential flow part of polarisation drift
      -Div_a_Grad_perp_XZ(E_adv * boussinesq_coeff, dPhidt)
      // Reynolds stress part of the polarisation drift
      // Note: The stream functions are negated so that upwinding works correctly
      + Div_n_bxGrad_f_B_XPPM(E_adv * 0.5 * Vort, -phi, false, poloidal_flows)
      - Div_a_Grad_perp_XZ(E_adv * 0.5 * boussinesq_coeff, vEdotGradPsi)
      + Div_n_bxGrad_f_B_XPPM(E_adv * DelpPhi_2B2, -(phi + psi), false, poloidal_flows);

  // Fluxes of ion densities and (if set) momentum
  // This requires a partition of the polarisation current between species.
}

void EvolveVorticityIonEnergy::finally(const Options& state) {
  //////////////////////////////////////////////////////////
  // Vorticity equation

  ddt(Vort) = DivJReynolds + DivJdia;

  //////////////////////////////////////////////////////////
  // Ion energy equation

  // Each species carries its own parallel enthalpy and kinetic energy.
  // Read these from the completed state so boundary values are available.
  DivEparallel = 0.0;
  forAllIons(state["species"], [&](const Options& species) {
    if (!species.isSet("velocity")) {
      return;
    }
    const Field3D N = get<Field3D>(species["density"]);
    const Field3D V = get<Field3D>(species["velocity"]);
    const Field3D P = get<Field3D>(species["pressure"]);
    const BoutReal AA = get<BoutReal>(species["AA"]);
    const Field3D energy = (5. / 2) * P + 0.5 * AA * N * SQ(V);
    const Field3D fastest_wave = state.isSet("fastest_wave")
                                     ? get<Field3D>(state["fastest_wave"])
                                     : sqrt(get<Field3D>(species["temperature"]) / AA);
    Field3D flow_ylow;
    DivEparallel += FV::Div_par_mod<hermes::Limiter>(energy, V, fastest_wave, flow_ylow);
  });

  Eion_src = sumIonEnergySources(state["species"]);

  ddt(Eion) = -(DivEExB + DivEReynolds + DivEelectric + DivEparallel + DivEpol)
              - phi_Div_Je + Eion_src;
}

void EvolveVorticityIonEnergy::outputVars(Options& state) {
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto rho_s0 = get<BoutReal>(state["rho_s0"]);

  BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

  state["Eion"].setAttributes({{"time_dimension", "t"},
                               {"units", "J/m^3"},
                               {"conversion", Pnorm},
                               {"standard_name", "E_ion"},
                               {"long_name", "Ion energy density"},
                               {"species", "all ions"},
                               {"source", "evolve_vorticity_ion_energy"}});

  set_with_attrs(state["Tion"], Ti,
                 {{"time_dimension", "t"},
                  {"units", "eV"},
                  {"conversion", Tnorm},
                  {"standard_name", "temperature"},
                  {"long_name", "Common ion temperature"},
                  {"species", "all ions"},
                  {"source", "evolve_vorticity_ion_energy"}});

  state["Vort"].setAttributes({{"time_dimension", "t"},
                               {"units", "C m^-3"},
                               {"conversion", SI::qe * Nnorm},
                               {"long_name", "vorticity"},
                               {"source", "evolve_vorticity_ion_energy"}});

  set_with_attrs(state["phi"], phi,
                 {{"time_dimension", "t"},
                  {"units", "V"},
                  {"conversion", Tnorm},
                  {"standard_name", "potential"},
                  {"long_name", "plasma potential"},
                  {"source", "evolve_vorticity_ion_energy"}});

  if (diagnose) {
    set_with_attrs(state["ion_pressure_potential"], psi,
                   {{"time_dimension", "t"}, {"units", "V"}, {"conversion", Tnorm}});
    set_with_attrs(
        state["ddt(Phi)"], dPhidt,
        {{"time_dimension", "t"}, {"units", "V s^-1"}, {"conversion", Tnorm * Omega_ci}});
    set_with_attrs(state["ion_KE_parallel"], ion_KE_parallel,
                   {{"time_dimension", "t"}, {"units", "J m^-3"}, {"conversion", Pnorm}});
    set_with_attrs(state["ion_KE_perp"], ion_KE_perp,
                   {{"time_dimension", "t"}, {"units", "J m^-3"}, {"conversion", Pnorm}});
    set_with_attrs(
        state["ddt(Eion)"], ddt(Eion),
        {{"time_dimension", "t"}, {"units", "W m^-3"}, {"conversion", Pnorm * Omega_ci}});
    set_with_attrs(
        state["DivEExB"], DivEExB,
        {{"time_dimension", "t"}, {"units", "W m^-3"}, {"conversion", Pnorm * Omega_ci}});
    set_with_attrs(
        state["DivEReynolds"], DivEReynolds,
        {{"time_dimension", "t"}, {"units", "W m^-3"}, {"conversion", Pnorm * Omega_ci}});
    set_with_attrs(
        state["DivEelectric"], DivEelectric,
        {{"time_dimension", "t"}, {"units", "W m^-3"}, {"conversion", Pnorm * Omega_ci}});
    set_with_attrs(
        state["DivEparallel"], DivEparallel,
        {{"time_dimension", "t"}, {"units", "W m^-3"}, {"conversion", Pnorm * Omega_ci}});
    set_with_attrs(
        state["DivEpol"], DivEpol,
        {{"time_dimension", "t"}, {"units", "W m^-3"}, {"conversion", Pnorm * Omega_ci}});
    set_with_attrs(
        state["phi_Div_Je"], phi_Div_Je,
        {{"time_dimension", "t"}, {"units", "W m^-3"}, {"conversion", Pnorm * Omega_ci}});
    set_with_attrs(
        state["Eion_src"], Eion_src,
        {{"time_dimension", "t"}, {"units", "W m^-3"}, {"conversion", Pnorm * Omega_ci}});
  }
}

void EvolveVorticityIonEnergy::applyPressureBoundary(Field3D& P) const {
  if (neumann_boundary_average_z) {
    // Take Z (usually toroidal) average and apply as X (radial) boundary condition
    if (mesh->firstX()) {
      for (int j = mesh->ystart; j <= mesh->yend; j++) {
        BoutReal Pavg = 0.0; // Average P in Z
        for (int k = 0; k < mesh->LocalNz; k++) {
          Pavg += P(mesh->xstart, j, k);
        }
        Pavg /= mesh->LocalNz;

        // Apply boundary condition
        for (int k = 0; k < mesh->LocalNz; k++) {
          P(mesh->xstart - 1, j, k) = (2. * Pavg) - P(mesh->xstart, j, k);
          P(mesh->xstart - 2, j, k) = P(mesh->xstart - 1, j, k);
        }
      }
    }

    if (mesh->lastX()) {
      for (int j = mesh->ystart; j <= mesh->yend; j++) {
        BoutReal Pavg = 0.0; // Average P in Z
        for (int k = 0; k < mesh->LocalNz; k++) {
          Pavg += P(mesh->xend, j, k);
        }
        Pavg /= mesh->LocalNz;

        for (int k = 0; k < mesh->LocalNz; k++) {
          P(mesh->xend + 1, j, k) = (2. * Pavg) - P(mesh->xend, j, k);
          P(mesh->xend + 2, j, k) = P(mesh->xend + 1, j, k);
        }
      }
    }
  } else {
    P.applyBoundary("neumann");
  }
}

Field3D EvolveVorticityIonEnergy::sumIonEnergySources(const Options& allspecies) const {
  Field3D src = 0.0;
  forAllIons(allspecies, [&](const Options& species) {
    if (species.isSet("energy_source")) {
      src += get<Field3D>(species["energy_source"]); // For diagnostic output
    }
    if (species.isSet("velocity")) {
      const Field3D V = get<Field3D>(species["velocity"]);
      if (species.isSet("momentum_source")) {
        src += V * get<Field3D>(species["momentum_source"]);
      }
      if (species.isSet("density_source")) {
        const BoutReal AA = get<BoutReal>(species["AA"]);
        src -= 0.5 * AA * SQ(V) * get<Field3D>(species["density_source"]);
      }
    }
  });
  return src;
}

// Perpendicular kinetic energy
// This uses the same discretisation at Div_a_Grad_perp_XZ
Field3D EvolveVorticityIonEnergy::calculateKEperp(const Field3D& phi,
                                                  const Field3D& psi) const {

  const Coordinates* coord = phi.getCoordinates();

  const auto& J = coord->J();
  const auto& g11 = coord->g11();
  const auto& g33 = coord->g33();
  const auto& dx = coord->dx();
  const auto& dz = coord->dz();

  const Field3D a = boussinesq_coeff;
  const Field3D f = phi + psi;

  Field3D result{zeroFrom(f)};

  BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
    const BoutReal volume = J[i] * dx[i] * dz[i];

    // X cell faces
    const auto ixp = i.xp();
    const auto ixm = i.xm();

    // Area of cell edges
    const BoutReal area_xp = 0.5 * (J[i] * g11[i] * dz[i] + J[ixp] * g11[ixp] * dz[ixp]);
    const BoutReal area_xm = 0.5 * (J[i] * g11[i] * dz[i] + J[ixm] * g11[ixm] * dz[ixm]);

    // Coefficient `a` at cell edges
    const BoutReal coeff_xp = 0.5 * (a[i] + a[ixp]);
    const BoutReal coeff_xm = 0.5 * (a[i] + a[ixm]);

    // Gradient at cell edges
    const BoutReal grad2_xp = 2. * SQ(f[ixp] - f[i]) / (dx[i] + dx[ixp]);
    const BoutReal grad2_xm = 2. * SQ(f[i] - f[ixm]) / (dx[i] + dx[ixm]);

    // Flux in z
    const auto izp = i.zp();
    const auto izm = i.zm();

    const BoutReal area_zp = 0.5 * (J[i] * g33[i] * dx[i] + J[izp] * g33[izp] * dx[izp]);
    const BoutReal area_zm = 0.5 * (J[i] * g33[i] * dx[i] + J[izm] * g33[izm] * dx[izm]);

    const BoutReal coeff_zp = 0.5 * (a[i] + a[izp]);
    const BoutReal coeff_zm = 0.5 * (a[i] + a[izm]);

    const BoutReal grad2_zp = 2. * SQ(f[izp] - f[i]) / (dz[i] + dz[izp]);
    const BoutReal grad2_zm = 2. * SQ(f[i] - f[izm]) / (dz[i] + dz[izm]);

    result[i] = (+(area_xp * coeff_xp * grad2_xp) + (area_xm * coeff_xm * grad2_xm)
                 + (area_zp * coeff_zp * grad2_zp) + (area_zm * coeff_zm * grad2_zm))
                / (4. * volume);
  }

  return result;
}
