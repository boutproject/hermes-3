
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
#include <utility>

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

// Each mixed term comes from the same corner bilinear form. In particular,
// averaging transverse derivatives independently at faces would not give a
// symmetric operator on a nonuniform mesh.
bout::ion_energy::PerpendicularFaceFluxes
geometricPerpendicularFaceFluxes(const Field3D& a, const Field3D& f, const Ind3D& i) {
  const auto* coord = f.getCoordinates();
  const auto& J = coord->J();
  const auto& g11 = coord->g11();
  const auto& g33 = coord->g33();
  const auto& g13 = coord->g13();
  const auto& dx = coord->dx();
  const auto& dz = coord->dz();

  const auto xflux = [&](const Ind3D& l) {
    const auto r = l.xp();
    return 0.5 * (a[l] + a[r]) * (J[l] * g11[l] * dz[l] + J[r] * g11[r] * dz[r])
           * (f[r] - f[l]) / (dx[l] + dx[r]);
  };
  const auto zflux = [&](const Ind3D& l) {
    const auto r = l.zp();
    return 0.5 * (a[l] + a[r]) * (J[l] * g33[l] * dx[l] + J[r] * g33[r] * dx[r])
           * (f[r] - f[l]) / (dz[l] + dz[r]);
  };
  const auto corner = [&](const Ind3D& l) {
    const auto x = l.xp();
    const auto z = l.zp();
    const auto xz = x.zp();
    const BoutReal coefficient =
        0.0625 * (a[l] + a[x] + a[z] + a[xz])
        * (J[l] * g13[l] + J[x] * g13[x] + J[z] * g13[z] + J[xz] * g13[xz]);
    // Centre-to-centre lengths cancel against the corner area in the
    // integrated cross flux. These are averaged differences, not derivatives.
    return std::pair<BoutReal, BoutReal>{0.5 * coefficient * (f[x] - f[l] + f[xz] - f[z]),
                                         0.5 * coefficient
                                             * (f[z] - f[l] + f[xz] - f[x])};
  };
  const auto pp = corner(i);
  const auto mp = corner(i.xm());
  const auto pm = corner(i.zm());
  const auto mm = corner(i.xm().zm());
  return {xflux(i) + 0.5 * (pp.second + pm.second),
          xflux(i.xm()) + 0.5 * (mp.second + mm.second),
          zflux(i) + 0.5 * (pp.first + mp.first),
          zflux(i.zm()) + 0.5 * (pm.first + mm.first)};
}
} // namespace

namespace bout::ion_energy {
PerpendicularFaceFluxes perpendicularFaceFluxes(const Field3D& a, const Field3D& f,
                                                const Ind3D& i, bool bndry_flux) {
  auto flux = geometricPerpendicularFaceFluxes(a, f, i);
  if (!bndry_flux && !f.getMesh()->periodicX) {
    if (f.getMesh()->firstX() && i.x() == f.getMesh()->xstart) {
      flux.xm = 0.0;
    }
    if (f.getMesh()->lastX() && i.x() == f.getMesh()->xend) {
      flux.xp = 0.0;
    }
  }
  return flux;
}

Field3D Div_a_Grad_perp_XZ(const Field3D& a, const Field3D& f, bool bndry_flux) {
  ASSERT0(a.getLocation() == f.getLocation());
  Field3D result{zeroFrom(f)};
  const auto* coord = f.getCoordinates();
  BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
    const auto flux = perpendicularFaceFluxes(a, f, i, bndry_flux);
    // A zero X derivative alone does not eliminate the g13 normal flux.
    result[i] = (flux.xp - flux.xm + flux.zp - flux.zm)
                / (coord->J()[i] * coord->dx()[i] * coord->dz()[i]);
  }
  return result;
}

Field3D kineticEnergyPerp(const Field3D& a, const Field3D& f) {
  ASSERT0(a.getLocation() == f.getLocation());
  Field3D result{zeroFrom(f)};
  const auto* coord = f.getCoordinates();
  BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
    const auto flux = perpendicularFaceFluxes(a, f, i);
    // Half of each face's 1/2 * delta(f) * flux belongs to this cell.
    result[i] = ((f[i.xp()] - f[i]) * flux.xp + (f[i] - f[i.xm()]) * flux.xm
                 + (f[i.zp()] - f[i]) * flux.zp + (f[i] - f[i.zm()]) * flux.zm)
                / (4. * coord->J()[i] * coord->dx()[i] * coord->dz()[i]);
  }
  return result;
}
} // namespace bout::ion_energy

using bout::ion_energy::Div_a_Grad_perp_XZ;

EvolveVorticityIonEnergy::EvolveVorticityIonEnergy(std::string name, Options& alloptions,
                                                   Solver* solver)
    : NamedComponent(name, {readOnly("species:{all_species}:charge"),
                            readIfSet("species:{charged}:velocity", Regions::Interior),
                            readIfSet("species:{charged}:momentum", Regions::Interior),

                            // Electron inputs
                            readOnly("species:e:pressure", Regions::Interior),
                            readIfSet("species:e:AA"),
                            // Electron outputs
                            readWrite("species:e:energy_source"),
                            // Ion inputs
                            readOnly("species:{ions}:density", Regions::Interior),
                            readOnly("species:{ions}:AA"),
                            // Ion outputs
                            readWrite("species:{ions}:pressure"),
                            readWrite("species:{ions}:temperature"),
                            readWrite("species:{ions}:density_source"),
                            readWrite("species:{ions}:momentum_source"),
                            // Fields
                            readWrite("fields:vorticity"), readWrite("fields:phi")}) {

  auto& options = alloptions[name];

  phi_boundary_relax = options["phi_boundary_relax"]
                           .doc("Relax x boundaries of phi towards Neumann?")
                           .withDefault<bool>(false);
  if (phi_boundary_relax) {
    if (mesh->xstart < 2) {
      throw BoutException("phi_boundary_relax requires two radial guard cells");
    }
    auto& evolve_boundary = Options::root()["Vort"]["evolve_bndry"];
    if (evolve_boundary.isSet() && !evolve_boundary.as<bool>()) {
      throw BoutException("phi_boundary_relax requires Vort:evolve_bndry=true");
    }
    // Set an effective value, rather than a different default from Solver::add.
    evolve_boundary.force(true, "evolve_vorticity_ion_energy");
    phi_boundary_timescale = options["phi_boundary_timescale"]
                                 .doc("Timescale for phi boundary relaxation [seconds]")
                                 .withDefault(1e-4)
                             / get<BoutReal>(alloptions["units"]["seconds"]);
    if (!std::isfinite(phi_boundary_timescale) || phi_boundary_timescale <= 0.0) {
      throw BoutException("phi_boundary_timescale must be finite and positive");
    }
    // Own a no-op boundary configuration before moving the original Vort
    // boundary operators into the physical-vorticity field. Swap only the
    // FieldData metadata: solver pointers to Vort and ddt(Vort) stay valid.
    const std::string storage_name = name + "_boundary_storage";
    Options::root()[storage_name]["bndry_all"].force(std::string{"none"});
    vorticity.setBoundary(storage_name);
  }

  // Evolve the vorticity and ion energy in time
  solver->add(Vort, "Vort");
  solver->add(Eion, "Eion");
  if (phi_boundary_relax) {
    swap(static_cast<FieldData&>(Vort), static_cast<FieldData&>(vorticity));
    ddt(Vort).copyBoundary(Vort);
    initialiseBoundaryState();
  }

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

  if (phi_boundary_relax) {
    phi_core_averagey = options["phi_core_averagey"]
                            .doc("Average phi core boundary in Y?")
                            .withDefault<bool>(false)
                        and mesh->periodicY(mesh->xstart);

    phiSolver->setInnerBoundaryFlags(INVERT_SET);
    phiSolver->setOuterBoundaryFlags(INVERT_SET);
  }

  // Fixed point iteration settings
  phi_max_iterations = options["phi_max_iterations"]
                           .doc("Maximum number of finite-volume potential iterations")
                           .withDefault(100);
  phi_atol =
      options["phi_atol"].doc("phi convergence absolute tolerance [eV]").withDefault(1.0)
      / Tnorm;
  phi_rtol =
      options["phi_rtol"].doc("phi convergence relative tolerance").withDefault(1e-3);
  vort_atol = options["vort_atol"]
                  .doc("Absolute tolerance on the normalised forward vorticity residual")
                  .withDefault(1e-10);

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
  vorticity = copy(Vort);
  if (phi_boundary_relax) {
    // Assignment preserves vorticity's independently owned boundary rules.
    vorticity.applyBoundary();
    mesh->communicate(vorticity);
    applyPotentialBoundary();
  }

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
  Field3D mass_density = 0.0;
  forAllIons(allspecies, [&](GuardedOptions species) {
    const BoutReal AA = get<BoutReal>(species["AA"]);
    const BoutReal Z = get<BoutReal>(species["charge"]);
    auto N = GET_NOBOUNDARY(Field3D, species["density"]);
    Ni_sum += N;
    mass_charge_density += (AA / Z) * N;
    mass_density += AA * N;

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
  if (phi_boundary_relax) {
    // Start from this stage's energy rather than a previous pressure iterate.
    // Pressure guards must be consistent before the first Phi correction.
    Pi_sum = floor((2. / 3) * (Eion - ion_KE_parallel), 0.0);
    applyPressureBoundary(Pi_sum);
  }
  psi = pressure_weight * Pi_sum;
  const BoutReal vort_tolerance = vort_atol + phi_rtol * max(abs(Vort), true);

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

    mesh->communicate(phi, psi);
    applyPotentialBoundary();
    // The spectral solver is only an approximate inverse. Defect correction
    // converges to the same finite-volume operator used for kinetic energy.
    phi = correctPotential(Vort, phi + psi) - psi;
    applyPotentialBoundary();

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
    mesh->communicate(phi, psi);
    applyPotentialBoundary();

    phi_err = max((abs(phi - phi_last)) / (abs(phi) + phi_atol), true);
    const BoutReal psi_err = max(abs(psi - psi_last) / (abs(psi) + phi_atol), true);
    phi_err = std::max(phi_err, psi_err);
    const BoutReal vort_err =
        max(abs(Vort - Div_a_Grad_perp_XZ(boussinesq_coeff, phi + psi)), true);
    ++iteration;

    //output.write("Iteration {} phi_err {}\n", iteration, phi_err);

    if (phi_err <= phi_rtol && vort_err <= vort_tolerance) {
      break;
    }
    if (iteration >= phi_max_iterations) {
      throw BoutException("Potential recovery did not converge in {} iterations: "
                          "potential error {}, vorticity residual {} (tolerance {})",
                          iteration, phi_err, vort_err, vort_tolerance);
    }
  } while (true);
  mesh->communicate(phi, Pi_sum, psi);
  applyPotentialBoundary();

  // Free parallel boundaries, as in the vorticity component. The spectral
  // inverse does not supply these, but electron compression uses Grad(phi).
  if (mesh->ystart > 0) {
    Field3D phi_fa = toFieldAligned(phi);
    const int lower_neighbour = std::min(mesh->ystart + 1, mesh->yend);
    const int upper_neighbour = std::max(mesh->yend - 1, mesh->ystart);
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int z = 0; z < mesh->LocalNz; ++z) {
        phi_fa(r.ind, mesh->ystart - 1, z) =
            2 * phi_fa(r.ind, mesh->ystart, z) - phi_fa(r.ind, lower_neighbour, z);
      }
    }
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int z = 0; z < mesh->LocalNz; ++z) {
        phi_fa(r.ind, mesh->yend + 1, z) =
            2 * phi_fa(r.ind, mesh->yend, z) - phi_fa(r.ind, upper_neighbour, z);
      }
    }
    phi = fromFieldAligned(phi_fa);
  }
  applyPotentialBoundary();

  auto fields = state["fields"];
  set(fields["vorticity"], vorticity);
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
  // Use the same X-Z face operator as potential recovery and polarisation.
  DivJReynolds = -Div_n_bxGrad_f_B_XPPM(0.5 * vorticity, phi, false, poloidal_flows)
                 - Div_a_Grad_perp_XZ(Field3D{0.5 * boussinesq_coeff}, vEdotGradPsi)
                 - Div_n_bxGrad_f_B_XPPM(DelpPhi_2B2, phi + psi, false, poloidal_flows);

  // Divergence of Reynolds stress energy flux
  DivEReynolds = Div_a_Grad_perp_XZ(0.5 * boussinesq_coeff * vEdotGradPsi, phi + psi);

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
    const Field3D NV = GET_NOBOUNDARY(Field3D, electrons["momentum"]);
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
    const Field3D NV = GET_NOBOUNDARY(Field3D, species["momentum"]);
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

  // The total electron energy exchange is phi * Div(J_e). evolve_pressure
  // evolves thermal energy, so only ExB compression belongs in energy_source.
  // Parallel electrical work is already in the electron momentum equation.
  phi_Div_Je = phi * (DivJdia_e + Div_par(Jpar_e));
  electron_exb_work = -(Jdia_e * Grad(phi));
  add(electrons["energy_source"], electron_exb_work);

  //////////////////////////////////////////////////////////
  // Polarisation drift

  // Auxiliary current-continuity potential. Its homogeneous boundary is
  // independent of the relaxing phi boundary and does not require ddt(psi).
  const Field3D dVortdt = DivJReynolds + DivJdia;
  const BoutReal dVort_tolerance = vort_atol + phi_rtol * max(abs(dVortdt), true);
  for (int step = 0; step < phi_max_iterations; ++step) {
    mesh->communicate(chi);
    chi = correctPotential(dVortdt, chi);
    mesh->communicate(chi);
    const BoutReal residual =
        max(abs(dVortdt - Div_a_Grad_perp_XZ(boussinesq_coeff, chi)), true);
    if (residual <= dVort_tolerance) {
      break;
    }
    if (step + 1 == phi_max_iterations) {
      throw BoutException("Polarisation potential did not converge in {} iterations: "
                          "vorticity residual {} (tolerance {})",
                          phi_max_iterations, residual, dVort_tolerance);
    }
  }

  // Set Neumann boundaries; divPol also suppresses the g13 radial flux.
  // Solving for dPhidt with two Neumann boundaries is ill conditioned
  if (!phi_boundary_relax) {
    chi.applyBoundary("neumann_o2");
  }
  mesh->communicate(chi);

  // Closed polarisation transport, distinct from the full current solve.
  const auto divPol = [&](const Field3D& weight) -> Field3D {
    return
        // Potential flow part of polarisation drift
        -Div_a_Grad_perp_XZ(weight * boussinesq_coeff, chi, false)
        // Reynolds stress part of the polarisation drift
        // Note: The stream functions are negated so that upwinding works correctly
        + Div_n_bxGrad_f_B_XPPM(0.5 * vorticity, -phi, false, poloidal_flows, false,
                                &weight)
        - Div_a_Grad_perp_XZ(weight * 0.5 * boussinesq_coeff, vEdotGradPsi, false)
        + Div_n_bxGrad_f_B_XPPM(DelpPhi_2B2, -(phi + psi), false, poloidal_flows, false,
                                &weight);
  };
  DivJpol = divPol(Field3D{1.0});

  // Jps = (A_s n_s / rho) Jpol, Gamma_ps = Jps / Z_s. Use the actual
  // mass density, so the species currents sum to Jpol above the floor.
  mass_density = floor(mass_density, average_atomic_mass * density_floor);
  Field3D E_adv = copy(phi);
  polarisation_sources = Options{};
  forAllIons(allspecies, [&](GuardedOptions species) {
    const BoutReal AA = get<BoutReal>(species["AA"]);
    const BoutReal Z = get<BoutReal>(species["charge"]);
    const Field3D N = GET_NOBOUNDARY(Field3D, species["density"]);
    const Field3D P = GET_NOBOUNDARY(Field3D, species["pressure"]);
    const Field3D coefficient = (AA / Z) / mass_density;
    const Field3D density_source = -divPol(coefficient * N);
    add(species["density_source"], density_source);
    polarisation_sources[species.name()]["density"] = density_source;
    E_adv += coefficient * (5. / 2) * P;
    if (species.isSet("velocity")) {
      const Field3D V = GET_NOBOUNDARY(Field3D, species["velocity"]);
      E_adv += coefficient * 0.5 * AA * N * SQ(V);
    }
    if (species.isSet("momentum")) {
      const Field3D NV = GET_NOBOUNDARY(Field3D, species["momentum"]);
      const Field3D momentum_source = -divPol(coefficient * NV);
      add(species["momentum_source"], momentum_source);
      polarisation_sources[species.name()]["momentum"] = momentum_source;
    }
  });
  // Sum_s[(5/2 p_s + K_parallel_s) A_s/(Z_s rho)] + phi.
  DivEpol = divPol(E_adv);
}

void EvolveVorticityIonEnergy::finally(const Options& state) {
  //////////////////////////////////////////////////////////
  // Vorticity equation

  if (phi_boundary_relax) {
    ddt(Vort) = 0.0;
    BOUT_FOR(i, Vort.getRegion("RGN_NOBNDRY")) {
      ddt(Vort)[i] = DivJReynolds[i] + DivJdia[i];
    }
    boundaryRHS();
  } else {
    ddt(Vort) = DivJReynolds + DivJdia;
  }

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

  if (phi_boundary_relax) {
    set(state["Vort"], vorticity);
  }
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
    set_with_attrs(state["DivJpol"], DivJpol,
                   {{"time_dimension", "t"},
                    {"units", "A m^-3"},
                    {"long_name", "Closed polarisation transport current divergence"},
                    {"conversion", SI::qe * Nnorm * Omega_ci}});
    set_with_attrs(state["ion_pressure_potential"], psi,
                   {{"time_dimension", "t"}, {"units", "V"}, {"conversion", Tnorm}});
    set_with_attrs(state["ddt(Phi)"], chi,
                   {{"time_dimension", "t"},
                    {"units", "V s^-1"},
                    {"conversion", Tnorm * Omega_ci},
                    {"long_name", "Auxiliary polarisation current potential"},
                    {"description",
                     "L(chi)=physical ddt(Vort); approximate, not exact ddt(phi+psi)"}});
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
        state["electron_exb_work"], electron_exb_work,
        {{"time_dimension", "t"}, {"units", "W m^-3"}, {"conversion", Pnorm * Omega_ci}});
    set_with_attrs(
        state["Eion_src"], Eion_src,
        {{"time_dimension", "t"}, {"units", "W m^-3"}, {"conversion", Pnorm * Omega_ci}});
  }
}

void EvolveVorticityIonEnergy::initialiseBoundaryState() {
  BOUT_FOR(i, Vort.getRegion("RGN_BNDRY")) { Vort[i] = 0.0; }
  if (mesh->periodicX) {
    return;
  }
  const auto initialise = [&](int interior, int guard, int begin, int end) {
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      BoutReal midpoint = 0.0;
      for (int z = 0; z < mesh->LocalNz; ++z) {
        midpoint += 0.5 * (phi(guard, y, z) + phi(interior, y, z));
      }
      midpoint /= mesh->LocalNz;
      for (int x = begin; x < end; ++x) {
        for (int z = 0; z < mesh->LocalNz; ++z) {
          Vort(x, y, z) = midpoint;
        }
      }
    }
  };
  if (mesh->firstX()) {
    initialise(mesh->xstart, mesh->xstart - 1, 0, mesh->xstart);
  }
  if (mesh->lastX()) {
    initialise(mesh->xend, mesh->xend + 1, mesh->xend + 1, mesh->LocalNx);
  }
}

void EvolveVorticityIonEnergy::applyPotentialBoundary() {
  if (!phi_boundary_relax || mesh->periodicX) {
    return;
  }
  const auto apply = [&](int interior, int guard, int begin, int end) {
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      BoutReal midpoint = 0.0;
      for (int z = 0; z < mesh->LocalNz; ++z) {
        midpoint += Vort(guard, y, z);
      }
      midpoint /= mesh->LocalNz;
      for (int z = 0; z < mesh->LocalNz; ++z) {
        phi(guard, y, z) = 2.0 * midpoint - phi(interior, y, z);
        for (int x = begin; x < end; ++x) {
          phi(x, y, z) = phi(guard, y, z);
        }
      }
    }
  };
  if (mesh->firstX()) {
    apply(mesh->xstart, mesh->xstart - 1, 0, mesh->xstart - 1);
  }
  if (mesh->lastX()) {
    apply(mesh->xend, mesh->xend + 1, mesh->xend + 2, mesh->LocalNx);
  }
}

void EvolveVorticityIonEnergy::boundaryRHS() {
  if (mesh->periodicX) {
    return;
  }
  BoutReal core_average = 0.0;
  if (mesh->firstX() && phi_core_averagey) {
    BoutReal local[2]{
        0.0, static_cast<BoutReal>((mesh->yend - mesh->ystart + 1) * mesh->LocalNz)};
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      for (int z = 0; z < mesh->LocalNz; ++z) {
        local[0] += phi(mesh->xstart, y, z);
      }
    }
    BoutReal global[2];
    MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, mesh->getYcomm(0));
    core_average = global[0] / global[1];
  }
  const auto relax = [&](int interior, int begin, int end, bool average_y) {
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      BoutReal target = core_average;
      if (!average_y) {
        target = 0.0;
        for (int z = 0; z < mesh->LocalNz; ++z) {
          target += phi(interior, y, z);
        }
        target /= mesh->LocalNz;
      }
      for (int x = begin; x < end; ++x) {
        for (int z = 0; z < mesh->LocalNz; ++z) {
          ddt(Vort)(x, y, z) = (target - Vort(x, y, z)) / phi_boundary_timescale;
        }
      }
    }
  };
  if (mesh->firstX()) {
    relax(mesh->xstart, 0, mesh->xstart, phi_core_averagey);
  }
  if (mesh->lastX()) {
    relax(mesh->xend, mesh->xend + 1, mesh->LocalNx, false);
  }
}

void EvolveVorticityIonEnergy::restartVars(Options& state) {
  // The first call is the postInit hook, after solver variables were loaded.
  if (restore_phi_on_restart) {
    restore_phi_on_restart = false;
    if (state.isSet("phi")) {
      phi = state["phi"].as<Field3D>();
    }
    if (phi_boundary_relax) {
      const int version = state["ion_energy_boundary_state_version"].withDefault(0);
      if (version < 0 || version > 1) {
        throw BoutException("Unsupported ion energy boundary state version {}", version);
      }
      if (version == 0) {
        initialiseBoundaryState();
      }
    }
  }
  if (phi_boundary_relax) {
    set_with_attrs(state["ion_energy_boundary_state_version"], 1,
                   {{"description", "Restart Vort physical radial boundary cells store "
                                    "phi midpoints in normalised potential units"}});
  }
  set_with_attrs(state["phi"], phi,
                 {{"long_name", "plasma potential"}, {"source", "vorticity"}});
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
      src += get<Field3D>(species["energy_source"]); // Thermal energy only
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
      // Polarisation advection is already present in the total-energy flux.
      // Remove its local kinetic-energy conversion using the completed V.
      if (polarisation_sources.isSection(species.name())) {
        const auto& pol = polarisation_sources[species.name()];
        if (pol.isSet("momentum")) {
          src -= V * get<Field3D>(pol["momentum"]);
        }
        if (pol.isSet("density")) {
          const BoutReal AA = get<BoutReal>(species["AA"]);
          src += 0.5 * AA * SQ(V) * get<Field3D>(pol["density"]);
        }
      }
    }
  });
  return src;
}

// Use exactly the same face fluxes as the forward vorticity operator.
Field3D EvolveVorticityIonEnergy::calculateKEperp(const Field3D& phi,
                                                  const Field3D& psi) const {
  return bout::ion_energy::kineticEnergyPerp(boussinesq_coeff, phi + psi);
}

Field3D EvolveVorticityIonEnergy::correctPotential(const Field3D& rhs,
                                                   const Field3D& potential) const {
  // Homogeneous correction boundaries preserve the separately imposed phi
  // boundaries, including INVERT_SET values supplied by a restart.
  Field3D residual{zeroFrom(rhs)};
  const Field3D forward = Div_a_Grad_perp_XZ(boussinesq_coeff, potential);
  BOUT_FOR(i, residual.getRegion("RGN_NOBNDRY")) {
    residual[i] = (rhs[i] - forward[i]) / boussinesq_coeff[i];
  }
  const Field3D approximate = phiSolver->solve(residual, zeroFrom(potential));
  // Spectral solvers can leave Y guard cells uninitialised. Only copy the
  // X-Z planes they solve; communication supplies inter-process Y guards.
  Field3D correction{zeroFrom(potential)};
  BOUT_FOR(i, correction.getRegion("RGN_NOY")) { correction[i] = approximate[i]; }
  mesh->communicate(correction);
  return potential + correction;
}
