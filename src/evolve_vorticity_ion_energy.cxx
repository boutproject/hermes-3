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
#include <bout/sys/range.hxx>
#include <bout/vecops.hxx>

#include "../include/component.hxx"
#include "../include/div_ops.hxx"
#include "../include/evolve_vorticity_ion_energy.hxx"
#include "../include/guarded_options.hxx"
#include "../include/hermes_build_config.hxx"
#include "../include/hermes_utils.hxx"
#include "../include/permissions.hxx"

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <utility>

using bout::globals::mesh;

namespace {
/// Update every physical radial guard column, preserving interior and Y guards.
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

/// Average a local X-Y line over Z for wall-state and geometry calculations.
///
/// @param field Read-only potential, storage field or 2D/3D metric field.
/// @param x Local radial cell index.
/// @param y Local poloidal cell index.
/// @returns Arithmetic mean over all LocalNz cells, in the input field's units.
///          A 2D field is treated as constant in Z. Uses the global mesh and
///          assumes a complete local Z line; no communication or mutation occurs.
template <typename Field>
BoutReal meanZ(const Field& field, int x, int y) {
  BoutReal mean = 0.0;
  for (int z = 0; z < mesh->LocalNz; ++z) {
    mean += field(x, y, z);
  }
  return mean / mesh->LocalNz;
}

/// Average per-Y wall values over the connected flux surface with J dy weights.
///
/// @param values Read-only wall values indexed by local Y, with LocalNy entries.
///               Only physical Y rows contribute to the average.
/// @param interior Local index of the first interior radial cell at the wall.
/// @param guard Local index of its adjacent radial guard cell.
/// @returns Weighted surface mean in the input values' units, identical on all
///          ranks in mesh->getYcomm(interior). The face weight is the mean of
///          guard/interior J times interior dy, using Z-averaged geometry.
///
/// This MPI collective must be called by every rank in that Y communicator.
/// Callers restrict it to a first-X core surface and assume positive weights.
BoutReal fluxSurfaceAverage(const std::vector<BoutReal>& values, int interior,
                            int guard) {
  const auto* coord = mesh->getCoordinates();
  BoutReal local[2]{0.0, 0.0};
  for (int y = mesh->ystart; y <= mesh->yend; ++y) {
    const BoutReal weight =
        0.5 * (meanZ(coord->J(), guard, y) + meanZ(coord->J(), interior, y))
        * meanZ(coord->dy(), interior, y);
    local[0] += weight * values[static_cast<size_t>(y)];
    local[1] += weight;
  }
  BoutReal global[2];
  MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, mesh->getYcomm(interior));
  return global[0] / global[1];
}

/// Broadcast per-Y wall values into all physical radial guard columns and Z cells.
///
/// @param[in,out] field Destination storage, RHS or diagnostic field. Only radial
///                     guards at physical Y rows on this rank's X walls are written;
///                     interior cells and Y boundary cells retain their values.
/// @param inner Read-only inner wall values in field's units, indexed by local Y.
///              Must contain LocalNy entries when this rank owns the inner wall.
/// @param outer Read-only outer wall values, with the same convention as inner.
///
/// Uses the global mesh without communication. Periodic X leaves field unchanged.
void writeRadialValues(Field3D& field, const std::vector<BoutReal>& inner,
                       const std::vector<BoutReal>& outer) {
  if (mesh->periodicX) {
    return;
  }
  const auto write = [&](int begin, int end, const std::vector<BoutReal>& values) {
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      for (int x = begin; x < end; ++x) {
        for (int z = 0; z < mesh->LocalNz; ++z) {
          field(x, y, z) = values[static_cast<size_t>(y)];
        }
      }
    }
  };
  if (mesh->firstX()) {
    write(0, mesh->xstart, inner);
  }
  if (mesh->lastX()) {
    write(mesh->xend + 1, mesh->LocalNx, outer);
  }
}

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
    return 0.5 * (a[l] + a[r]) * ((J[l] * g11[l] * dz[l]) + (J[r] * g11[r] * dz[r]))
           * (f[r] - f[l]) / (dx[l] + dx[r]);
  };
  const auto zflux = [&](const Ind3D& l) {
    const auto r = l.zp();
    return 0.5 * (a[l] + a[r]) * ((J[l] * g33[l] * dx[l]) + (J[r] * g33[r] * dx[r]))
           * (f[r] - f[l]) / (dz[l] + dz[r]);
  };
  const auto corner = [&](const Ind3D& l) {
    const auto x = l.xp();
    const auto z = l.zp();
    const auto xz = x.zp();
    const BoutReal coefficient =
        0.0625 * (a[l] + a[x] + a[z] + a[xz])
        * ((J[l] * g13[l]) + (J[x] * g13[x]) + (J[z] * g13[z]) + (J[xz] * g13[xz]));
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
void applyWallPressureBoundary(Field3D& f) {
  auto* localmesh = f.getMesh();
  if (localmesh->periodicX) {
    return;
  }
  if (localmesh->getNZPE() != 1 || localmesh->zstart != 0
      || localmesh->zend + 1 != localmesh->LocalNz) {
    throw BoutException("radial_boundary requires complete periodic Z lines");
  }
  // Calculate the mean once for each wall row, then copy the first guard outward.
  const Field2D mean = DC(f);
  fillRadialGuards(
      f, [&](int x, int, int, int y, int z) { return 2.0 * mean(x, y) - f(x, y, z); });
}

void applyVorticityWallBoundary(Field3D& w, bool midpoint_zero) {
  fillRadialGuards(w, [&](int x, int, int, int y, int z) {
    return midpoint_zero ? -w(x, y, z) : w(x, y, z);
  });
}

BoutReal wallInwardDerivative(BoutReal phi_midpoint, BoutReal phi0, BoutReal dx0) {
  // The stored midpoint is the arithmetic mean of a mirrored guard/interior
  // pair, so their centred face derivative is exact for quadratic profiles.
  return 2.0 * (phi0 - phi_midpoint) / dx0;
}

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
    const auto flux = perpendicularFaceFluxes(a, f, i, true);
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

  const Options& units = alloptions["units"];
  const BoutReal Tnorm = units["eV"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();
  const BoutReal Bnorm = units["Tesla"];
  const BoutReal Lnorm = units["meters"];

  using bout::ion_energy::PolarisationACBoundary;
  using bout::ion_energy::RadialBoundary;
  radial_model.walls = !mesh->periodicX;
  radial_model.mode = options["radial_boundary"]
                          .doc("Radial wall model: fixed, free_slip or no_slip")
                          .withDefault(RadialBoundary::fixed);
  radial_model.chi_ac =
      options["polarisation_ac_boundary"]
          .doc("Auxiliary potential AC wall condition: neumann or dirichlet")
          .withDefault(PolarisationACBoundary::neumann);
  radial_model.evolve_wall_potential = radial_model.mode == RadialBoundary::no_slip;
  radial_model.vorticity_wall_zero = radial_model.mode == RadialBoundary::free_slip;
  radial_model.viscous_wall_flux = radial_model.vorticity_wall_zero;
  radial_model.chi_wall_flux = radial_model.chi_ac == PolarisationACBoundary::dirichlet;
  radial_model.chi_flags = radial_model.chi_wall_flux ? 0 : INVERT_AC_GRAD;
  radial_model.core_flux_surface_average =
      radial_model.walls && mesh->firstX() && mesh->periodicY(mesh->xstart);

  // These flags are owned by the wall model, including when X is periodic.
  // Reject before creating a solver, which would otherwise consume the setting.
  for (const auto* section : {"laplacian", "polarisation_laplacian"}) {
    for (const auto* key : {"inner_boundary_flags", "outer_boundary_flags"}) {
      if (options[section][key].isSet()) {
        throw BoutException("{}:{} is controlled by radial_boundary", section, key);
      }
    }
  }
  auto& vort_options = Options::root()["Vort"];
  for (const auto* key : {"bndry_xin", "bndry_xout"}) {
    if (vort_options[key].isSet()) {
      throw BoutException("Vort:{} is controlled by radial_boundary", key);
    }
  }
  auto& evolve_boundary = vort_options["evolve_bndry"];
  if (evolve_boundary.isSet() && !evolve_boundary.as<bool>()) {
    throw BoutException("radial_boundary requires Vort:evolve_bndry=true");
  }
  if (radial_model.evolve_wall_potential) {
    const BoutReal kappa_SI =
        options["phi_boundary_velocity"]
            .doc("Velocity relaxing the mean wall-normal E-field to zero [m/s]. "
                 "For boundary-cell physical width h and an old relaxation "
                 "timescale tau, kappa = h/(2*tau) matches at that resolution.")
            .as<BoutReal>();
    if (!std::isfinite(kappa_SI) || kappa_SI <= 0.0) {
      throw BoutException("radial_boundary=no_slip requires finite positive "
                          "phi_boundary_velocity");
    }
    radial_model.phi_boundary_velocity = kappa_SI / (Lnorm * Omega_ci);
  }
  if (radial_model.walls) {
    if (mesh->xstart < 2) {
      throw BoutException("radial_boundary requires two radial guard cells");
    }
    evolve_boundary.force(true, "evolve_vorticity_ion_energy");
    wall_rate_inner.resize(static_cast<size_t>(mesh->LocalNy), 0.0);
    wall_rate_outer.resize(static_cast<size_t>(mesh->LocalNy), 0.0);
    // Own a no-op boundary configuration before moving the original Vort
    // boundary operators into the physical-vorticity field. Swap only the
    // FieldData metadata: solver pointers to Vort and ddt(Vort) stay valid.
    const std::string storage_name = name + "_boundary_storage";
    Options::root()[storage_name]["bndry_all"].force(std::string{"none"});
    vorticity.setBoundary(storage_name);
  }
  output_info.write(
      "{}: radial_boundary={}, walls={}, polarisation_ac_boundary={}, "
      "evolve_wall_potential={}, vorticity_wall_zero={}, viscous_wall_flux={}, "
      "chi_wall_flux={}, core_flux_surface_average={}, "
      "phi_boundary_velocity(normalised)={}, phi_flags={}, chi_flags={}\n",
      name, toString(radial_model.mode), radial_model.walls,
      toString(radial_model.chi_ac), radial_model.evolve_wall_potential,
      radial_model.vorticity_wall_zero, radial_model.viscous_wall_flux,
      radial_model.chi_wall_flux, radial_model.core_flux_surface_average,
      radial_model.phi_boundary_velocity, radial_model.phi_flags, radial_model.chi_flags);

  // Evolve the vorticity and ion energy in time
  solver->add(Vort, "Vort");
  solver->add(Eion, "Eion");

  if (radial_model.walls) {
    swap(static_cast<FieldData&>(Vort), static_cast<FieldData&>(vorticity));
    ddt(Vort).copyBoundary(Vort);
    initialiseBoundaryState();
  }

  const auto* coord = mesh->getCoordinates();

  average_atomic_mass = options["average_atomic_mass"]
                            .doc("Weighted average atomic mass, for polarisation current "
                                 "(Boussinesq approximation)")
                            .withDefault<BoutReal>(2.0); // Deuterium

  poloidal_flows =
      options["poloidal_flows"].doc("Include poloidal ExB flow").withDefault<bool>(false);

  density_floor = options["density_floor"].doc("Minimum density floor").withDefault(1e-7);

  // Match vorticity's physical units, scalar coefficient and guards.
  include_viscosity = options.isSet("viscosity");
  viscosity =
      options["viscosity"].doc("Kinematic viscosity [m^2/s]").withDefault<BoutReal>(0.0)
      / (Lnorm * Lnorm * Omega_ci);
  viscosity.applyBoundary("dirichlet");
  // Keep the physical coefficient at radial faces. Dirichlet coefficient guards
  // would cancel nu at the wall and suppress even the free-slip viscous flux.
  if (radial_model.walls) {
    fillRadialGuards(viscosity,
                     [&](int x, int, int, int y, int z) { return viscosity(x, y, z); });
  }
  mesh->communicate(viscosity);
  if (include_viscosity) {
    setPermissions(readIfSet("species:e:density", Regions::Interior));
  }

  // Add phi to restart files so that the value in the boundaries
  // is restored on restart. This is done even when phi is not evolving,
  // so that phi can be saved and re-loaded

  // Coefficients for Boussinesq solve
  boussinesq_coeff = average_atomic_mass / SQ(coord->Bxy());

  phiSolver = Laplacian::create(&options["laplacian"]);
  phiSolver->setCoefC(boussinesq_coeff);

  if (radial_model.walls) {
    auto& chi_options = options["polarisation_laplacian"];
    if (!chi_options["type"].isSet()) {
      chi_options["type"] = options["laplacian"]["type"].as<std::string>();
    }
    const int chi_global_flags = chi_options["global_flags"].withDefault(0);
    if (chi_global_flags & (INVERT_ZERO_DC | INVERT_KX_ZERO)) {
      throw BoutException("polarisation_laplacian must retain the DC solution");
    }
    chiSolver = Laplacian::create(&chi_options);
    chiSolver->setCoefC(boussinesq_coeff);
    chiSolver->setInnerBoundaryFlags(radial_model.chi_flags);
    chiSolver->setOuterBoundaryFlags(radial_model.chi_flags);
  }

  if (radial_model.walls) {
    phiSolver->setInnerBoundaryFlags(radial_model.phi_flags);
    phiSolver->setOuterBoundaryFlags(radial_model.phi_flags);
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

  polarisation_corrections =
      options["polarisation_corrections"]
          .doc("Finite-volume defect corrections after the spectral chi solve")
          .withDefault(2);
  if (polarisation_corrections < 0) {
    throw BoutException("polarisation_corrections must be >= 0");
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

void EvolveVorticityIonEnergy::acceptWarmStart() {
  if (!have_solution) {
    return; // No converged recovery yet: keep the initial or restart guess
  }
  phi_guess = copy(phi);
  ion_KE_perp_guess = copy(ion_KE_perp);
}

void EvolveVorticityIonEnergy::transform_impl(GuardedOptions& state) {
  mesh->communicate(Vort, Eion);
  vorticity = copy(Vort);

  // Start from the guess stored at the last output state, not from whatever
  // the previous RHS call (possibly a Jacobian probe or rejected stage) left.
  // This keeps the RHS a deterministic function of the state within a step.
  phi = copy(phi_guess);

  if (radial_model.walls) {
    // Assignment preserves vorticity's independently owned boundary rules.
    vorticity.applyBoundary();
    bout::ion_energy::applyVorticityWallBoundary(vorticity,
                                                 radial_model.vorticity_wall_zero);
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
  // Initial pressure: this state's E and K_|| with K⊥ from the last output
  // state. Equivalent to reusing the previous Pi_sum corrected for the change in
  // (E - K_||), so only the change in K⊥ is left for the iteration. Before the
  // first output the guess is K⊥ = 0.
  // Pressure guards must be consistent before the first Phi correction.
  Pi_sum = floor((2. / 3) * (Eion - ion_KE_parallel - ion_KE_perp_guess), 0.0);
  bout::ion_energy::applyWallPressureBoundary(Pi_sum);
  psi = pressure_weight * Pi_sum;
  bout::ion_energy::applyWallPressureBoundary(psi);
  const BoutReal vort_tolerance = vort_atol + (phi_rtol * max(abs(Vort), true));

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

    bout::ion_energy::applyWallPressureBoundary(Pi_sum);
    psi = pressure_weight * Pi_sum;
    bout::ion_energy::applyWallPressureBoundary(psi);
    mesh->communicate(phi, psi);
    applyPotentialBoundary();

    phi_err = max((abs(phi - phi_last)) / (abs(phi) + phi_atol), true);
    const BoutReal psi_err = max(abs(psi - psi_last) / (abs(psi) + phi_atol), true);
    phi_err = std::max(phi_err, psi_err);
    const BoutReal vort_err =
        max(abs(Vort - Div_a_Grad_perp_XZ(boussinesq_coeff, phi + psi, true)), true);
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
  phi_iterations = iteration;
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
            (2 * phi_fa(r.ind, mesh->ystart, z)) - phi_fa(r.ind, lower_neighbour, z);
      }
    }
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int z = 0; z < mesh->LocalNz; ++z) {
        phi_fa(r.ind, mesh->yend + 1, z) =
            (2 * phi_fa(r.ind, mesh->yend, z)) - phi_fa(r.ind, upper_neighbour, z);
      }
    }
    phi = fromFieldAligned(phi_fa);
  }
  applyPotentialBoundary();
  have_solution = true; // phi and ion_KE_perp are converged for this state
  if (radial_model.walls) {
    computeWallRates();
  }

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
  vEdotGradPsi.applyBoundary("neumann_o2"); // Retain the Y treatment
  bout::ion_energy::applyVorticityWallBoundary(vEdotGradPsi, true); // Q_wall = 0

  // delp2(phi) term
  Field3D DelpPhi_2B2 = 0.5 * boussinesq_coeff * Delp2(phi);
  DelpPhi_2B2.applyBoundary("neumann_o2"); // Retain the Y treatment
  if (radial_model.walls) {
    // free_o2 linear extrapolation, independent of mesh boundary labels.
    fillRadialGuards(DelpPhi_2B2, [&](int x, int step, int layer, int y, int z) {
      return DelpPhi_2B2(x, y, z)
             + layer * (DelpPhi_2B2(x, y, z) - DelpPhi_2B2(x - step, y, z));
    });
  }

  mesh->communicate(vEdotGradPsi, DelpPhi_2B2);

  // Divergence of Reynolds stress term DivJReynolds = -Div(G)
  // Use the same X-Z face operator as potential recovery and polarisation.
  DivJReynolds =
      -Div_n_bxGrad_f_B_XPPM(0.5 * vorticity, phi, false, poloidal_flows)
      - Div_a_Grad_perp_XZ(Field3D{0.5 * boussinesq_coeff}, vEdotGradPsi, false)
      - Div_n_bxGrad_f_B_XPPM(DelpPhi_2B2, phi + psi, false, poloidal_flows);

  // Divergence of Reynolds stress energy flux
  {
    const Field3D Phi = phi + psi;
    DivEReynolds = Div_a_Grad_perp_XZ(0.5 * boussinesq_coeff * vEdotGradPsi, Phi, false);
  }

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
  // Note: This requires the sheath boundary conditions

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
  // Kinematic viscosity

  if (include_viscosity) {
    calculateViscosity();
    // Match vorticity's mass-density partition, including electrons when
    // their mass and density are available. Ion heating is already internal
    // to Eion; only transfer the electron share out of the ion reservoir.
    if (electrons.isSet("AA") && electrons.isSet("density")) {
      const Field3D electron_mass_density =
          get<BoutReal>(electrons["AA"]) * GET_NOBOUNDARY(Field3D, electrons["density"]);
      electron_viscous_heating = viscous_heating * electron_mass_density
                                 / floor(mass_density + electron_mass_density,
                                         average_atomic_mass * density_floor);
      add(electrons["energy_source"], electron_viscous_heating);
    }
  }

  //////////////////////////////////////////////////////////
  // Polarisation drift

  // Auxiliary current-continuity potential, independent of phi relaxation.
  // DC Dirichlet fixes the constant mode; the AC enum controls the auxiliary
  // guards and whether its normal potential-flow current enters transport.
  const Field3D dVortdt = DivJReynolds + DivJdia + DivJviscous;
  solvePolarisationPotential(dVortdt);

  // Explicit wall closure for transport, distinct from the full current solve.
  const auto divPol = [&](const Field3D& weight) -> Field3D {
    return
        // Potential flow part of polarisation drift
        -Div_a_Grad_perp_XZ(weight * boussinesq_coeff, chi, radial_model.chi_wall_flux)
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
    if (diagnose) {
      const Field3D particle_coefficient = coefficient * N * boussinesq_coeff;
      // Outward particle loss per cell volume, supported only next to X walls.
      // Its volume integral is the signed particle flux through those walls.
      polarisation_sources[species.name()]["wall_loss"] =
          Div_a_Grad_perp_XZ(particle_coefficient, chi, false)
          - Div_a_Grad_perp_XZ(particle_coefficient, chi, radial_model.chi_wall_flux);
    }
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

  if (radial_model.walls) {
    ddt(Vort) = 0.0;
    BOUT_FOR(i, Vort.getRegion("RGN_NOBNDRY")) {
      ddt(Vort)[i] = DivJReynolds[i] + DivJdia[i] + DivJviscous[i];
    }
    writeWallRates();
  } else {
    ddt(Vort) = DivJReynolds + DivJdia + DivJviscous;
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

  ddt(Eion) =
      -(DivEExB + DivEReynolds + DivEelectric + DivEparallel + DivEpol + DivEviscous)
      - phi_Div_Je - electron_viscous_heating + Eion_src;
}

void EvolveVorticityIonEnergy::calculateViscosity() {
  const Field3D nu{viscosity};
  const Field3D Phi = phi + psi;
  const auto* coord = mesh->getCoordinates();
  DivJviscous = zeroFrom(vorticity);
  DivEviscous = zeroFrom(vorticity);
  viscous_heating = zeroFrom(vorticity);
  electron_viscous_heating = 0.0;
  BOUT_FOR(i, vorticity.getRegion("RGN_NOBNDRY")) {
    const auto flux = bout::ion_energy::perpendicularFaceFluxes(
        nu, vorticity, i, radial_model.viscous_wall_flux);
    const BoutReal inverse_volume =
        1. / (coord->J()[i] * coord->dx()[i] * coord->dz()[i]);
    DivJviscous[i] = (flux.xp - flux.xm + flux.zp - flux.zm) * inverse_volume;
    // Weight the common current at faces, rather than interpolating nu*Phi.
    // In particular a zero viscosity face must also have zero energy flux.
    DivEviscous[i] =
        0.5
        * ((Phi[i] + Phi[i.xp()]) * flux.xp - (Phi[i] + Phi[i.xm()]) * flux.xm
           + (Phi[i] + Phi[i.zp()]) * flux.zp - (Phi[i] + Phi[i.zm()]) * flux.zm)
        * inverse_volume;
    // Discrete product rule: Qnu = Phi*DivJviscous - DivEviscous.
    // Ion total energy carries the flux, not an additional Qnu source.
    viscous_heating[i] =
        -0.5
        * ((Phi[i.xp()] - Phi[i]) * flux.xp + (Phi[i] - Phi[i.xm()]) * flux.xm
           + (Phi[i.zp()] - Phi[i]) * flux.zp + (Phi[i] - Phi[i.zm()]) * flux.zm)
        * inverse_volume;
  }
}

void EvolveVorticityIonEnergy::outputVars(Options& state) {
  // Outputs follow an RHS evaluation at the accepted output state.
  // Update the recovery guess here so it remains fixed between outputs.
  acceptWarmStart();

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

  if (radial_model.walls) {
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
    set_with_attrs(state["DivJReynolds"], DivJReynolds,
                   {{"time_dimension", "t"},
                    {"units", "A m^-3"},
                    {"conversion", SI::qe * Nnorm * Omega_ci},
                    {"long_name", "Reynolds current divergence"}});
    for (const auto& [species, sources] : polarisation_sources.getChildren()) {
      set_with_attrs(state["particle_pol_wall_loss_" + species],
                     sources["wall_loss"].as<Field3D>(),
                     {{"time_dimension", "t"},
                      {"units", "m^-3 s^-1"},
                      {"conversion", Nnorm * Omega_ci},
                      {"long_name", "Outward polarisation particle loss per volume"},
                      {"species", species}});
    }
    set_with_attrs(
        state["ion_chi"], chi,
        {{"time_dimension", "t"},
         {"units", "V s^-1"},
         {"conversion", Tnorm * Omega_ci},
         {"long_name", "Auxiliary polarisation current potential"},
         {"description", "DC wall value is d(phi_b)/dt during relaxation; "
                         "the wall rate of the pressure potential is neglected"}});
    if (radial_model.walls) {
      Field3D wall_rate = zeroFrom(phi);
      writeRadialValues(wall_rate, wall_rate_inner, wall_rate_outer);
      set_with_attrs(state["phi_wall_rate"], wall_rate,
                     {{"time_dimension", "t"},
                      {"units", "V s^-1"},
                      {"conversion", Tnorm * Omega_ci},
                      {"long_name", "Radial midpoint potential rate"},
                      {"description", "Wall rates in physical radial guard cells; "
                                      "zero in interior and Y boundary cells"}});
    }
    if (include_viscosity) {
      set_with_attrs(state["kinematic_viscosity"], viscosity,
                     {{"units", "m^2 s^-1"},
                      {"conversion", SQ(rho_s0) * Omega_ci},
                      {"long_name", "Kinematic viscosity"}});
      set_with_attrs(state["kinematic_viscous_heating"], viscous_heating,
                     {{"time_dimension", "t"},
                      {"units", "W m^-3"},
                      {"conversion", Pnorm * Omega_ci},
                      {"long_name", "Discrete kinematic viscous heating"}});
      set_with_attrs(state["electron_viscous_heating"], electron_viscous_heating,
                     {{"time_dimension", "t"},
                      {"units", "W m^-3"},
                      {"conversion", Pnorm * Omega_ci},
                      {"long_name", "Electron share of kinematic viscous heating"}});
      set_with_attrs(state["DivJviscous"], DivJviscous,
                     {{"time_dimension", "t"},
                      {"units", "A m^-3"},
                      {"conversion", SI::qe * Nnorm * Omega_ci},
                      {"long_name", "Kinematic viscous current divergence"}});
      set_with_attrs(state["DivEviscous"], DivEviscous,
                     {{"time_dimension", "t"},
                      {"units", "W m^-3"},
                      {"conversion", Pnorm * Omega_ci},
                      {"long_name", "Kinematic viscous energy flux divergence"}});
    }
    set_with_attrs(state["DivJpol"], DivJpol,
                   {{"time_dimension", "t"},
                    {"units", "A m^-3"},
                    {"long_name", "Closed polarisation transport current divergence"},
                    {"conversion", SI::qe * Nnorm * Omega_ci}});
    set_with_attrs(state["ion_pressure_potential"], psi,
                   {{"time_dimension", "t"}, {"units", "V"}, {"conversion", Tnorm}});
    set_with_attrs(
        state["ddt(Phi)"], chi,
        {{"time_dimension", "t"},
         {"units", "V s^-1"},
         {"conversion", Tnorm * Omega_ci},
         {"long_name", "Auxiliary polarisation current potential"},
         {"description", "Fixed-correction approximation to L(chi)=physical ddt(Vort); "
                         "nonperiodic X: midpoint Dirichlet DC, "
                         "radial Neumann AC; "
                         "approximate, not exact ddt(phi+psi)"}});
    set_with_attrs(state["ion_phi_iterations"], phi_iterations,
                   {{"time_dimension", "t"},
                    {"description", "Iterations used in the last coupled pressure and "
                                    "potential recovery"}});
    set_with_attrs(state["ion_chi_iterations"], chi_iterations,
                   {{"time_dimension", "t"},
                    {"description", "Auxiliary chi solver calls: one initial solve plus "
                                    "polarisation_corrections defect corrections"}});
    set_with_attrs(state["ion_chi_residual"], chi_residual,
                   {{"time_dimension", "t"},
                    {"description", "Maximum absolute normalised finite-volume residual "
                                    "of L(chi)=physical ddt(Vort) over physical cells"}});
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
  if (!radial_model.walls) {
    return;
  }
  BOUT_FOR(i, Vort.getRegion("RGN_BNDRY")) { Vort[i] = 0.0; }
  const auto initialise = [&](int interior, int guard, int begin, int end,
                              bool average_y) {
    std::vector<BoutReal> midpoints(static_cast<size_t>(mesh->LocalNy), 0.0);
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      midpoints[static_cast<size_t>(y)] =
          0.5 * (meanZ(phi, guard, y) + meanZ(phi, interior, y));
    }
    if (average_y) {
      const BoutReal average = fluxSurfaceAverage(midpoints, interior, guard);
      std::fill(midpoints.begin(), midpoints.end(), average);
    }
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      for (int x = begin; x < end; ++x) {
        for (int z = 0; z < mesh->LocalNz; ++z) {
          Vort(x, y, z) = midpoints[static_cast<size_t>(y)];
        }
      }
    }
  };
  if (mesh->firstX()) {
    initialise(mesh->xstart, mesh->xstart - 1, 0, mesh->xstart,
               radial_model.core_flux_surface_average);
  }
  if (mesh->lastX()) {
    initialise(mesh->xend, mesh->xend + 1, mesh->xend + 1, mesh->LocalNx, false);
  }
}

void EvolveVorticityIonEnergy::applyPotentialBoundary() {
  if (!radial_model.walls) {
    return;
  }
  const auto apply = [&](int interior, int guard, int begin, int end, bool average_y) {
    std::vector<BoutReal> midpoints(static_cast<size_t>(mesh->LocalNy), 0.0);
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      midpoints[static_cast<size_t>(y)] = meanZ(Vort, guard, y);
    }
    if (average_y) {
      const BoutReal average = fluxSurfaceAverage(midpoints, interior, guard);
      std::fill(midpoints.begin(), midpoints.end(), average);
    }
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      for (int z = 0; z < mesh->LocalNz; ++z) {
        phi(guard, y, z) = 2.0 * midpoints[static_cast<size_t>(y)] - phi(interior, y, z);
        for (int x = begin; x < end; ++x) {
          phi(x, y, z) = phi(guard, y, z);
        }
      }
    }
  };
  if (mesh->firstX()) {
    apply(mesh->xstart, mesh->xstart - 1, 0, mesh->xstart - 1,
          radial_model.core_flux_surface_average);
  }
  if (mesh->lastX()) {
    apply(mesh->xend, mesh->xend + 1, mesh->xend + 2, mesh->LocalNx, false);
  }
}

void EvolveVorticityIonEnergy::computeWallRates() {
  // Reset every RHS evaluation; no rate survives from a previous solver state.
  std::fill(wall_rate_inner.begin(), wall_rate_inner.end(), 0.0);
  std::fill(wall_rate_outer.begin(), wall_rate_outer.end(), 0.0);
  if (!radial_model.walls || !radial_model.evolve_wall_potential) {
    return;
  }
  const auto* coord = mesh->getCoordinates();
  const auto compute = [&](int interior, int guard, bool average_y,
                           std::vector<BoutReal>& rates) {
    BoutReal core_midpoint = 0.0;
    if (average_y) {
      std::vector<BoutReal> midpoints(static_cast<size_t>(mesh->LocalNy), 0.0);
      for (int y = mesh->ystart; y <= mesh->yend; ++y) {
        midpoints[static_cast<size_t>(y)] = meanZ(Vort, guard, y);
      }
      core_midpoint = fluxSurfaceAverage(midpoints, interior, guard);
    }
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      const BoutReal midpoint = average_y ? core_midpoint : meanZ(Vort, guard, y);
      const BoutReal derivative = bout::ion_energy::wallInwardDerivative(
          midpoint, meanZ(phi, interior, y), meanZ(coord->dx(), interior, y));
      // For a 2D metric the mean g13 contribution vanishes on periodic Z.
      // With BOUT_USE_METRIC_3D we use Z-averaged geometry, averaging g11
      // before taking the square root to estimate the DC normal derivative.
      const BoutReal normal_metric =
          sqrt(0.5 * (meanZ(coord->g11(), guard, y) + meanZ(coord->g11(), interior, y)));
      rates[static_cast<size_t>(y)] =
          radial_model.phi_boundary_velocity * normal_metric * derivative;
    }
    if (average_y) {
      const BoutReal average = fluxSurfaceAverage(rates, interior, guard);
      // Y guard entries stay zero, including the DC values passed to chi.
      std::fill(rates.begin() + mesh->ystart, rates.begin() + mesh->yend + 1, average);
    }
  };
  if (mesh->firstX()) {
    compute(mesh->xstart, mesh->xstart - 1, radial_model.core_flux_surface_average,
            wall_rate_inner);
  }
  if (mesh->lastX()) {
    compute(mesh->xend, mesh->xend + 1, false, wall_rate_outer);
  }
}

void EvolveVorticityIonEnergy::writeWallRates() {
  writeRadialValues(ddt(Vort), wall_rate_inner, wall_rate_outer);
}

void EvolveVorticityIonEnergy::restartVars(Options& state) {
  // The first call is the postInit hook, after solver variables were loaded.
  if (restore_phi_on_restart) {
    restore_phi_on_restart = false;
    if (state.isSet("phi")) {
      phi = state["phi"].as<Field3D>();
      phi_guess = copy(phi);
    }
    if (radial_model.walls) {
      const int version = state["ion_energy_boundary_state_version"].withDefault(0);
      if (version < 0 || version > 1) {
        throw BoutException("Unsupported ion energy boundary state version {}", version);
      }
      if (version == 0) {
        initialiseBoundaryState();
      }
    }
  }
  if (radial_model.walls) {
    set_with_attrs(state["ion_energy_boundary_state_version"], 1,
                   {{"description", "Restart Vort physical radial boundary cells store "
                                    "phi midpoints in normalised potential units"}});
  }
  set_with_attrs(state["phi"], phi,
                 {{"long_name", "plasma potential"}, {"source", "vorticity"}});
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
  return correctPotential(rhs, potential, *phiSolver);
}

Field3D EvolveVorticityIonEnergy::correctPotential(const Field3D& rhs,
                                                   const Field3D& potential,
                                                   Laplacian& approximate_inverse) const {
  // Homogeneous correction boundaries preserve the separately imposed phi
  // boundaries, including INVERT_SET values supplied by a restart.
  Field3D residual{zeroFrom(rhs)};
  const Field3D forward = Div_a_Grad_perp_XZ(boussinesq_coeff, potential, true);
  BOUT_FOR(i, residual.getRegion("RGN_NOBNDRY")) {
    residual[i] = (rhs[i] - forward[i]) / boussinesq_coeff[i];
  }
  const Field3D approximate = approximate_inverse.solve(residual, zeroFrom(potential));
  // Spectral solvers can leave Y guard cells uninitialised. Only copy the
  // X-Z planes they solve; communication supplies inter-process Y guards.
  Field3D correction{zeroFrom(potential)};
  BOUT_FOR(i, correction.getRegion("RGN_NOY")) { correction[i] = approximate[i]; }
  mesh->communicate(correction);
  return potential + correction;
}

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

void EvolveVorticityIonEnergy::solvePolarisationPotential(const Field3D& rhs) {
  // Preserve the existing periodic inverse/gauge; physical X uses the
  // separate mixed-boundary inverse, matching applyPolarisationBoundary.
  Laplacian& inverse = radial_model.walls ? *chiSolver : *phiSolver;
  // The DC wall value approximates d(Phi)/dt by d(phi_b)/dt. The pressure
  // potential's wall rate d(psi_bar_b)/dt is only known in finally(), so it is
  // neglected here. Corrections retain homogeneous spectral boundaries; the
  // inhomogeneous value enters through every FV residual, including the first.
  chi = bout::ion_energy::fixedPolarisationPotential(
      boussinesq_coeff, rhs, inverse, polarisation_corrections, radial_model.chi_ac,
      radial_model.walls ? &wall_rate_inner : nullptr,
      radial_model.walls ? &wall_rate_outer : nullptr);
  chi_iterations = 1 + polarisation_corrections;

  // Diagnostic only: residual against the finite-volume operator
  chi_residual = max(abs(rhs - Div_a_Grad_perp_XZ(boussinesq_coeff, chi, true)), true);
}
