
#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/constants.hxx>
#include <bout/derivs.hxx>
#include <bout/difops.hxx>
#include <bout/field.hxx>
#include <bout/field3d.hxx>
#include <bout/fv_ops.hxx>
#include <bout/globals.hxx>
#include <bout/options.hxx>
#include <bout/output.hxx>
#include <bout/output_bout_types.hxx>
#include <bout/solver.hxx>

#include "../include/component.hxx"
#include "../include/div_ops.hxx"
#include "../include/evolve_momentum.hxx"
#include "../include/guarded_options.hxx"
#include "../include/hermes_build_config.hxx"
#include "../include/hermes_utils.hxx"
#include "../include/permissions.hxx"
#include "../include/viscosity_precon.hxx"

#include <algorithm>
#include <cmath>
#include <string>

using bout::globals::mesh;

EvolveMomentum::EvolveMomentum(std::string name, Options& alloptions, Solver* solver)
    : NamedComponent(name, {readOnly("species:{name}:AA"),
                            readOnly("species:{name}:density", Regions::Interior),
                            readWrite("species:{name}:{outputs}")}),
      name(name) {

  // Evolve the momentum in time
  solver->add(NV, std::string("NV") + name);

  auto& options = alloptions[name];

  density_floor = options["density_floor"].doc("Minimum density floor").withDefault(1e-7);

  const BoutReal temperature_floor =
      options["temperature_floor"]
          .doc("Low temperature scale for low_T_diffuse_perp")
          .withDefault<BoutReal>(0.1)
      / get<BoutReal>(alloptions["units"]["eV"]);

  low_n_diffuse_perp = options["low_n_diffuse_perp"]
                           .doc("Perpendicular diffusion at low density")
                           .withDefault<bool>(false);

  pressure_floor = density_floor * temperature_floor;

  low_p_diffuse_perp = options["low_p_diffuse_perp"]
                           .doc("Perpendicular diffusion at low pressure")
                           .withDefault<bool>(false);

  bndry_flux = options["bndry_flux"]
                   .doc("Allow flows through radial boundaries")
                   .withDefault<bool>(true);

  poloidal_flows =
      options["poloidal_flows"].doc("Include poloidal ExB flow").withDefault<bool>(true);

  hyper_z = options["hyper_z"].doc("Hyper-diffusion in Z").withDefault(-1.0);

  V.setBoundary(std::string("V") + name);

  diagnose =
      options["diagnose"].doc("Output additional diagnostics?").withDefault<bool>(false);

  use_div_par_fvv = options["use_div_par_fvv"]
                        .doc("Use Div_par_fvv instead of Div_par\nOnly for MMS tests")
                        .withDefault<bool>(use_div_par_fvv);
  if ((not use_div_par_fvv) and (not alloptions["solver"]["mms"].as<bool>())) {
    throw BoutException("use_div_par_fvv is only for MMS tests");
  }

  fix_momentum_boundary_flux =
      options["fix_momentum_boundary_flux"]
          .doc("Fix Y boundary momentum flux to boundary midpoint value?")
          .withDefault<bool>(false);

  // Set to zero so set for output
  momentum_source = 0.0;
  NV_err = 0.0;

  substitutePermissions("name", {name});
  substitutePermissions("outputs", {"velocity", "momentum"});
}

void EvolveMomentum::transform_impl(GuardedOptions& state) {
  mesh->communicate(NV);

  if (NV.isFci()) {
    NV.applyParallelBoundary();
  }

  auto species = state["species"][name];

  // Not using density boundary condition
  Field3DParallel N = getNoBoundary<Field3D>(species["density"]);

  const Field3DParallel Nlim = softFloor(N, density_floor);

  const BoutReal AA = get<BoutReal>(species["AA"]); // Atomic mass

  V = Field3DParallel(NV / (AA * Nlim));
  V.applyBoundary();

  if (NV.isFci()) {
    ASSERT2(V.hasParallelSlices());
  }
  set(species["velocity"], V);

  NV_solver = NV;                   // Save the momentum as calculated by the solver
  NV = Field3DParallel(AA * N * V); // Re-calculate consistent with V and N
  // Note: Now NV and NV_solver will differ when N < density_floor
  NV_err = NV - NV_solver; // This is used in the finally() function

  if (NV.isFci()) {
    ASSERT2(NV.hasParallelSlices());
  }

  set(species["momentum"], NV);
}

void EvolveMomentum::finally(const Options& state) {

  const auto& species = state["species"][name];
  const BoutReal AA = get<BoutReal>(species["AA"]);

  // Get updated momentum with boundary conditions
  // Note: This momentum may be modified by electromagnetic terms
  NV = get<Field3D>(species["momentum"]);

  // Get the species density
  const Field3D N = get<Field3D>(species["density"]);
  // Apply a floor to the density
  const Field3D Nlim = softFloor(N.asField3DParallel(), density_floor);

  // Typical wave speed used for numerical diffusion
  Field3DParallel fastest_wave;
  if (state.isSet("fastest_wave")) {
    fastest_wave = get<Field3D>(state["fastest_wave"]);
  } else {
    const Field3D T = get<Field3D>(species["temperature"]);
    fastest_wave = sqrt(T / AA);
  }

  // Parallel flow
  V = get<Field3D>(species["velocity"]);

  if (state.isSection("fields") and state["fields"].isSet("phi")
      and species.isSet("charge")) {

    const BoutReal Z = get<BoutReal>(species["charge"]);
    if (Z != 0.0) {
      // Electrostatic potential set and species has charge
      // -> include ExB flow and parallel force

      const Field3D phi = get<Field3D>(state["fields"]["phi"]);

      ddt(NV) = -Div_n_bxGrad_f_B_XPPM(NV, phi, bndry_flux, poloidal_flows,
                                       true); // ExB drift

      // Parallel electric field
      // Force density = - Z N ∇ϕ
      ddt(NV) -= Z * N * Grad_par(phi);

      if (state["fields"].isSet("Apar")) {
        // Include a correction term for electromagnetic simulations
        const Field3D Apar = get<Field3D>(state["fields"]["Apar"]);

        const Field3D density_source = species.isSet("density_source")
                                           ? get<Field3D>(species["density_source"])
                                           : zeroFrom(Apar);

        Field3D dummy;

        // This is Z * Apar * dn/dt, keeping just leading order terms
        Field3D dndt = density_source
                       - FV::Div_par_mod<hermes::Limiter>(N, V, fastest_wave, dummy)
                       - Div_n_bxGrad_f_B_XPPM(N, phi, bndry_flux, poloidal_flows, true);
        if (low_n_diffuse_perp) {
          dndt += Div_Perp_Lap_FV_Index(
              density_floor / softFloor(N, 1e-3 * density_floor), N);
        }
        ddt(NV) += Z * Apar * dndt;
      }
      if (state["fields"].isSet("Apar_flutter")) {
        // Magnetic flutter term
        const Field3D Apar_flutter = get<Field3D>(state["fields"]["Apar_flutter"]);

        // Using the approximation for small delta-B/B
        // b dot Grad(phi) = Grad_par(phi) + [phi, Apar]
        ddt(NV) -= Z * N * bracket(phi, Apar_flutter, BRACKET_ARAKAWA);
      }
    } else {
      ddt(NV) = 0.0;
    }
  } else {
    ddt(NV) = 0.0;
  }

  // Note:
  //  - Density floor should be consistent with calculation of V
  //    otherwise energy conservation is affected
  //  - using the same operator as in density and pressure equations doesn't work
  if (use_div_par_fvv) {
    ddt(NV) -= AA
               * FV::Div_par_fvv<hermes::Limiter>(Nlim, V, fastest_wave,
                                                  fix_momentum_boundary_flux);
  } else {
    ddt(NV) -= Div_par(NV * V);
  }

  // Parallel pressure gradient
  if (species.isSet("pressure")) {
    const Field3D P = get<Field3D>(species["pressure"]);
    if (NV.isFci()) {
      ASSERT2(P.hasParallelSlices());
    }
    ddt(NV) -= Grad_par(P.asField3DParallel());
  }

  if (state.isSection("fields") and state["fields"].isSet("Apar_flutter")) {
    // Magnetic flutter term
    const Field3D Apar_flutter = get<Field3D>(state["fields"]["Apar_flutter"]);
    ddt(NV) -= Div_n_g_bxGrad_f_B_XZ(NV, V, -Apar_flutter);

    if (species.isSet("pressure")) {
      const Field3D P = get<Field3D>(species["pressure"]);
      ddt(NV) -= bracket(P, Apar_flutter, BRACKET_ARAKAWA);
    }
  }

  if (species.isSet("low_n_coeff")) {
    // Low density parallel diffusion
    const Field3D low_n_coeff = get<Field3D>(species["low_n_coeff"]);
    ddt(NV) += FV::Div_par_K_Grad_par(low_n_coeff * V, N)
               + FV::Div_par_K_Grad_par(low_n_coeff * Nlim, V);
  }

  if (low_n_diffuse_perp) {
    ddt(NV) +=
        Div_Perp_Lap_FV_Index(density_floor / softFloor(N, 1e-3 * density_floor), NV);
  }

  if (low_p_diffuse_perp) {
    const Field3D Plim =
        softFloor(get<Field3D>(species["pressure"]), 1e-3 * pressure_floor);
    ddt(NV) += Div_Perp_Lap_FV_Index(pressure_floor / Plim, NV);
  }

  if (hyper_z > 0.) {
    auto* coord = N.getCoordinates();
    ddt(NV) -= hyper_z * SQ(SQ(coord->dz())) * D4DZ4(NV);
  }

  // Other sources/sinks
  if (species.isSet("momentum_source")) {
    momentum_source = get<Field3D>(species["momentum_source"]);
    ddt(NV) += momentum_source;
  }

  // If N < density_floor then NV and NV_solver may differ
  // -> Add term to force NV_solver towards NV
  // Note: This correction is calculated in transform()
  //       because NV may be modified by electromagnetic terms
  ddt(NV) += NV_err;

  // Scale time derivatives
  if (state.isSet("scale_timederivs")) {
    ddt(NV) *= get<Field3D>(state["scale_timederivs"]);
  }

#if CHECKLEVEL >= 1
  for (const auto& i : NV.getRegion("RGN_NOBNDRY")) {
    if (!std::isfinite(ddt(NV)[i])) {
      throw BoutException("ddt(NV{}) non-finite at {}\n", name, i);
    }
  }
#endif

  if (diagnose) {
    // Save flows if they are set

    if (species.isSet("momentum_flow_xlow")) {
      flow_xlow = get<Field3D>(species["momentum_flow_xlow"]);
    }
    if (species.isSet("momentum_flux_ylow")) {
      flow_ylow = get<Field3D>(species["momentum_flow_ylow"]);
    }
  }
  // Restore NV to the value returned by the solver
  // so that restart files contain the correct values
  // Note: Copy boundary condition so dump file has correct boundary.
  NV_solver.setBoundaryTo(NV);
  NV = NV_solver;
}

void EvolveMomentum::precon(const Options& state, BoutReal gamma) {
  namespace pc = hermes::viscosity_precon;
  // First, check whether the relevant viscosity component's eta has been set.
  // If it hasn't (i.e. if the input file does not include viscosity in the model,
  // or it does but precondition=false), this precon function will exit immediately.
  // This is how the preconditioner for momentum is controlled by options in the
  // braginskii_ion/electron_viscosity section in the input file.
  const auto frozen = pc::frozen.find(name);
  if (frozen == pc::frozen.end()) {
    return;
  }

  if (state.isSet("scale_timederivs")) {
    throw BoutException("The viscosity preconditioner does not support scale_timederivs");
  }
  if constexpr (bout::build::use_metric_3d) {
    // InvertParDiv has no implementation with 3D metrics (including FCI)
    return;
  }
  if (!viscosity_inverse) {
    viscosity_inverse = InvertParDiv::create();
  }
  const auto& species = state["species"][name];
  const auto* coord = NV.getCoordinates();

  const BoutReal AA = get<BoutReal>(species["AA"]);
  const Field3D N = softFloor(get<Field3D>(species["density"]), density_floor);
  const Field2D B = DC(coord->Bxy());
  const Field2D sqrtB = sqrt(B);

  // With g = sqrtB dV and dNV = AA N dV, the linearised viscous force
  // sqrtB Div_par((eta/B) Grad_par(sqrtB dV)) gives
  //   (AA N / B) g - gamma Div_par((eta / B) Grad_par g) = r / sqrtB.
  // Coefficients are toroidal averages, as in the conduction preconditioner.
  const Field2D A = AA * DC(N) / B;
  const Field2D Bcoef = -gamma * DC(frozen->second.eta) / B;

  // Sheath-face closure. The sheath boundary sets V in the target guard cells, so
  // the face flux b (g_guard - g), with g = sqrtB dV and g_guard = closure rho g
  // (rho = sqrtB_guard / sqrtB), becomes the extra end-cell diagonal
  // -b (1 - closure rho), averaged over z; InvertParDiv's own boundary rows keep zero
  // flux. The recovered dNV still uses the mass term A.
  const Field3D V = get<Field3D>(species["velocity"]);
  const Field3D closure = sheathClosure(V);
  Field2D A_sh{0.0, mesh};
  A_sh.allocate();
  const Field2D J = DC(coord->J());
  const Field2D dy = DC(coord->dy());
  const Field2D g_22 = DC(coord->g_22());
  auto face = [&](int x, int y, int yg) {
    if (x < mesh->xstart or x > mesh->xend) {
      return;
    }
    const BoutReal b = 0.5 * (Bcoef(x, y) + Bcoef(x, yg)) * 0.5 * (J(x, y) + J(x, yg))
                       / (0.5 * (dy(x, y) + dy(x, yg)) * 0.5 * (g_22(x, y) + g_22(x, yg))
                          * J(x, y) * dy(x, y));
    const BoutReal rho = sqrtB(x, yg) / sqrtB(x, y);
    BoutReal factor = 0.0;
    for (int z = 0; z < mesh->LocalNz; z++) {
      factor += 1.0 - closure(x, yg, z) * rho;
    }
    A_sh(x, y) -= b * factor / mesh->LocalNz;
  };
  for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
    face(r.ind, mesh->ystart, mesh->ystart - 1);
  }
  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
    face(r.ind, mesh->yend, mesh->yend + 1);
  }
  viscosity_inverse->setCoefA(A + A_sh);
  viscosity_inverse->setCoefB(Bcoef);

  Field3D rhs = ddt(NV) / sqrtB;
  rhs.applyBoundary("neumann");
  // Use the local density (i.e. non-toroidally averaged) to reconstruct
  // the approximate residual
  ddt(NV) = (AA * N / sqrtB) * viscosity_inverse->solve(rhs);
}

Field3D EvolveMomentum::sheathClosure(const Field3D& V) const {
  // Account for the smooth Bohm boundary condition in the viscsotiy preconditioner
  const auto smooth = hermes::viscosity_precon::sheath_closure.find(name);
  if (smooth != hermes::viscosity_precon::sheath_closure.end()) {
    return smooth->second;
  }
  // The sheath boundary sets V_guard = 2 v_sheath - V in field-aligned coordinates,
  // with v_sheath independent of V, so dV_guard = -dV; ions leaving supersonically
  // get V_guard = V, so dV_guard = dV.
  const Field3D V_fa = toFieldAligned(V);
  Field3D closure{zeroFrom(V_fa)};
  auto classify = [&](int x, int yi, int yg) {
    BoutReal scale = 0.0;
    for (int z = 0; z < mesh->LocalNz; z++) {
      scale = std::max({scale, std::abs(V_fa(x, yi, z)), std::abs(V_fa(x, yg, z))});
    }
    for (int z = 0; z < mesh->LocalNz; z++) {
      closure(x, yg, z) =
          std::abs(V_fa(x, yg, z) - V_fa(x, yi, z)) <= 1e-8 * scale ? 1.0 : -1.0;
    }
  };
  for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
    classify(r.ind, mesh->ystart, mesh->ystart - 1);
  }
  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
    classify(r.ind, mesh->yend, mesh->yend + 1);
  }
  return closure;
}

void EvolveMomentum::outputVars(Options& state) {
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto Cs0 = get<BoutReal>(state["Cs0"]);

  state[std::string("NV") + name].setAttributes(
      {{"time_dimension", "t"},
       {"units", "kg / m^2 / s"},
       {"conversion", SI::Mp * Nnorm * Cs0},
       {"standard_name", "momentum"},
       {"long_name", name + " parallel momentum"},
       {"species", name},
       {"source", "evolve_momentum"}});

  if (diagnose) {
    set_with_attrs(state[std::string("V") + name], V,
                   {{"time_dimension", "t"},
                    {"units", "m / s"},
                    {"conversion", Cs0},
                    {"long_name", name + " parallel velocity"},
                    {"standard_name", "velocity"},
                    {"species", name},
                    {"source", "evolve_momentum"}});

    set_with_attrs(state[std::string("ddt(NV") + name + std::string(")")], ddt(NV),
                   {{"time_dimension", "t"},
                    {"units", "kg m^-2 s^-2"},
                    {"conversion", SI::Mp * Nnorm * Cs0 * Omega_ci},
                    {"long_name", std::string("Rate of change of ") + name + " momentum"},
                    {"species", name},
                    {"source", "evolve_momentum"}});

    set_with_attrs(state[std::string("SNV") + name], momentum_source,
                   {{"time_dimension", "t"},
                    {"units", "kg m^-2 s^-2"},
                    {"conversion", SI::Mp * Nnorm * Cs0 * Omega_ci},
                    {"standard_name", "momentum source"},
                    {"long_name", name + " momentum source"},
                    {"species", name},
                    {"source", "evolve_momentum"}});

    // If fluxes have been set then add them to the output
    auto rho_s0 = get<BoutReal>(state["rho_s0"]);

    if (flow_xlow.isAllocated()) {
      set_with_attrs(
          state[fmt::format("mf{}_tot_xlow", name)], flow_xlow,
          {{"time_dimension", "t"},
           {"units", "N"},
           {"conversion", rho_s0 * SQ(rho_s0) * SI::Mp * Nnorm * Cs0 * Omega_ci},
           {"standard_name", "momentum flow"},
           {"long_name", name + " momentum flow in X. Note: May be incomplete."},
           {"species", name},
           {"source", "evolve_momentum"}});
    }
    if (flow_ylow.isAllocated()) {
      set_with_attrs(
          state[fmt::format("mf{}_tot_ylow", name)], flow_ylow,
          {{"time_dimension", "t"},
           {"units", "N"},
           {"conversion", rho_s0 * SQ(rho_s0) * SI::Mp * Nnorm * Cs0 * Omega_ci},
           {"standard_name", "momentum flow"},
           {"long_name", name + " momentum flow in Y. Note: May be incomplete."},
           {"species", name},
           {"source", "evolve_momentum"}});
    }
  }
}
