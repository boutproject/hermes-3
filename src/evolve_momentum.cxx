
#include <bout/derivs.hxx>
#include <bout/difops.hxx>
#include <bout/constants.hxx>
#include <bout/fv_ops.hxx>
#include <bout/output_bout_types.hxx>

#include "../include/evolve_momentum.hxx"
#include "../include/div_ops.hxx"
#include "../include/hermes_build_config.hxx"

namespace {
BoutReal floor(BoutReal value, BoutReal min) {
  if (value < min)
    return min;
  return value;
}
}

using bout::globals::mesh;

EvolveMomentum::EvolveMomentum(std::string name, Options &alloptions, Solver *solver) : name(name) {
  AUTO_TRACE();
  
  // Evolve the momentum in time
  solver->add(NV, std::string("NV") + name);

  auto& options = alloptions[name];

  density_floor = options["density_floor"].doc("Minimum density floor").withDefault(1e-7);

  low_n_diffuse_perp = options["low_n_diffuse_perp"]
                           .doc("Perpendicular diffusion at low density")
                           .withDefault<bool>(false);

  pressure_floor = density_floor * (1./get<BoutReal>(alloptions["units"]["eV"]));

  low_p_diffuse_perp = options["low_p_diffuse_perp"]
                           .doc("Perpendicular diffusion at low pressure")
                           .withDefault<bool>(false);

  bndry_flux = options["bndry_flux"]
                      .doc("Allow flows through radial boundaries")
                      .withDefault<bool>(true);

  poloidal_flows = options["poloidal_flows"]
                       .doc("Include poloidal ExB flow")
                       .withDefault<bool>(true);

  hyper_z = options["hyper_z"].doc("Hyper-diffusion in Z").withDefault(-1.0);

  V.setBoundary(std::string("V") + name);

  diagnose = options["diagnose"]
    .doc("Output additional diagnostics?")
    .withDefault<bool>(false);

  fix_momentum_boundary_flux = options["fix_momentum_boundary_flux"]
    .doc("Fix Y boundary momentum flux to boundary midpoint value?")
    .withDefault<bool>(false);

  // Set to zero so set for output
  momentum_source = 0.0;
}

void EvolveMomentum::transform(Options &state) {
  AUTO_TRACE();
  mesh->communicate(NV);

  auto& species = state["species"][name];

  // Not using density boundary condition
  auto N = getNoBoundary<Field3D>(species["density"]);
  Field3D Nlim = floor(N, density_floor);
  BoutReal AA = get<BoutReal>(species["AA"]); // Atomic mass

  V = NV / (AA * Nlim);
  V.applyBoundary();
  set(species["velocity"], V);

  NV_solver = NV; // Save the momentum as calculated by the solver
  NV = AA * N * V; // Re-calculate consistent with V and N
  // Note: Now NV and NV_solver will differ when N < density_floor
  NV_err = NV - NV_solver; // This is used in the finally() function
  set(species["momentum"], NV);
}

void EvolveMomentum::finally(const Options &state) {
  AUTO_TRACE();

  auto& species = state["species"][name];
  BoutReal AA = get<BoutReal>(species["AA"]);

  // Get updated momentum with boundary conditions
  // Note: This momentum may be modified by electromagnetic terms
  NV = get<Field3D>(species["momentum"]);

  // Get the species density
  Field3D N = get<Field3D>(species["density"]);
  // Apply a floor to the density
  Field3D Nlim = floor(N, density_floor);

  // Typical wave speed used for numerical diffusion
  Field3D fastest_wave;
  if (state.isSet("fastest_wave")) {
    fastest_wave = get<Field3D>(state["fastest_wave"]);
  } else {
    Field3D T = get<Field3D>(species["temperature"]);
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

        const Field3D density_source = species.isSet("density_source") ?
          get<Field3D>(species["density_source"])
          : zeroFrom(Apar);

        Field3D dummy;
        
        // This is Z * Apar * dn/dt, keeping just leading order terms
        Field3D dndt = density_source
          - FV::Div_par_mod<hermes::Limiter>(N, V, fastest_wave, dummy)
          - Div_n_bxGrad_f_B_XPPM(N, phi, bndry_flux, poloidal_flows, true)
          ;
        if (low_n_diffuse_perp) {
          dndt += Div_Perp_Lap_FV_Index(density_floor / floor(N, 1e-3 * density_floor), N,
                                        bndry_flux);
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
  ddt(NV) -= AA * FV::Div_par_fvv<hermes::Limiter>(Nlim, V, fastest_wave, fix_momentum_boundary_flux);
  
  // Parallel pressure gradient
  if (species.isSet("pressure")) {
    Field3D P = get<Field3D>(species["pressure"]);
    ddt(NV) -= Grad_par(P);
  }

  if (state.isSection("fields") and state["fields"].isSet("Apar_flutter")) {
    // Magnetic flutter term
    const Field3D Apar_flutter = get<Field3D>(state["fields"]["Apar_flutter"]);
    ddt(NV) -= Div_n_g_bxGrad_f_B_XZ(NV, V, -Apar_flutter);

    if (species.isSet("pressure")) {
      Field3D P = get<Field3D>(species["pressure"]);
      ddt(NV) -= bracket(P, Apar_flutter, BRACKET_ARAKAWA);
    }
  }

  if (species.isSet("low_n_coeff")) {
    // Low density parallel diffusion
    Field3D low_n_coeff = get<Field3D>(species["low_n_coeff"]);
    ddt(NV) += FV::Div_par_K_Grad_par(low_n_coeff * V, N) + FV::Div_par_K_Grad_par(low_n_coeff * Nlim, V);
  }

  if (low_n_diffuse_perp) {
    ddt(NV) += Div_Perp_Lap_FV_Index(density_floor / floor(N, 1e-3 * density_floor), NV, true);
  }

  if (low_p_diffuse_perp) {
    Field3D Plim = floor(get<Field3D>(species["pressure"]), 1e-3 * pressure_floor);
    ddt(NV) += Div_Perp_Lap_FV_Index(pressure_floor / Plim, NV, true);
  }

  if (hyper_z > 0.) {
    auto* coord = N.getCoordinates();
    ddt(NV) -= hyper_z * SQ(SQ(coord->dz)) * D4DZ4(NV);
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
  for (auto& i : NV.getRegion("RGN_NOBNDRY")) {
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

void EvolveMomentum::outputVars(Options &state) {
  AUTO_TRACE();
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
      set_with_attrs(state[fmt::format("mf{}_tot_xlow", name)], flow_xlow,
                   {{"time_dimension", "t"},
                    {"units", "N"},
                    {"conversion", rho_s0 * SQ(rho_s0) * SI::Mp * Nnorm * Cs0 * Omega_ci},
                    {"standard_name", "momentum flow"},
                    {"long_name", name + " momentum flow in X. Note: May be incomplete."},
                    {"species", name},
                    {"source", "evolve_momentum"}});
    }
    if (flow_ylow.isAllocated()) {
      set_with_attrs(state[fmt::format("mf{}_tot_ylow", name)], flow_ylow,
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
