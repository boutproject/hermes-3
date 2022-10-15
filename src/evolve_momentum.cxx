
#include <derivs.hxx>
#include <difops.hxx>
#include <bout/constants.hxx>
#include <bout/fv_ops.hxx>
#include <bout/output_bout_types.hxx>

#include "../include/evolve_momentum.hxx"
#include "../include/div_ops.hxx"

namespace {
BoutReal floor(BoutReal value, BoutReal min) {
  if (value < min)
    return min;
  return value;
}
}

using bout::globals::mesh;

namespace FV {
  template<typename CellEdges = MC>
  const Field3D Div_par_fvv(const Field3D &f_in, const Field3D &v_in,
                            const Field3D &wave_speed_in, bool fixflux=true) {

    ASSERT1(areFieldsCompatible(f_in, v_in));
    ASSERT1(areFieldsCompatible(f_in, wave_speed_in));

    Mesh* mesh = f_in.getMesh();

    CellEdges cellboundary;

    /// Ensure that f, v and wave_speed are field aligned
    Field3D f = toFieldAligned(f_in, "RGN_NOX");
    Field3D v = toFieldAligned(v_in, "RGN_NOX");
    Field3D wave_speed = toFieldAligned(wave_speed_in, "RGN_NOX");

    Coordinates *coord = f_in.getCoordinates();

    Field3D result{zeroFrom(f)};
    
    // Only need one guard cell, so no need to communicate fluxes
    // Instead calculate in guard cells to preserve fluxes
    int ys = mesh->ystart-1;
    int ye = mesh->yend+1;

    for (int i = mesh->xstart; i <= mesh->xend; i++) {

      if (!mesh->firstY(i) || mesh->periodicY(i)) {
        // Calculate in guard cell to get fluxes consistent between processors
        ys = mesh->ystart - 1;
      } else {
        // Don't include the boundary cell. Note that this implies special
        // handling of boundaries later
        ys = mesh->ystart;
      }

      if (!mesh->lastY(i) || mesh->periodicY(i)) {
        // Calculate in guard cells
        ye = mesh->yend + 1;
      } else {
        // Not in boundary cells
        ye = mesh->yend;
      }

      for (int j = ys; j <= ye; j++) {
        // Pre-calculate factors which multiply fluxes

        // For right cell boundaries
        BoutReal common_factor = (coord->J(i, j) + coord->J(i, j + 1)) /
          (sqrt(coord->g_22(i, j)) + sqrt(coord->g_22(i, j + 1)));
        
        BoutReal flux_factor_rc = common_factor / (coord->dy(i, j) * coord->J(i, j));
        BoutReal flux_factor_rp = common_factor / (coord->dy(i, j + 1) * coord->J(i, j + 1));

        // For left cell boundaries
        common_factor = (coord->J(i, j) + coord->J(i, j - 1)) /
          (sqrt(coord->g_22(i, j)) + sqrt(coord->g_22(i, j - 1)));

        BoutReal flux_factor_lc = common_factor / (coord->dy(i, j) * coord->J(i, j));
        BoutReal flux_factor_lm = common_factor / (coord->dy(i, j - 1) * coord->J(i, j - 1));
        
        for (int k = 0; k < mesh->LocalNz; k++) {

          ////////////////////////////////////////////
          // Reconstruct f at the cell faces
          // This calculates s.R and s.L for the Right and Left
          // face values on this cell
          
          // Reconstruct f at the cell faces
          Stencil1D s;
          s.c = f(i, j, k);
          s.m = f(i, j - 1, k);
          s.p = f(i, j + 1, k);
          
          cellboundary(s); // Calculate s.R and s.L

          // Reconstruct v at the cell faces
          Stencil1D sv;
          sv.c = v(i, j, k);
          sv.m = v(i, j - 1, k);
          sv.p = v(i, j + 1, k);
          
          cellboundary(sv);
          
          ////////////////////////////////////////////
          // Right boundary

          // Calculate velocity at right boundary (y+1/2)
          BoutReal vpar = 0.5 * (v(i, j, k) + v(i, j + 1, k));
          BoutReal flux;

          if (mesh->lastY(i) && (j == mesh->yend) && !mesh->periodicY(i)) {
            // Last point in domain

            BoutReal bndryval = 0.5 * (s.c + s.p);
            if (fixflux) {
              // Use mid-point to be consistent with boundary conditions
              flux = bndryval * vpar * vpar;
            } else {
              // Handle pathological case of diverging flows
              flux = bndryval * vpar * BOUTMIN(vpar, BOUTMAX(v(i, j, k), 0.0));
            }
          } else {
            
            // Maximum wave speed in the two cells
            BoutReal amax = BOUTMAX(wave_speed(i, j, k), wave_speed(i, j + 1, k));

            if (vpar > amax) {
              // Supersonic flow out of this cell
              flux = s.R * vpar * sv.R;
            } else if (vpar < -amax) {
              // Supersonic flow into this cell
              flux = 0.0;
            } else {
              // Subsonic flow, so a mix of right and left fluxes
              flux = s.R * 0.5 * (vpar + amax) * sv.R;
            }
          }
          
          result(i, j, k) += flux * flux_factor_rc;
          result(i, j + 1, k) -= flux * flux_factor_rp;

          ////////////////////////////////////////////
          // Calculate at left boundary
          
          vpar = 0.5 * (v(i, j, k) + v(i, j - 1, k));

          if (mesh->firstY(i) && (j == mesh->ystart) && !mesh->periodicY(i)) {
            // First point in domain
            BoutReal bndryval = 0.5 * (s.c + s.m);
            if (fixflux) {
              // Use mid-point to be consistent with boundary conditions
              flux = bndryval * vpar * vpar;
            } else {
              // Add flux due to difference in boundary values
              flux = s.L * vpar * sv.L - wave_speed(i, j, k) * (s.L * sv.L - bndryval * vpar);
            }
          } else {
            
            // Maximum wave speed in the two cells
            BoutReal amax = BOUTMAX(wave_speed(i, j, k), wave_speed(i, j - 1, k));

            if (vpar < -amax) {
              // Supersonic out of this cell
              flux = s.L * vpar * sv.L;
            } else if (vpar > amax) {
              // Supersonic into this cell
              flux = 0.0;
            } else {
              flux = s.L * 0.5 * (vpar - amax) * sv.L;
            }
          }
          
          result(i, j, k) -= flux * flux_factor_lc;
          result(i, j - 1, k) += flux * flux_factor_lm;
          
        }
      }
    }
    return fromFieldAligned(result, "RGN_NOBNDRY");
  }
  
}

EvolveMomentum::EvolveMomentum(std::string name, Options &alloptions, Solver *solver) : name(name) {
  AUTO_TRACE();
  
  // Evolve the momentum in time
  solver->add(NV, std::string("NV") + name);

  auto& options = alloptions[name];

  density_floor = options["density_floor"].doc("Minimum density floor").withDefault(1e-5);

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
}

void EvolveMomentum::transform(Options &state) {
  AUTO_TRACE();
  mesh->communicate(NV);

  auto& species = state["species"][name];

  set(species["momentum"], NV);

  // Not using density boundary condition
  Field3D N = floor(getNoBoundary<Field3D>(species["density"]), density_floor);
  BoutReal AA = get<BoutReal>(species["AA"]); // Atomic mass

  V = NV / (AA * N);
  V.applyBoundary();
  set(species["velocity"], V);
}

void EvolveMomentum::finally(const Options &state) {
  AUTO_TRACE();

  auto& species = state["species"][name];

  // Get updated momentum with boundary conditions
  NV = get<Field3D>(species["momentum"]);

  // Get the species density
  Field3D N = get<Field3D>(species["density"]);

  if (state.isSection("fields") and state["fields"].isSet("phi")
      and species.isSet("charge")) {

    BoutReal Z = get<BoutReal>(species["charge"]);
    if (Z != 0.0) {
      // Electrostatic potential set and species has charge
      // -> include ExB flow and parallel force

      Field3D phi = get<Field3D>(state["fields"]["phi"]);

      ddt(NV) = -Div_n_bxGrad_f_B_XPPM(NV, phi, bndry_flux, poloidal_flows,
                                       true); // ExB drift

      // Parallel electric field
      // Force density = - Z N ∇ϕ
      ddt(NV) -= Z * N * Grad_par(phi);
    }
  } else {
    ddt(NV) = 0.0;
  }

  // Parallel flow
  V = get<Field3D>(species["velocity"]);

  // Typical wave speed used for numerical diffusion
  Field3D T = get<Field3D>(species["temperature"]);
  BoutReal AA = get<BoutReal>(species["AA"]);
  Field3D sound_speed = sqrt(T / AA);

  // Note: Density floor should be consistent with calculation of V
  //       otherwise energy conservation is affected
  ddt(NV) -= FV::Div_par_fvv(floor(N, density_floor), V, sound_speed);

  // Parallel pressure gradient
  if (species.isSet("pressure")) {
    Field3D P = get<Field3D>(species["pressure"]);
    ddt(NV) -= Grad_par(P);
  }

  if (species.isSet("low_n_coeff")) {
    // Low density parallel diffusion
    Field3D low_n_coeff = get<Field3D>(species["low_n_coeff"]);
    ddt(NV) += FV::Div_par_K_Grad_par(low_n_coeff * V, N) + FV::Div_par_K_Grad_par(low_n_coeff * floor(N, density_floor), V);
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

#if CHECKLEVEL >= 1
  for (auto& i : NV.getRegion("RGN_NOBNDRY")) {
    if (!std::isfinite(ddt(NV)[i])) {
      throw BoutException("ddt(NV{}) non-finite at {}\n", name, i);
    }
  }
#endif
}

void EvolveMomentum::outputVars(Options &state) {
  AUTO_TRACE();
  // Normalisations
  auto Nnorm = state["Nnorm"].as<BoutReal>();
  auto Omega_ci = state["Omega_ci"].as<BoutReal>();
  auto Cs0 = state["Cs0"].as<BoutReal>();

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
  }
}
