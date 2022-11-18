
#include <derivs.hxx>
#include <difops.hxx>
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

  if (options["diagnose"]
          .doc("Output additional diagnostics?")
          .withDefault<bool>(false)) {
    bout::globals::dump.addRepeat(V, std::string("V") + name);

    bout::globals::dump.addRepeat(ddt(NV), std::string("ddt(NV") + name + std::string(")"));
    bout::globals::dump.addRepeat(momentum_source, std::string("SNV") + name);
    momentum_source = 0.0;
  }
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

