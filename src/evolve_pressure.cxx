
#include <bout/constants.hxx>
#include <bout/fv_ops.hxx>
#include <difops.hxx>

#include "../include/div_ops.hxx"
#include "../include/evolve_pressure.hxx"

using bout::globals::mesh;

namespace {
Ind3D indexAt(const Field3D& f, int x, int y, int z) {
  int ny = f.getNy();
  int nz = f.getNz();
  return Ind3D{(x * ny + y) * nz + z, ny, nz};
}
} // namespace

EvolvePressure::EvolvePressure(std::string name, Options& alloptions, Solver* solver)
    : name(name) {
  AUTO_TRACE();

  solver->add(P, std::string("P") + name);

  auto& options = alloptions[name];

  bndry_flux = options["bndry_flux"]
                   .doc("Allow flows through radial boundaries")
                   .withDefault<bool>(true);

  poloidal_flows =
      options["poloidal_flows"].doc("Include poloidal ExB flow").withDefault<bool>(true);

  thermal_conduction = options["thermal_conduction"]
                           .doc("Include parallel heat conduction?")
                           .withDefault<bool>(true);

  p_div_v = options["p_div_v"]
                .doc("Use p*Div(v) form? Default, false => v * Grad(p) form")
                .withDefault<bool>(false);

  if (options["diagnose"]
          .doc("Save additional output diagnostics")
          .withDefault<bool>(false)) {
    if (thermal_conduction) {
      bout::globals::dump.addRepeat(kappa_par, std::string("kappa_par_") + name);
    }
    bout::globals::dump.addRepeat(T, std::string("T") + name);

    bout::globals::dump.addRepeat(ddt(P), std::string("ddt(P") + name + std::string(")"));
    bout::globals::dump.addRepeat(Sp, std::string("SP") + name);
    Sp = 0.0;
  }

  const Options& units = alloptions["units"];
  const BoutReal Nnorm = units["inv_meters_cubed"];
  const BoutReal Tnorm = units["eV"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();

  source = alloptions[std::string("P") + name]["source"]
               .doc(std::string("Source term in ddt(P") + name
                    + std::string("). Units [N/m^2/s]"))
               .withDefault(Field3D(0.0))
           / (SI::qe * Nnorm * Tnorm * Omega_ci);
}

void EvolvePressure::transform(Options& state) {
  AUTO_TRACE();

  mesh->communicate(P);

  auto& species = state["species"][name];

  // Check for low pressures, ensure that Pi >= 0
  P = floor(P, 0.0);

  set(species["pressure"], P);

  // Calculate temperature
  // Not using density boundary condition
  N = getNoBoundary<Field3D>(species["density"]);
  T = P / floor(N, 1e-5);
  T.applyBoundary("neumann");

  set(species["temperature"], T);
}

void EvolvePressure::finally(const Options& state) {
  AUTO_TRACE();

  /// Get the section containing this species
  const auto& species = state["species"][name];

  if (state.isSection("fields") and state["fields"].isSet("phi")) {
    // Electrostatic potential set -> include ExB flow

    Field3D phi = get<Field3D>(state["fields"]["phi"]);

    ddt(P) = -Div_n_bxGrad_f_B_XPPM(P, phi, bndry_flux, poloidal_flows, true);
  } else {
    ddt(P) = 0.0;
  }

  if (species.isSet("velocity")) {
    Field3D V = get<Field3D>(species["velocity"]);

    // Typical wave speed used for numerical diffusion
    Field3D sound_speed;
    if (state.isSet("sound_speed")) {
      Field3D sound_speed = get<Field3D>(state["sound_speed"]);
    } else {
      Field3D T = get<Field3D>(species["temperature"]);
      sound_speed = sqrt(T);
    }

    if (p_div_v) {
      // Use the P * Div(V) form
      ddt(P) -= FV::Div_par(P, V, sound_speed);

      // Work done. This balances energetically a term in the momentum equation
      ddt(P) -= (2. / 3) * P * Div_par(V);

    } else {
      // Use V * Grad(P) form
      // Note: A mixed form has been tried (on 1D neon example)
      //       -(4/3)*FV::Div_par(P,V) + (1/3)*(V * Grad_par(P) - P * Div_par(V))
      //       Caused heating of charged species near sheath like p_div_v
      ddt(P) -= (5. / 3) * FV::Div_par(P, V, sound_speed);

      ddt(P) += (2. / 3) * V * Grad_par(P);
    }
  }

  // Parallel heat conduction
  if (thermal_conduction) {

    // Calculate ion collision times
    const Field3D tau = 1. / get<Field3D>(species["collision_frequency"]);
    const BoutReal AA = get<BoutReal>(species["AA"]); // Atomic mass

    // Parallel heat conduction
    // Braginskii expression for parallel conduction
    // kappa ~ n * v_th^2 * tau
    //
    // Note: Coefficient is slightly different for electrons (3.16) and ions (3.9)
    kappa_par = 3.9 * P * tau / AA;

    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        auto i = indexAt(kappa_par, r.ind, mesh->ystart, jz);
        auto im = i.ym();
        kappa_par[im] = kappa_par[i];
      }
    }
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        auto i = indexAt(kappa_par, r.ind, mesh->yend, jz);
        auto ip = i.yp();
        kappa_par[ip] = kappa_par[i];
      }
    }

    // Note: Flux through boundary turned off, because sheath heat flux
    // is calculated and removed separately
    ddt(P) += (2. / 3) * FV::Div_par_K_Grad_par(kappa_par, T, false);
  }

  //////////////////////
  // Other sources

  Sp = source;
  if (species.isSet("energy_source")) {
    Sp += (2. / 3) * get<Field3D>(species["energy_source"]); // For diagnostic output
  }
  ddt(P) += Sp;

#if CHECK > 1
  bout::checkFinite(ddt(P), std::string("ddt P") + name, "RGN_NOBNDRY");
#endif
}
