#pragma once
#ifndef TOKAMAK_CORE_H
#define TOKAMAK_CORE_H

#include <bout/constants.hxx>
#include "component.hxx"

/// Allows to set boundary conditions on the core in a tokamak simulation
struct TokamakCore : public Component {

  /// Comments on inputs will go here
  ///  - <name> (e.g. "d+")
  ///    - input        comment

  TokamakCore(std::string name, Options& alloptions, Solver*) : name(name) {
    const auto& units = alloptions["units"];
    const BoutReal Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
    const BoutReal Tnorm = get<BoutReal>(units["eV"]);
    const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();
    const BoutReal rho_s0 = get<BoutReal>(units["meters"]);
    
    /// NOTE ON NORMALISATION AND UNITS - Mike Kryjak 20/11/2023
    // Density normalisation is technically in particles per m3 and not 1/m3, which is why rho_s0**3 != 1/Nnorm
    // Taking the example of the core particle source, we want it to be normalised to 1/(Nnorm*Omega_ci).
    // However, the volume has been normalised to rho_s0**3, so it must be taken out:
    // S_norm = S_si / (Nnorm * Omega_ci)
    // S_norm = flow_norm / volume_norm
    // -> S_si / (Nnorm * Omega_ci) = flow_norm * (rho_s0**3 / volume_si)
    // -> flow_si = flow_norm * rho_s0**3 * Nnorm * Omega_ci
    // 
    // As per the above, there must be a rho_s0**3 * Nnorm somewhere in the calculation. However, nowhere else in the code
    // applies this correction to the volume. To stay in line with this, I am applying it to the flow inputs instead.

    Options& options = alloptions[name];

    power = options["core_power"]
                               .doc("Inflow of energy in the core [W]. -1 = off")
                               .withDefault<BoutReal>(-1.0) / (SI::qe * Tnorm * Omega_ci * Nnorm * rho_s0*rho_s0*rho_s0);

    particle_flow = options["core_particle_flow"]
                               .doc("Inflow of particles in the core [s^-1]. -1 = off")
                               .withDefault<BoutReal>(-1.0) / (Omega_ci * Nnorm * rho_s0*rho_s0*rho_s0);

    diagnose = options["diagnose"]
                   .doc("Output additional diagnostics?")
                   .withDefault<bool>(false);
  }

  /// Inputs
  ///  - <name>
  ///    - density
  ///
  /// Outputs
  ///
  ///  - <name>
  ///    - density_source
  ///
  void transform(Options& state) override;

  void outputVars(Options& state) override {
    AUTO_TRACE();

    auto Nnorm = get<BoutReal>(state["Nnorm"]);
    auto Tnorm = get<BoutReal>(state["Tnorm"]);
    auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
    auto rho_s0 = get<BoutReal>(state["rho_s0"]);
    const BoutReal Pnorm = SI::qe * Tnorm * Nnorm; 

    if (diagnose) {

      if (power > 0) {
        set_with_attrs(state[std::string("E") + name + std::string("_core_src")], energy_source,
                    {{"time_dimension", "t"},
                      {"units", "W m^-3"},
                      {"conversion", Nnorm * Tnorm * SI::qe * Omega_ci},
                      {"standard_name", "energy source"},
                      {"long_name", name + " core energy source"},
                      {"species", name},
                      {"source", "tokamak_core"}});
      }

      if (particle_flow > 0) {
        set_with_attrs(state[std::string("S") + name + std::string("_core_src")], density_source,
                    {{"time_dimension", "t"},
                      {"units", "m^-3 s^-1"},
                      {"conversion", Nnorm * Omega_ci},
                      {"standard_name", "density source"},
                      {"long_name", name + " core number density source"},
                      {"species", name},
                      {"source", "tokamak_core"}});
      }

    }
  }


private:
  std::string name; ///< The species name

  BoutReal power;            ///< Core energy inflow in [W]
  Field3D energy_source;     ///< Core energy source in [W/m^3]
  BoutReal particle_flow;    ///< Core particle inflow in [s^-1]
  Field3D density_source;   ///< Core particle source in [m^-3/s]
  BoutReal core_volume_local, core_volume;      ///< Volume of first core ring in [m^3] for one processor, and for whole domain

  Field2D loopmarker; /// Debug only

  bool diagnose; ///< Output diagnostic information?
};

namespace {
RegisterComponent<TokamakCore> register_tokamak_core("tokamak_core");
}

#endif // TOKAMAK_CORE
