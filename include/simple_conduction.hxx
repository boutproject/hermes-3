#pragma once
#ifndef SIMPLE_CONDUCTION_H
#define SIMPLE_CONDUCTION_H

#include <bout/constants.hxx>
#include <bout/fv_ops.hxx>

#include "component.hxx"

/// Simplified models of parallel heat conduction
///
/// Intended mainly for testing.
///
/// Expressions taken from:
/// https://farside.ph.utexas.edu/teaching/plasma/lectures1/node35.html
struct SimpleConduction : public Component {
  SimpleConduction(std::string name, Options& alloptions, Solver*) : name(name) {
    auto& units = alloptions["units"];
    Tnorm = units["eV"];
    Nnorm = units["inv_meters_cubed"];
    BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();

    auto& options = alloptions[name];

    BoutReal kappa_coefficient =
        options["kappa_coefficient"]
            .doc("Numerical coefficient in parallel heat conduction. Default is 3.16 for "
                 "electrons, 3.9 otherwise")
            .withDefault((name == "e") ? 3.16 : 3.9);

    if (name == "e") {
      // Electrons
      kappa0 = kappa_coefficient * 6 * sqrt(2 * SI::Mp) * pow(PI * SI::qe * Tnorm, 1.5)
               * SQ(SI::e0) / (SQ(SQ(SI::qe)) * Nnorm) * Omega_ci;
    } else {
      // Ions
      kappa0 = kappa_coefficient * 12 * sqrt(SI::Mp) * pow(PI * SI::qe * Tnorm, 1.5)
               * SQ(SI::e0) / (SQ(SQ(SI::qe)) * Nnorm) * Omega_ci;
    }

    // Fix the temperature in the heat conduction coefficients
    temperature = options["conduction_temperature"]
                      .doc("Fix temperature in the calculation of the Coulomb log and "
                           "conduction coefficient. < 0 for not fixed")
                      .withDefault<BoutReal>(-1.0)
                  / Tnorm;

    density = options["conduction_density"]
                  .doc("Fix density in the calculation of the Coulomb log and conduction "
                       "coefficient. < 0 for not fixed")
                  .withDefault<BoutReal>(-1.0)
              / Nnorm;

    boundary_flux = options["conduction_boundary_flux"]
      .doc("Allow heat conduction through sheath boundaries?")
      .withDefault<bool>(false);
  }

private:
  std::string name;      ///< Name of the species e.g. "e"
  BoutReal kappa0;       ///< Pre-calculated constant in heat conduction coefficient
  BoutReal Nnorm, Tnorm; ///< Normalisation coefficients

  BoutReal temperature; ///< Fix temperature if > 0
  BoutReal density;     ///< Fix density if > 0

  bool boundary_flux;   ///< Allow flux through sheath boundaries?

  void transform_impl(GuardedOptions& state) override {
    auto species = state["species"][name];

    // Species time-evolving temperature
    Field3D T = GET_NOBOUNDARY(Field3D, species["temperature"]);

    Field3D Tcoef = T; // Temperature used in coefficients.
    if (temperature > 0.0) {
      Tcoef = temperature; // Fixed
    }
    Field3D N;
    if (density > 0.0) {
      N = density;
    } else {
      N = GET_NOBOUNDARY(Field3D, species["density"]);
    }
    auto AA = get<BoutReal>(species["AA"]);

    Field3D coulomb_log = 6.6 - 0.5 * log(N * Nnorm / 1e20) + 1.5 * log(Tcoef * Tnorm);

    // Parallel heat conduction coefficient
    Field3D kappa_par = kappa0 * pow(Tcoef, 2.5) / (coulomb_log * sqrt(AA));

    // Note: Flux through boundary turned off by default, because sheath heat flux
    // is calculated and removed separately
    Field3D DivQ = FV::Div_par_K_Grad_par(kappa_par, T, boundary_flux);

    add(species["energy_source"], DivQ);
  }
};

namespace {
RegisterComponent<SimpleConduction>
    registercomponentsimpleconduction("simple_conduction");
}

#endif // SIMPLE_CONDUCTION_H
