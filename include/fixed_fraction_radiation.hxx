#pragma once
#ifndef FIXED_FRACTION_RADIATION_H
#define FIXED_FRACTION_RADIATION_H

#include <bout/constants.hxx>
#include "component.hxx"

namespace {
  /// Carbon in coronal equilibrium
  /// From I.H.Hutchinson Nucl. Fusion 34 (10) 1337 - 1348 (1994)
  ///
  /// Input: Electron temperature [eV]
  /// Output: Radiation cooling curve [Wm^3]
  ///         (multiply by Ne * Ni to get W/m^3)
  struct HutchinsonCarbon {
    BoutReal curve(BoutReal Te) {
      if (Te < 0.0) {
        return 0.0;
      }
      return 2e-31*pow(Te/10., 3) / (1. + pow(Te/10., 4.5));
    }
  };

  /// Below fits are from ADAS data by Mike Kryjak 03/06/2023 and supersede above deprecated fits
  /// Data generated using radas https://github.com/cfs-energy/radas
  /// Chose N = 1E20m-3 and tau = 0.5ms based on Moulton, 2021 (DOI: 10.1088/1741-4326/abe4b2)
  /// Those values are applicable in the ITER scrape-off layer but may not be valid in other conditions.
  /// The fits are 10 coefficient polynomials fitted in log-log space like the AMJUEL database in EIRENE

  /// Argon
  struct Argon_adas{
    BoutReal curve(BoutReal Te) {
      BoutReal logT = log(Te);
      BoutReal log_out = 0;

      if (Te >= 1.5 and Te <= 1500) {
        log_out = log_out
        -8.4541e+01 * pow(logT, 0)
        +1.5773e+01 * pow(logT, 1)
        -1.5426e+01 * pow(logT, 2)
        +1.4941e+01 * pow(logT, 3)
        -1.0482e+01 * pow(logT, 4)
        +5.0092e+00 * pow(logT, 5)
        -1.6003e+00 * pow(logT, 6)
        +3.2946e-01 * pow(logT, 7)
        -4.1404e-02 * pow(logT, 8)
        +2.8706e-03 * pow(logT, 9)
        -8.3889e-05 * pow(logT, 10);
        return exp(log_out);

      } else if (Te < 1.5) {
        return 1.9535e-35;   
      } else {
        return 1.2265e-32;
      }
    }
  };

  /// Neon
  struct Neon_adas{
    BoutReal curve(BoutReal Te) {
      BoutReal logT = log(Te);
      BoutReal log_out = 0;

      if (Te >= 2 and Te <= 1000) {
        log_out = log_out
        -8.2148e+01 * pow(logT, 0)
        +1.2893e+01 * pow(logT, 1)
        -4.7427e+01 * pow(logT, 2)
        +7.4522e+01 * pow(logT, 3)
        -5.7571e+01 * pow(logT, 4)
        +2.5738e+01 * pow(logT, 5)
        -7.1276e+00 * pow(logT, 6)
        +1.2429e+00 * pow(logT, 7)
        -1.3294e-01 * pow(logT, 8)
        +7.9737e-03 * pow(logT, 9)
        -2.0549e-04 * pow(logT, 10);
        return exp(log_out);

      } else if (Te < 2) {
        return 6.3530e-36;   
      } else {
        return 1.1789e-32;
      }
    }
  };


  /// Nitrogen
  struct Nitrogen_adas{
    BoutReal curve(BoutReal Te) {
      BoutReal logT = log(Te);
      BoutReal log_out = 0;

      if (Te >= 2 and Te <= 500) {
        log_out = log_out
        -5.0165e+01 * pow(logT, 0)
        -1.3575e+02 * pow(logT, 1)
        +2.7351e+02 * pow(logT, 2)
        -2.9211e+02 * pow(logT, 3)
        +1.9012e+02 * pow(logT, 4)
        -7.9516e+01 * pow(logT, 5)
        +2.1776e+01 * pow(logT, 6)
        -3.8833e+00 * pow(logT, 7)
        +4.3473e-01 * pow(logT, 8)
        -2.7768e-02 * pow(logT, 9)
        +7.7272e-04 * pow(logT, 10);
        return exp(log_out);

      } else if (Te < 2) {
        return 4.3484e-34;   
      } else {
        return 8.1110e-33;
      }
    }
  };


  /// Carbon
  struct Carbon_adas{
    BoutReal curve(BoutReal Te) {
      BoutReal logT = log(Te);
      BoutReal log_out = 0;

      if (Te >= 1 and Te <= 500) {
        log_out = log_out
        -7.8784e+01 * pow(logT, 0)
        +1.5533e+00 * pow(logT, 1)
        +1.6590e+01 * pow(logT, 2)
        -3.2380e+01 * pow(logT, 3)
        +3.1278e+01 * pow(logT, 4)
        -1.7483e+01 * pow(logT, 5)
        +5.9139e+00 * pow(logT, 6)
        -1.2297e+00 * pow(logT, 7)
        +1.5400e-01 * pow(logT, 8)
        -1.0680e-02 * pow(logT, 9)
        +3.1566e-04 * pow(logT, 10);
        return exp(log_out);

      } else if (Te < 1) {
        return 6.0062e-35;   
      } else {
        return 4.5306e-33;
      }
    }
  };
}

/// Set ion densities from electron densities
///
template <typename CoolingCurve>
struct FixedFractionRadiation : public Component {
  /// Inputs
  /// - <name>
  ///   - fraction
  FixedFractionRadiation(std::string name, Options &alloptions, Solver *UNUSED(solver)) : name(name) {
    auto& options = alloptions[name];

    fraction = options["fraction"]
      .doc("Impurity ion density as fraction of electron density")
      .withDefault(0.0);

    diagnose = options["diagnose"]
      .doc("Output radiation diagnostic?")
      .withDefault<bool>(false);

    // Get the units
    auto& units = alloptions["units"];
    Tnorm = get<BoutReal>(units["eV"]);
    Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
    FreqNorm = 1. / get<BoutReal>(units["seconds"]);
  }

  /// Required inputs
  ///
  /// - species
  ///   - e
  ///     - density
  ///     - temperature
  ///
  /// Sets the electron energy loss
  ///
  /// - species
  ///   - e
  ///     - energy_source
  ///
  void transform(Options &state) override {
    auto& electrons = state["species"]["e"];
    // Don't need boundary cells
    const Field3D Ne = GET_NOBOUNDARY(Field3D, electrons["density"]);
    const Field3D Te = GET_NOBOUNDARY(Field3D, electrons["temperature"]);

    radiation = cellAverage(
                            [&](BoutReal ne, BoutReal te) {
                              if (ne < 0.0 or te < 0.0) {
                                return 0.0;
                              }
                              // Fixed fraction ions
                              const BoutReal ni = fraction * ne;
                              // cooling in Wm^3 so normalise.
                              // Note factor of qe due to Watts rather than eV
                              return ne * ni * cooling.curve(te * Tnorm) * Nnorm /
                                (SI::qe * Tnorm * FreqNorm);
                            },
                            Ne.getRegion("RGN_NOBNDRY"))(Ne, Te);

    // Remove radiation from the electron energy source
    subtract(electrons["energy_source"], radiation);
  }

  void outputVars(Options& state) override {
    AUTO_TRACE();

    if (diagnose) {
      set_with_attrs(state[std::string("R") + name], radiation,
                     {{"time_dimension", "t"},
                      {"units", "W / m^3"},
                      {"conversion", SI::qe * Tnorm * Nnorm * FreqNorm},
                      {"long_name", std::string("Radiation cooling ") + name},
                      {"source", "fixed_fraction_radiation"}});
    }
  }
 private:
  std::string name;

  CoolingCurve cooling; ///< The cooling curve L(T) -> Wm^3
  BoutReal fraction; ///< Fixed fraction

  bool diagnose; ///< Output radiation diagnostic?
  Field3D radiation; ///< For output diagnostic

  // Normalisations
  BoutReal Tnorm, Nnorm, FreqNorm;
};

namespace {
  RegisterComponent<FixedFractionRadiation<HutchinsonCarbon>>
    registercomponentfixedfractionhutchinsonbcarbon("fixed_fraction_hutchinson_carbon");
    
  RegisterComponent<FixedFractionRadiation<Carbon_adas>>
    registercomponentfixedfractioncarbon("fixed_fraction_carbon");  

  RegisterComponent<FixedFractionRadiation<Nitrogen_adas>>
    registercomponentfixedfractionnitrogen("fixed_fraction_nitrogen");

  RegisterComponent<FixedFractionRadiation<Neon_adas>>
    registercomponentfixedfractionneon("fixed_fraction_neon");

  RegisterComponent<FixedFractionRadiation<Argon_adas>>
    registercomponentfixedfractionargon("fixed_fraction_argon");

}

#endif // FIXED_FRACTION_IONS_H
