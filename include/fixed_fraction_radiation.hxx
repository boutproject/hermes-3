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
        -8.7276e+01 * pow(logT, 0)
        +2.6028e+01 * pow(logT, 1)
        -3.1418e+01 * pow(logT, 2)
        +2.8683e+01 * pow(logT, 3)
        -1.7736e+01 * pow(logT, 4)
        +7.4835e+00 * pow(logT, 5)
        -2.1562e+00 * pow(logT, 6)
        +4.1118e-01 * pow(logT, 7)
        -4.8969e-02 * pow(logT, 8)
        +3.2708e-03 * pow(logT, 9)
        -9.3114e-05 * pow(logT, 10);
        return exp(log_out);

      } else if (Te < 1.5) {
        return 1.2316e-35;   
      } else if (Te > 1500) {
        return 1.2572e-32;
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
        -9.5485e+01 * pow(logT, 0)
        +6.8706e+01 * pow(logT, 1)
        -1.4171e+02 * pow(logT, 2)
        +1.6144e+02 * pow(logT, 3)
        -1.0674e+02 * pow(logT, 4)
        +4.3751e+01 * pow(logT, 5)
        -1.1490e+01 * pow(logT, 6)
        +1.9365e+00 * pow(logT, 7)
        -2.0259e-01 * pow(logT, 8)
        +1.1977e-02 * pow(logT, 9)
        -3.0591e-04 * pow(logT, 10);
        return exp(log_out);

      } else if (Te < 2) {
        return 7.0655e-36;   
      } else if (Te > 1000) {
        return 1.1600e-32;
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
        -5.6832e+01 * pow(logT, 0)
        -1.0690e+02 * pow(logT, 1)
        +2.2106e+02 * pow(logT, 2)
        -2.3902e+02 * pow(logT, 3)
        +1.5679e+02 * pow(logT, 4)
        -6.5878e+01 * pow(logT, 5)
        +1.8075e+01 * pow(logT, 6)
        -3.2221e+00 * pow(logT, 7)
        +3.6003e-01 * pow(logT, 8)
        -2.2931e-02 * pow(logT, 9)
        +6.3593e-04 * pow(logT, 10);
        return exp(log_out);

      } else if (Te < 2) {
        return 4.0959e-34;   
      } else if (Te > 500) {
        return 7.7710e-33;
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
    registercomponentfixedfractioncarbon("fixed_fraction_carbon");

  RegisterComponent<FixedFractionRadiation<Nitrogen_adas>>
    registercomponentfixedfractionnitrogen("fixed_fraction_nitrogen");

  RegisterComponent<FixedFractionRadiation<Neon_adas>>
    registercomponentfixedfractionneon("fixed_fraction_neon");

  RegisterComponent<FixedFractionRadiation<Argon_adas>>
    registercomponentfixedfractionargon("fixed_fraction_argon");

  // RegisterComponent<FixedFractionRadiation<RyokoArgon>>
  //   registercomponentfixedfractionargon("fixed_fraction_deprecated");
}

#endif // FIXED_FRACTION_IONS_H
