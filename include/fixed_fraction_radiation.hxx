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

  /// Nitrogen based cooling curve used in Lipschultz 2016
  struct LipschultzNitrogen {
    BoutReal curve(BoutReal Te) {
      if (Te > 1. and Te < 80.) {
        return 5.9e-34 * sqrt(Te - 1.0) * (80. - Te) / (1. + (3.1e-3 * SQ(Te - 1.)));
      }
      return 0.0;
    }
  };

  /// Neon based cooling curve produced by Matlab polynominal curve
  /// fitting "polyval" (Ryoko 2020 Nov)
  struct RyokoNeon {
    BoutReal curve(BoutReal Te) {
      if (Te >= 3 and Te <= 100) {
        return -2.0385e-40 * pow(Te, 5)
          + 5.4824e-38 * pow(Te, 4)
          - 5.1190E-36 * pow(Te, 3)
          + 1.7347E-34 * SQ(Te)
          -3.4151E-34 * Te
          -3.2798E-34;
      } else if (Te >=2 and Te < 3) {
        return 7e-35 * (Te - 2.) + 1e-35;
      } else if (Te >=1 and Te < 2) {
        return 1e-35 * (Te - 1.);
      }
      return 0.0;
    }
  };

  /// Argon based cooling curve produced by Matlab polynominal curve
  /// fitting "polyval" (Ryoko 2020 Nov)
  struct RyokoArgon {
    BoutReal curve(BoutReal Te) {
      if (Te >= 1.5 and Te <= 100) {
        return -4.9692e-48 * pow(Te, 10)
          + 2.8025e-45 * pow(Te, 9)
          - 6.7148e-43 * pow(Te, 8)
          + 8.8636e-41 * pow(Te, 7)
          - 6.9642e-39 * pow(Te, 6)
          + 3.2559e-37 * pow(Te, 5)
          - 8.3410e-36 * pow(Te, 4)
          + 8.6011e-35 * pow(Te, 3)
          + 1.9958e-34 * pow(Te, 2)
          + 4.9864e-34 * Te
          - 9.9412e-34;
      } else if (Te >= 1.0 and Te < 1.5) {
        return 5e-35 * (Te - 1.0);
      }
      return 0.0;
    }
  };

  /// Fit from ADAS data by Mike Kryjak 19/05/2023
  /// Not super accurate above 2000 - instead it converges on near-zero
  /// Extra care taken to avoid discontinuities or jumps, although there is a small one at 0.2eV.
  struct Argon_tau_0dot5ms{
    BoutReal curve(BoutReal Te) {
      BoutReal logT = log(Te);
      BoutReal log_out = 0;

      if (Te >= 0.2 and Te <= 2500) {
        log_out = log_out 
        -8.4367e+01 * pow(logT, 0)
        +1.1075e+01 * pow(logT, 1)
        -2.3092e+00 * pow(logT, 2)
        -1.2378e+00 * pow(logT, 3)
        +8.4987e-01 * pow(logT, 4)
        +5.6445e-02 * pow(logT, 5)
        -2.0179e-01 * pow(logT, 6)
        +7.4687e-02 * pow(logT, 7)
        -1.2541e-02 * pow(logT, 8)
        +1.0245e-03 * pow(logT, 9)
        -3.3029e-05 * pow(logT, 10);
        return exp(log_out);

      } else if (Te < 0.2) {
        return 0;    /// Already really near zero
      } else if (Te > 2500) {
        return 1.2856e-33;
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

  RegisterComponent<FixedFractionRadiation<LipschultzNitrogen>>
    registercomponentfixedfractionnitrogen("fixed_fraction_nitrogen");

  RegisterComponent<FixedFractionRadiation<RyokoNeon>>
    registercomponentfixedfractionneon("fixed_fraction_neon");

  RegisterComponent<FixedFractionRadiation<Argon_tau_0dot5ms>>
    registercomponentfixedfractionargon("fixed_fraction_argon");

  // RegisterComponent<FixedFractionRadiation<RyokoArgon>>
  //   registercomponentfixedfractionargon("fixed_fraction_deprecated");
}

#endif // FIXED_FRACTION_IONS_H
