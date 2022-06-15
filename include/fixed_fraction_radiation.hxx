#pragma once
#ifndef FIXED_FRACTION_RADIATION_H
#define FIXED_FRACTION_RADIATION_H

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
}

/// Set ion densities from electron densities
///
template <typename CoolingCurve>
struct FixedFractionRadiation : public Component {
  /// Inputs
  /// - <name>
  ///   - fraction
  FixedFractionRadiation(std::string name, Options &alloptions, Solver *UNUSED(solver)) {
    auto& options = alloptions[name];

    fraction = options["fraction"]
      .doc("Impurity ion density as fraction of electron density")
      .withDefault(0.0);

    if (options["diagnose"]
        .doc("Output radiation diagnostic?")
        .withDefault<bool>(false)) {
      bout::globals::dump.addRepeat(radiation, std::string("R") + name);
    }

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

 private:
  CoolingCurve cooling; ///< The cooling curve L(T) -> Wm^3
  BoutReal fraction; ///< Fixed fraction

  bool diagnose; ///< Output radiationdiagnostic?
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

  RegisterComponent<FixedFractionRadiation<RyokoArgon>>
    registercomponentfixedfractionargon("fixed_fraction_argon");
}

#endif // FIXED_FRACTION_IONS_H
