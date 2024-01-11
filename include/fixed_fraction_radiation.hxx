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
        -8.45410692e+01 * pow(logT, 0)
        +1.57727040e+01 * pow(logT, 1)
        -1.54264860e+01 * pow(logT, 2)
        +1.49409902e+01 * pow(logT, 3)
        -1.04815113e+01 * pow(logT, 4)
        +5.00924595e+00 * pow(logT, 5)
        -1.60029106e+00 * pow(logT, 6)
        +3.29455609e-01 * pow(logT, 7)
        -4.14036827e-02 * pow(logT, 8)
        +2.87063206e-03 * pow(logT, 9)
        -8.38888002e-05 * pow(logT, 10);
        return exp(log_out);

    } else if (Te < 1.5) {
        return 1.95353412e-35;   
    } else {
        return 1.22649600e-32;
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
        -8.21475117e+01 * pow(logT, 0)
        +1.28929854e+01 * pow(logT, 1)
        -4.74266289e+01 * pow(logT, 2)
        +7.45222324e+01 * pow(logT, 3)
        -5.75710722e+01 * pow(logT, 4)
        +2.57375965e+01 * pow(logT, 5)
        -7.12758563e+00 * pow(logT, 6)
        +1.24287546e+00 * pow(logT, 7)
        -1.32943407e-01 * pow(logT, 8)
        +7.97368445e-03 * pow(logT, 9)
        -2.05487897e-04 * pow(logT, 10);
        return exp(log_out);

    } else if (Te < 2) {
        return 6.35304113e-36;   
    } else {
        return 1.17894628e-32;
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
        -5.01649969e+01 * pow(logT, 0)
        -1.35749724e+02 * pow(logT, 1)
        +2.73509608e+02 * pow(logT, 2)
        -2.92109992e+02 * pow(logT, 3)
        +1.90120639e+02 * pow(logT, 4)
        -7.95164871e+01 * pow(logT, 5)
        +2.17762218e+01 * pow(logT, 6)
        -3.88334992e+00 * pow(logT, 7)
        +4.34730098e-01 * pow(logT, 8)
        -2.77683605e-02 * pow(logT, 9)
        +7.72720422e-04 * pow(logT, 10);
        return exp(log_out);

    } else if (Te < 2) {
        return 4.34835380e-34;   
    } else {
        return 8.11096182e-33;
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
        -7.87837896e+01 * pow(logT, 0)
        +1.55326376e+00 * pow(logT, 1)
        +1.65898194e+01 * pow(logT, 2)
        -3.23804546e+01 * pow(logT, 3)
        +3.12784663e+01 * pow(logT, 4)
        -1.74826039e+01 * pow(logT, 5)
        +5.91393245e+00 * pow(logT, 6)
        -1.22974105e+00 * pow(logT, 7)
        +1.54004499e-01 * pow(logT, 8)
        -1.06797106e-02 * pow(logT, 9)
        +3.15657594e-04 * pow(logT, 10);
        return exp(log_out);

    } else if (Te < 1) {
        return 6.00623928e-35;   
    } else {
        return 4.53057707e-33;
    }
    }
  };

  /// Argon simplified 1
  /// Based on the ADAS curve above but simplified as a linear interpolation
  /// between the LHS minimum, peak, bottom of RHS slope and final value at Te = 3000eV.
  /// RHS shoulder is preserved.
  /// Helpful for studying impact of cooling curve nonlinearity 
  struct Argon_simplified1{
    BoutReal curve(BoutReal Te) {
      
     if (Te < 0.52) { 
        return 0; 
    
     } else if (Te >= 0.52 and Te < 2.50) {
        return (1.953534e-35*(2.50 - Te) + 6.053680e-34*(Te - 0.52)) / (2.50 - 0.52);

     } else if (Te >= 2.50 and Te < 19.72) {
        return (6.053680e-34*(19.72 - Te) + 2.175237e-31*(Te - 2.50)) / (19.72 - 2.50);

     } else if (Te >= 19.72 and Te < 59.98) {
        return (2.175237e-31*(59.98 - Te) + 4.190918e-32*(Te - 19.72)) / (59.98 - 19.72);

     } else if (Te >= 59.98 and Te < 3000.00) {
        return (4.190918e-32*(3000.00 - Te) + 1.226496e-32*(Te - 59.98)) / (3000.00 - 59.98);

     } else {
        return 1.226496e-32;
     }
    }
  };

  /// Argon simplified 2
  /// Based on the ADAS curve above but simplified as a linear interpolation
  /// between the LHS minimum, peak, and the bottom of RHS slope.
  /// RHS shoulder is eliminated, radiation becomes 0 at 60eV.
  /// Helpful for studying impact of cooling curve nonlinearity 
  struct Argon_simplified2{
    BoutReal curve(BoutReal Te) {
      
     if (Te < 0.52) { 
        return 0; 
    
     } else if (Te >= 0.52 and Te < 2.50) {
        return (1.953534e-35*(2.50 - Te) + 6.053680e-34*(Te - 0.52)) / (2.50 - 0.52);

     } else if (Te >= 2.50 and Te < 19.72) {
        return (6.053680e-34*(19.72 - Te) + 2.175237e-31*(Te - 2.50)) / (19.72 - 2.50);

     } else if (Te >= 19.72 and Te < 59.98) {
        return (2.175237e-31*(59.98 - Te) + 0.000000e+00*(Te - 19.72)) / (59.98 - 19.72);

     } else {
        return 0.0;
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

    radiation_multiplier = options["R_multiplier"]
      .doc("Scale the radiation rate by this factor")
      .withDefault<BoutReal>(1.0);

    output<<std::string("\n\n****************************************************\n");
    output << name << radiation_multiplier;
    output<<std::string("\n****************************************************\n\n");

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
                              return ne * ni * cooling.curve(te * Tnorm) * radiation_multiplier * 
                              Nnorm / (SI::qe * Tnorm * FreqNorm);
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
  BoutReal radiation_multiplier; ///< Scale the radiation rate by this factor
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

  RegisterComponent<FixedFractionRadiation<Argon_simplified1>>
    registercomponentfixedfractionargonsimplified1("fixed_fraction_argon_simplified1");

  RegisterComponent<FixedFractionRadiation<Argon_simplified2>>
    registercomponentfixedfractionargonsimplified2("fixed_fraction_argon_simplified2");

}

#endif // FIXED_FRACTION_IONS_H
