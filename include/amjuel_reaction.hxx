#pragma once
#ifndef AMJUEL_REACTION_H
#define AMJUEL_REACTION_H

#include "component.hxx"
#include "integrate.hxx"

struct AmjuelReaction : public Component {
  AmjuelReaction(std::string name, Options& alloptions, Solver*) {
    // Get the units
    const auto& units = alloptions["units"];
    Tnorm = get<BoutReal>(units["eV"]);
    Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
    FreqNorm = 1. / get<BoutReal>(units["seconds"]);
  }

protected:
  BoutReal Tnorm, Nnorm, FreqNorm; // Normalisations

  BoutReal clip(BoutReal value, BoutReal min, BoutReal max) {
    if (value < min)
      return min;
    if (value > max)
      return max;
    return value;
  }

  /// Evaluate a double polynomial fit in n and T
  /// (page 20 of amjuel.pdf)
  ///
  ///  coefs[T][n]
  /// Input in units:
  ///     n in m^-3
  ///     T in eV
  ///
  /// Output in SI, units m^3/s, or eV m^3/s for energy loss
  template <size_t rows, size_t cols>
  BoutReal evaluate(const BoutReal (&coefs)[rows][cols], BoutReal T, BoutReal n) {

    // Enforce range of validity
    n = clip(n, 1e14, 1e22); // 1e8 - 1e16 cm^-3
    T = clip(T, 0.1, 1e4);

    BoutReal logntilde = log(n / 1e14); // Note: 1e8 cm^-3
    BoutReal logT = log(T);
    BoutReal result = 0.0;

    BoutReal logT_n = 1.0; // log(T) ** n
    for (size_t n = 0; n < cols; ++n) {
      BoutReal logn_m = 1.0; // log(ntilde) ** m
      for (size_t m = 0; m < rows; ++m) {
        result += coefs[n][m] * logn_m * logT_n;
        logn_m *= logntilde;
      }
      logT_n *= logT;
    }
    return exp(result) * 1e-6; // Note: convert cm^3 to m^3
  }

  /// Electron-driven reaction
  /// e + from_ion -> to_ion [ + e? + e?]
  ///
  /// Coefficients from Amjuel:
  ///  - rate_coefs        Double-polynomial log fit [T][n] for <σv>
  ///  - radiation_coefs   Double-polynomial log fit [T][n] for electron loss
  ///  - electron_heating          Energy added to electrons per reaction [eV]
  ///  - heavy_particle_frequency  Collisionality for the heavy particle (ion/neutral)     [s^-1]
  ///  - electron_frequency        Collisionality for the electron     [s^-1]
  ///  - reaction_rate             Density source rate in     [m^-3 s^-1]
  ///  - momentum_exchange         Momentum source rate     [kg m^-2 s^-1]
  ///  - energy_exchange           Energy source rate     [W m^-3]
  ///  - energy_loss               Electron energy loss (radiation or ionisation potential)     [W m^-3]
  ///  - rate_multiplier           User-set arbitrary rate multiplier [-]
  ///  - energy_loss               User-set arbitrary multiplier on the radiation rate associated with the reaction [-]
  template <size_t rows, size_t cols>
  void electron_reaction(Options& electron, Options& from_ion, Options& to_ion,
                         const BoutReal (&rate_coefs)[rows][cols],
                         const BoutReal (&radiation_coefs)[rows][cols],
                         BoutReal electron_heating,
                         Field3D &heavy_particle_frequency,
                         Field3D &electron_frequency,
                         Field3D &reaction_rate,
                         Field3D &momentum_exchange,
                         Field3D &energy_exchange,
                         Field3D &energy_loss,
                         BoutReal rate_multiplier,
                         BoutReal radiation_multiplier) {

    Field3D Ne = get<Field3D>(electron["density"]);
    Field3D Te = get<Field3D>(electron["temperature"]);

    Field3D N1 = get<Field3D>(from_ion["density"]);
    Field3D T1 = get<Field3D>(from_ion["temperature"]);
    Field3D V1 = get<Field3D>(from_ion["velocity"]);

    auto AA = get<BoutReal>(from_ion["AA"]);
    ASSERT1(AA == get<BoutReal>(to_ion["AA"]));

    Field3D V2 = get<Field3D>(to_ion["velocity"]);

    const BoutReal from_charge =
        from_ion.isSet("charge") ? get<BoutReal>(from_ion["charge"]) : 0.0;
    const BoutReal to_charge =
        to_ion.isSet("charge") ? get<BoutReal>(to_ion["charge"]) : 0.0;

    // Calculate reaction rate in [m^3 s^-1] using cell averaging. Optionally scale by multiplier
    reaction_rate = cellAverage(
        [&](BoutReal ne, BoutReal n1, BoutReal te) {
          return ne * n1 * evaluate(rate_coefs, te * Tnorm, ne * Nnorm) * Nnorm
                 / FreqNorm * rate_multiplier;
        },
        Ne.getRegion("RGN_NOBNDRY"))(Ne, N1, Te);

    // Same as reaction_rate but without the n1 factor: returns atom/ion collisionality in [s^-1]
    heavy_particle_frequency = cellAverage(
        [&](BoutReal ne, BoutReal n1, BoutReal te) {
          return ne * evaluate(rate_coefs, te * Tnorm, ne * Nnorm) * Nnorm
                 / FreqNorm * rate_multiplier;
        },
        Ne.getRegion("RGN_NOBNDRY"))(Ne, N1, Te);

    // Same as reaction_rate but without the n1 factor: returns electron collisionality in [s^-1]
    electron_frequency = cellAverage(
        [&](BoutReal ne, BoutReal n1, BoutReal te) {
          return n1 * evaluate(rate_coefs, te * Tnorm, ne * Nnorm) * Nnorm
                 / FreqNorm * rate_multiplier;
        },
        Ne.getRegion("RGN_NOBNDRY"))(Ne, N1, Te);

    // Add collision frequency to each species' total
    add(from_ion["collision_frequency"], heavy_particle_frequency);  // Neutral if iz, ion if rec
    add(electron["collision_frequency"], electron_frequency);

    // Particles
    // For ionisation, "from_ion" is the neutral and "to_ion" is the ion
    subtract(from_ion["density_source"], reaction_rate);
    add(to_ion["density_source"], reaction_rate);

    if (from_charge != to_charge) {
      // To ensure quasineutrality, add electron density source
      add(electron["density_source"], (to_charge - from_charge) * reaction_rate);
      if (electron.isSet("velocity")) {
        // Transfer of electron kinetic to thermal energy due to density source
        auto Ve = get<Field3D>(electron["velocity"]);
        auto Ae = get<BoutReal>(electron["AA"]);
        add(electron["energy_source"], 0.5 * Ae * (to_charge - from_charge) * reaction_rate * SQ(Ve));
      }
    }

    // Momentum
    momentum_exchange = reaction_rate * AA * V1;

    subtract(from_ion["momentum_source"], momentum_exchange);
    add(to_ion["momentum_source"], momentum_exchange);

    // Kinetic energy transfer to thermal energy
    //
    // This is a combination of three terms in the pressure
    // (thermal energy) equation:
    //
    // d/dt(3/2 p_1) = ... - F_12 v_1 + W_12 - (1/2) m R v_1^2 // From ion
    //
    // d/dt(3/2 p_2) = ... - F_21 v_2 + W_21 + (1/2) m R v_2^2 // to_ion
    //
    // where forces are F_21 = -F_12, energy transfer W_21 = -W_12,
    // and masses are equal so m_1 = m_2 = m
    //
    // Here the force F_21 = m R v_1   i.e. momentum_exchange
    //
    // As p_1 -> 0 the sum of the p_1 terms must go to zero.
    // Hence:
    //
    // W_12 = F_12 v_1 + (1/2) m R v_1^2
    //      = -R m v_1^2 + (1/2) m R v_1^2
    //      = - (1/2) m R v_1^2
    //
    // d/dt(3/2 p_2) = - m R v_1 v_2 + (1/2) m R v_1^2 + (1/2) m R v_2^2
    //               = (1/2) m R (v_1 - v_2)^2
    //

    add(to_ion["energy_source"], 0.5 * AA * reaction_rate * SQ(V1 - V2));

    // Ion thermal energy transfer
    energy_exchange = reaction_rate * (3. / 2) * T1;
    subtract(from_ion["energy_source"], energy_exchange);
    add(to_ion["energy_source"], energy_exchange);

    // Electron energy loss (radiation, ionisation potential)
    energy_loss = cellAverage(
        [&](BoutReal ne, BoutReal n1, BoutReal te) {
          return ne * n1 * evaluate(radiation_coefs, te * Tnorm, ne * Nnorm) * Nnorm
                 / (Tnorm * FreqNorm) * radiation_multiplier;
        },
        Ne.getRegion("RGN_NOBNDRY"))(Ne, N1, Te);

    // Loss is reduced by heating
    energy_loss -= (electron_heating / Tnorm) * reaction_rate;

    subtract(electron["energy_source"], energy_loss);
  }
};

#endif // AMJUEL_REACTION_H
