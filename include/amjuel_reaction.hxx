#pragma once
#ifndef AMJUEL_REACTION_H
#define AMJUEL_REACTION_H

#include "component.hxx"

struct AmjuelReaction : public Component {
  AmjuelReaction(std::string, Options &alloptions, Solver *) {
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
  /// Input in SI units: 
  ///     n in m^-3
  ///     T in eV
  ///
  /// Output in SI, units m^3/s
  template <size_t rows, size_t cols>
  BoutReal evaluate(const BoutReal (&coefs)[rows][cols], BoutReal T, BoutReal n) {

    // Enforce range of validity
    n = clip(n, 1e14, 1e22); // 1e8 - 1e16 cm^-3
    T = clip(T, 0.1, 1e4);
    
    BoutReal logntilde = log(n / 1e14); // Note: 1e8 cm^-3
    BoutReal logT = log(T);
    BoutReal result = 0.0;

    BoutReal logT_n = 1.0;  // log(T) ** n
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
};

#endif // AMJUEL_REACTION_H
