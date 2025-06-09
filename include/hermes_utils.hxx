#pragma once
#ifndef HERMES_UTILS_H
#define HERMES_UTILS_H

#include "bout/traits.hxx"
#include "bout/bout_enum_class.hxx"

inline BoutReal floor(BoutReal value, BoutReal min) {
  if (value < min)
    return min;
  return value;
}

/// Apply a smoothly varying "soft" floor to the value
/// The intention is to keep the RHS function differentiable
///
/// Note: This function cannot be used with min = 0!
///
/// Uses value + min * exp(-value / min)
///
/// Alternatives include:
///  - sqrt(value^2 + min^2)
///
inline BoutReal softFloor(BoutReal value, BoutReal min) {
  if (value < 0.0)
    return min;
  return value + min * exp(-value / min);
}

/// Apply a soft floor value \p f to a field \p var. Any value lower than
/// the floor is set to the floor.
///
/// @param[in] var  Variable to apply floor to
/// @param[in] f    The floor value. Must be > 0 (NOT zero)
/// @param[in] rgn  The region to calculate the result over
template <typename T, typename = bout::utils::EnableIfField<T>>
inline T softFloor(const T& var, BoutReal f, const std::string& rgn = "RGN_ALL") {
  checkData(var);
  T result {emptyFrom(var)};
  result.allocate();

  BOUT_FOR(d, var.getRegion(rgn)) {
    result[d] = softFloor(var[d], f);
  }

  return result;
}

template<typename T, typename = bout::utils::EnableIfField<T>>
inline T clamp(const T& var, BoutReal lo, BoutReal hi, const std::string& rgn = "RGN_ALL") {
  checkData(var);
  T result = copy(var);

  BOUT_FOR(d, var.getRegion(rgn)) {
    if (result[d] < lo) {
      result[d] = lo;
    } else if (result[d] > hi) {
      result[d] = hi;
    }
  }

  return result;
}

/// Enum that identifies the type of a species: electron, ion, neutral
BOUT_ENUM_CLASS(SpeciesType, electron, ion, neutral);

/// Identify species name string as electron, ion or neutral
inline SpeciesType identifySpeciesType(const std::string& species) {
  if (species == "e") {
    return SpeciesType::electron;
  } else if ((species == "i") or
             species.find(std::string("+")) != std::string::npos) {
    return SpeciesType::ion;
  }
  // Not electron or ion -> neutral
  return SpeciesType::neutral;
}

template<typename T, typename = bout::utils::EnableIfField<T>>
Ind3D indexAt(const T& f, int x, int y, int z) {
  int ny = f.getNy();
  int nz = f.getNz();
  return Ind3D{(x * ny + y) * nz + z, ny, nz};
}

#endif // HERMES_UTILS_H
