#pragma once
#ifndef HERMES_UTILS_H
#define HERMES_UTILS_H

#include "bout/bout_enum_class.hxx"

inline BoutReal floor(BoutReal value, BoutReal min) {
  if (value < min)
    return min;
  return value;
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

/// Enum that identifies the type of a species: electtron, ion, neutral
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
