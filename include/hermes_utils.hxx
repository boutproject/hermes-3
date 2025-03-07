#pragma once
#ifndef HERMES_UTILS_H
#define HERMES_UTILS_H

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


/// Identify species name string as electron, ion or neutral
inline std::string identifySpeciesType(const std::string& species) {

  std::string type = "";

  if (species == "e") {
    type = "electron";
  } else if (species.find(std::string("+")) != std::string::npos) {
    type = "ion";
  } else {
    type = "neutral";
  }

  return type;
}

template<typename T, typename = bout::utils::EnableIfField<T>>
Ind3D indexAt(const T& f, int x, int y, int z) {
  int ny = f.getNy();
  int nz = f.getNz();
  return Ind3D{(x * ny + y) * nz + z, ny, nz};
}

#endif // HERMES_UTILS_H
