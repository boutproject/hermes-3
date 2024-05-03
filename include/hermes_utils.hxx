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

template<typename T, typename = bout::utils::EnableIfField<T>>
Ind3D indexAt(const T& f, int x, int y, int z) {
  int ny = f.getNy();
  int nz = f.getNz();
  return Ind3D{(x * ny + y) * nz + z, ny, nz};
}

#endif // HERMES_UTILS_H

/// Function which returns true if any of a list of subtstrings is contained within string
inline bool containsAnySubstring(const std::string& mainString, const std::vector<std::string>& substrings) {
  for (const auto& subString : substrings) {
      if (mainString.find(subString) != std::string::npos) {
          return true;  // Found at least one substring
      }
  }
  return false;  // None of the substrings found
}

