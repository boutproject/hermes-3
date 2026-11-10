#pragma once
#ifndef HERMES_UTILS_H
#define HERMES_UTILS_H

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "bout/assert.hxx"
#include "bout/field3d.hxx"
#include "bout/traits.hxx"

inline BoutReal floor(BoutReal value, BoutReal min) { return std::max(value, min); }

/// Apply a smoothly varying "soft" floor to the value
/// The intention is to keep the RHS function differentiable
///
/// Note: This function cannot be used with min = 0!
inline BoutReal softFloor(BoutReal value, BoutReal min) {
  value = std::max(value, 0.0);
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
  T result{emptyFrom(var)};
  result.allocate();

  BOUT_FOR(d, var.getRegion(rgn)) { result[d] = softFloor(var[d], f); }

  if constexpr (std::is_base_of_v<Field3DParallel, T>) {
    if (var.isFci()) {
      for (size_t i = 0; i < var.numberParallelSlices(); ++i) {
        result.yup(i) =
            softFloor(var.yup(i), f,
                      fmt::format(FMT_STRING("RGN_YPAR_{:+d}"), static_cast<int>(i + 1)));
        result.ydown(i) = softFloor(
            var.ydown(i), f,
            fmt::format(FMT_STRING("RGN_YPAR_{:+d}"), -static_cast<int>(i + 1)));
      }
    }
  }

  return result;
}

template <typename T, typename = bout::utils::EnableIfField<T>>
inline T clamp(const T& var, BoutReal lo, BoutReal hi,
               const std::string& rgn = "RGN_ALL") {
  checkData(var);
  T result = copy(var);

  BOUT_FOR(d, var.getRegion(rgn)) {
    if (result[d] < lo) {
      result[d] = lo;
    } else if (result[d] > hi) {
      result[d] = hi;
    }
  }

  if constexpr (std::is_base_of_v<Field3DParallel, T>) {
    if (var.isFci()) {
      for (size_t i = 0; i < var.numberParallelSlices(); ++i) {
        result.yup(i) =
            clamp(var.yup(i), lo, hi,
                  fmt::format(FMT_STRING("RGN_YPAR_{:+d}"), static_cast<int>(i + 1)));
        result.ydown(i) =
            clamp(var.ydown(i), lo, hi,
                  fmt::format(FMT_STRING("RGN_YPAR_{:+d}"), -static_cast<int>(i + 1)));
      }
    }
  }

  return result;
}

/**
 * @brief Map `x` to a smooth sigmoid between `lo` and `hi`.
 * The transition is centred on (`lo` +`hi`) / 2.
 *
 * @param x The value to map
 * @param lo Lower bound, approached smoothly as x -> -inf
 * @param hi Upper bound, approached smoothly as x -> +inf
 * @param width Positive transition scale; smaller is closer to a hard step
 */
inline BoutReal sigmoid_map(BoutReal x, BoutReal lo, BoutReal hi, BoutReal width) {
  ASSERT2(width > 0.0);
  return lo + (hi - lo) / (1.0 + exp(-(x - 0.5 * (lo + hi)) / width));
}

/**
 * @brief Smoothly transition from 0 (ratio <= 1) to 1 (ratio >= 1 + width),
 * with continuous first and second derivatives.
 *
 * @param ratio A value expected to vary around 1 (e.g. density / density_threshold)
 * @param width Positive transition scale (above ratio = 1)
 */
inline BoutReal smooth_gate(BoutReal ratio, BoutReal width) {
  ASSERT2(width > 0.0);
  const BoutReal x = std::clamp((ratio - 1.0) / width, 0.0, 1.0);
  return x * x * x * (10.0 + x * (-15.0 + 6.0 * x));
}

/**
 * @brief Util function to get the keys of a std::string => T map
 *
 * @tparam T the type of the map values
 * @param map the map
 * @return std::vector<std::string> vector of keys
 */
template <typename T>
static inline std::vector<std::string> str_keys(const std::map<std::string, T>& map) {
  std::vector<std::string> keys;
  std::transform(map.begin(), map.end(), std::back_inserter(keys),
                 [](const std::pair<std::string, T>& pair) { return pair.first; });
  return keys;
}

template <typename T, typename = bout::utils::EnableIfField<T>>
Ind3D indexAt(const T& f, int x, int y, int z) {
  int ny = f.getNy();
  int nz = f.getNz();
  return Ind3D{(x * ny + y) * nz + z, ny, nz};
}

#endif // HERMES_UTILS_H

/// Function which returns true if any of a list of subtstrings is contained within string
inline bool containsAnySubstring(const std::string& mainString,
                                 const std::vector<std::string>& substrings) {
  for (const auto& subString : substrings) {
    if (mainString.find(subString) != std::string::npos) {
      return true; // Found at least one substring
    }
  }
  return false; // None of the substrings found
}

/// Takes a string representing a collision, e.g. d+_e_coll
/// Splits it using underscores and finds species, e.g. d+ and e
/// Does partial match against these, e.g. True if species1 = d+ or + and species2 = e
/// Used across all processes that require collisions for identifying the right ones.
inline bool collisionSpeciesMatch(std::string input, const std::string& species1,
                                  const std::string& species2,
                                  const std::string& reaction, const std::string& mode) {
  // Split the input string into substrings using underscore as delimiter
  std::vector<std::string> substrings;
  size_t pos = 0;
  std::string token;
  while ((pos = input.find('_')) != std::string::npos) {
    token = input.substr(0, pos);
    substrings.push_back(token);
    input.erase(0, pos + 1); // Erase string up to the current underscore
  }
  substrings.push_back(input); // Add the last substring after the last underscore

  bool species1_found = false;
  bool species2_found = false;
  bool reaction_found = false;

  // output << std::string("\n############################\n");

  if (mode == "partial") {
    if (substrings[0].find(species1) != std::string::npos) {
      species1_found = true;
    }

    if (substrings[1].find(species2) != std::string::npos) {
      species2_found = true;
    }

    if (substrings[2].find(reaction) != std::string::npos) {
      reaction_found = true;
    }

  } else if (mode == "exact") {
    // output << substrings[0] << std::string(" and ") << substrings[1] << std::string("
    // -- Testing for ") << species1 << std::string(" and ") << species2;
    if (substrings[0] == species1) {
      species1_found = true;
      // output << std::string(" <-- ") << species1 << std::string(" FOUND");
    }

    if (substrings[1] == species2) {
      species2_found = true;
      // output << std::string(" <-- ") << species2 << std::string(" FOUND");
    }

    if (substrings[2] == reaction) {
      reaction_found = true;
      // output << std::string(" <-- ") << species2 << std::string(" FOUND");
    }

    // else if (mode == "species") {

    // }

    // output << std::endl;
  } else {
    throw BoutException("Collision species match mode must be 'exact' or 'partial'");
  }

  // output << std::string("############################\n");

  // Check if the first species matches species1 and the second species matches species2
  return (species1_found && species2_found && reaction_found);
}

/// Helper class used for compile-time concatenation of string
/// literals. Taken from https://stackoverflow.com/a/62823211
template <const std::string_view&... Strs>
struct join {
  // Join all strings into a single std::array of chars
  static constexpr auto impl() noexcept {
    constexpr std::size_t len = (Strs.size() + ... + 0);
    std::array<char, len + 1> arr{};
    size_t i = 0;
    auto append = [&i, &arr](const auto& s) mutable {
      for (auto c : s) {
        arr[i++] = c;
      }
    };
    (append(Strs), ...);
    arr[len] = 0;
    return arr;
  }
  // Give the joined string static storage
  static constexpr auto arr = impl();
  // View as std::string_view
  static constexpr std::string_view value{arr.data(), arr.size() - 1};
};

// Helper to get the value out
template <const std::string_view&... Strs>
static constexpr auto join_v = join<Strs...>::value;
