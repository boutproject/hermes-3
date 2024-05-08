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


/// Takes a string representing a collision, e.g. d+_e_coll
/// Splits it using underscores and finds species, e.g. d+ and e
/// Does partial match against these, e.g. True if species1 = d+ or + and species2 = e
/// Used across all processes that require collisions for identifying the right ones.
inline bool collisionSpeciesMatch(std::string input, const std::string& species1, const std::string& species2, const std::string& reaction, const std::string& mode) {
  // Split the input string into substrings using underscore as delimiter
  std::vector<std::string> substrings;
  size_t pos = 0;
  std::string token;
  while ((pos = input.find('_')) != std::string::npos) {
      token = input.substr(0, pos);
      substrings.push_back(token);
      input.erase(0, pos + 1);  // Erase string up to the current underscore
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
    // output << substrings[0] << std::string(" and ") << substrings[1] << std::string("      -- Testing for ") << species1 << std::string(" and ") << species2;
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

