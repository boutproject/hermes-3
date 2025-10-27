#include "../include/permissions.hxx"

Permissions::Permissions(std::initializer_list<std::pair<std::string, AccessRights>> data) :
    variable_permissions() {
  for (const auto& [varname, access] : data) {
    setAccess(varname, access);
  }
}

void Permissions::setAccess(const std::string& variable, const AccessRights& rights) {
  
  variable_permissions[variable] = applyLowerPermissions(rights);
}

std::string replaceAll(const std::string& str, const std::string& from, const std::string& to) {
  std::string result = str;
  if(from.empty()) return result;
  size_t start_pos = 0;
  while((start_pos = result.find(from, start_pos)) != std::string::npos) {
    result.replace(start_pos, from.length(), to);
    start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
  }
  return result;  
}

void Permissions::substitute(const std::string& label,
                             const std::vector<std::string>& substitutions) {
  for (const auto [varname, access] : variable_permissions) {
    const std::string pattern = "{" + label + "}";
    if (varname.find(pattern)  == std::string::npos) continue;
    variable_permissions.erase(varname);
    for (const std::string& val : substitutions) {
      variable_permissions[replaceAll(varname, pattern, val)] = access;
    }
  }
}

Permissions::AccessRights Permissions::bestMatchRights(const std::string& variable) const {
  Permissions::AccessRights best_candidate = {Permissions::Nowhere, Permissions::Nowhere, Permissions::Nowhere};
  int max_len = 0;
  for (const auto& [varname, rights] : variable_permissions) {
    if (varname == variable) {
      return rights;
    }
    if (varname.back() == ':' and varname.size() > max_len
        and variable.find(varname) == 0) {
      max_len = varname.size();
      best_candidate = rights;
    }
  }
  return best_candidate;
}

bool Permissions::canAccess(const std::string& variable,
                            PermissionTypes permission, Regions region) const {
  return (bestMatchRights(variable)[permission] & region) == region;
}

Permissions::PermissionTypes Permissions::getHighestPermission(const std::string & variable, Permissions::Regions region) const {
  if (region == Nowhere) return None;
  AccessRights rights = bestMatchRights(variable);
  int i = Read;
  while (i < PERMISSION_TYPES_END and (rights[i] & region) == region) {
    i++;
  }
  return static_cast<PermissionTypes>(i - 1);
}

std::vector<std::pair<std::string, Permissions::Regions>>
Permissions::getVariablesWithPermission(PermissionTypes permission,
                                        bool highestOnly) const {
  std::vector<std::pair<std::string, Permissions::Regions>> result;
  if (highestOnly and permission < PERMISSION_TYPES_END - 1) {
    for (const auto& [varname, rights] : variable_permissions) {
      auto regions = rights[permission];
      auto perm_in_regions = static_cast<Regions>(rights[permission] & ~rights[permission+1]);
      if (perm_in_regions != Nowhere) result.emplace_back(varname, perm_in_regions);
    }
  } else {
    for (const auto& [varname, rights] : variable_permissions) {
      auto regions = rights[permission];
      if (regions != Nowhere) result.emplace_back(varname, regions);
    }
  }
  return result;
}

Permissions::AccessRights Permissions::applyLowerPermissions(const AccessRights& rights) {
  AccessRights result(rights);
  for (int i = Read; i < PERMISSION_TYPES_END; i++) {
    result[i] = rights[i];
    // Higher permissions imply lower permissions
    for (int j = Read; j < i; j++) {
      result[j] = static_cast<Regions>(result[j] | rights[i]);
    }
  }
  return result;
}
