#include "../include/permissions.hxx"

const std::map<Permissions::Regions, std::string> Permissions::fundamental_regions = {
    {Permissions::Interior, "Interior"}, {Permissions::Boundaries, "Boundaries"}};

Permissions::Permissions(std::initializer_list<std::pair<std::string, AccessRights>> data)
    : variable_permissions() {
  for (const auto& [varname, access] : data) {
    setAccess(varname, access);
  }
}

void Permissions::setAccess(const std::string& variable, const AccessRights& rights) {

  variable_permissions[variable] = applyLowerPermissions(rights);
}

std::string replaceAll(const std::string& str, const std::string& from,
                       const std::string& to) {
  std::string result = str;
  if (from.empty())
    return result;
  size_t start_pos = 0;
  while ((start_pos = result.find(from, start_pos)) != std::string::npos) {
    result.replace(start_pos, from.length(), to);
    start_pos +=
        to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
  }
  return result;
}

void Permissions::substitute(const std::string& label,
                             const std::vector<std::string>& substitutions) {
  for (const auto [varname, access] : variable_permissions) {
    const std::string pattern = "{" + label + "}";
    if (varname.find(pattern) == std::string::npos)
      continue;
    variable_permissions.erase(varname);
    for (const std::string& val : substitutions) {
      variable_permissions[replaceAll(varname, pattern, val)] = access;
    }
  }
}

std::pair<std::string, Permissions::AccessRights>
Permissions::bestMatchRights(const std::string& variable) const {
  Permissions::AccessRights best_candidate = {Permissions::Nowhere, Permissions::Nowhere,
                                              Permissions::Nowhere};
  std::string best_candidate_name = "";
  int max_len = 0;
  for (const auto& [varname, rights] : variable_permissions) {
    if (varname == variable) {
      return {varname, rights};
    }
    if (varname.size() > max_len and variable.find(varname + ":") == 0) {
      max_len = varname.size();
      best_candidate = rights;
      best_candidate_name = varname;
    }
  }
  return {best_candidate_name, best_candidate};
}

std::pair<bool, std::string> Permissions::canAccess(const std::string& variable,
                                                    PermissionTypes permission,
                                                    Regions region) const {
  auto [match_name, match_rights] = bestMatchRights(variable);
  if ((match_rights[permission] & region) == region) {
    return {true, match_name};
  } else {
    return {false, ""};
  }
}

std::pair<Permissions::PermissionTypes, std::string>
Permissions::getHighestPermission(const std::string& variable,
                                  Permissions::Regions region) const {
  if (region == Nowhere)
    return {None, ""};
  auto [varname, rights] = bestMatchRights(variable);
  int i = ReadIfSet;
  while (i < PERMISSION_TYPES_END and (rights[i] & region) == region) {
    i++;
  }
  return {static_cast<PermissionTypes>(i - 1), varname};
}

std::map<std::string, Permissions::Regions>
Permissions::getVariablesWithPermission(PermissionTypes permission,
                                        bool highestOnly) const {
  std::map<std::string, Permissions::Regions> result;
  if (highestOnly and permission < PERMISSION_TYPES_END - 1) {
    for (const auto& [varname, rights] : variable_permissions) {
      auto regions = rights[permission];
      auto perm_in_regions =
          static_cast<Regions>(rights[permission] & ~rights[permission + 1]);
      if (perm_in_regions != Nowhere)
        result.emplace(varname, perm_in_regions);
    }
  } else {
    for (const auto& [varname, rights] : variable_permissions) {
      auto regions = rights[permission];
      if (regions != Nowhere)
        result.emplace(varname, regions);
    }
  }
  return result;
}

Permissions::AccessRights Permissions::applyLowerPermissions(const AccessRights& rights) {
  AccessRights result(rights);
  for (int i = ReadIfSet; i < PERMISSION_TYPES_END; i++) {
    result[i] = rights[i];
    // Higher permissions imply lower permissions
    for (int j = ReadIfSet; j < i; j++) {
      result[j] = static_cast<Regions>(result[j] | rights[i]);
    }
  }
  return result;
}

std::string Permissions::regionNames(const Regions regions) {
  std::string result;
  for (auto & [region, name] : fundamental_regions) {
    if ((regions & region) == region) {
      if (result.size() > 0) result += ", ";
      result += name;
    }
  }
  return result;
}

std::pair<std::string, Permissions::AccessRights> readIfSet(std::string varname) {
  return {varname,
          {Permissions::AllRegions, Permissions::Nowhere, Permissions::Nowhere,
           Permissions::Nowhere}};
}

std::pair<std::string, Permissions::AccessRights> readOnly(std::string varname) {
  return {varname,
          {Permissions::Nowhere, Permissions::AllRegions, Permissions::Nowhere,
           Permissions::Nowhere}};
}

std::pair<std::string, Permissions::AccessRights> readWrite(std::string varname) {
  return {varname,
          {Permissions::Nowhere, Permissions::Nowhere, Permissions::AllRegions,
           Permissions::Nowhere}};
}

std::pair<std::string, Permissions::AccessRights> writeFinal(std::string varname) {
  return {varname,
          {Permissions::Nowhere, Permissions::Nowhere, Permissions::Nowhere,
           Permissions::AllRegions}};
}

std::pair<std::string, Permissions::AccessRights> writeBoundary(std::string varname) {
  return {varname,
          {Permissions::Nowhere, Permissions::Interior, Permissions::Nowhere,
           Permissions::Boundaries}};
}
