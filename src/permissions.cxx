#include <regex>

#include "../include/permissions.hxx"
#include "bout/boutexception.hxx"
#include <fmt/core.h>
#include <fmt/ranges.h>

// TODO: It might be useful to add an optional condition which must be
// met for conditions to apply. So a variable is only read or written
// if another variable is already set, for example. Could potentially
// have conditions based on values as well. Could put in boolean
// logic. This would allow finer-grained access control in Components,
// but is it worth the hassle?
//
// struct Condition {
//   virtual bool eval() const = 0;
// }
//
// struct IsSetCondition : Condition {
//   std::string varname
// }
//
// struct TrueCondition;
// struct AndCondition;
// struct OrCondition;
//
// There should be ways to do this with templates to allow greater inlining
//

/// Return a set of access rights where the lower permissions have
/// been updated so that they reflect higher permissions (e.g., read
/// permission will be set in all cases where write permission was
/// set).
Permissions::AccessRights applyLowerPermissions(const Permissions::AccessRights& rights) {
  Permissions::AccessRights result(rights);
  for (int i = static_cast<int>(PermissionTypes::ReadIfSet);
       i < static_cast<int>(PermissionTypes::END); i++) {
    // Higher permissions imply lower permissions
    for (int j = static_cast<int>(PermissionTypes::ReadIfSet); j < i; j++) {
      result[j] = result[j] | rights[i];
    }
  }
  return result;
}

const std::map<Regions, std::string> Permissions::fundamental_regions = {
    {Regions::Interior, "Interior"}, {Regions::Boundaries, "Boundaries"}};

Permissions::Permissions(std::initializer_list<VarRights> data) : variable_permissions() {
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
  for (auto it = variable_permissions.begin(); it != variable_permissions.end();) {
    const auto [varname, access] = *it;
    const std::string pattern = "{" + label + "}";
    if (varname.find(pattern) == std::string::npos) {
      it++;
      continue;
    }
    it = variable_permissions.erase(it);
    for (const std::string& val : substitutions) {
      const std::string newname = replaceAll(varname, pattern, val);
      // Do not overwrite permissiosn that are already set
      if (variable_permissions.count(newname) == 0) {
        variable_permissions[newname] = access;
      }
    }
  }
}

Permissions::VarRights Permissions::bestMatchRights(const std::string& variable) const {
  auto match = variable_permissions.find(variable);
  if (match != variable_permissions.end()) {
    return {match->first, match->second};
  }
  Permissions::AccessRights best_candidate = {Regions::Nowhere, Regions::Nowhere,
                                              Regions::Nowhere};
  std::string best_candidate_name = "";
  int max_len = 0;
  for (const auto& [varname, rights] : variable_permissions) {
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
  if ((match_rights[static_cast<int>(permission)] & region) == region) {
    return {true, match_name};
  } else {
    return {false, ""};
  }
}

std::pair<PermissionTypes, std::string>
Permissions::getHighestPermission(const std::string& variable, Regions region) const {
  if (region == Regions::Nowhere)
    return {PermissionTypes::None, ""};
  auto [varname, rights] = bestMatchRights(variable);
  int i = static_cast<int>(PermissionTypes::ReadIfSet);
  while (i < static_cast<int>(PermissionTypes::END) and (rights[i] & region) == region) {
    i++;
  }
  return {static_cast<PermissionTypes>(i - 1), varname};
}

std::map<std::string, Regions>
Permissions::getVariablesWithPermission(PermissionTypes permission,
                                        bool highestOnly) const {
  if (permission == PermissionTypes::None) {
    throw BoutException("Can not return information on variables with no permission.");
  }
  std::map<std::string, Regions> result;
  if (highestOnly
      and static_cast<int>(permission) < static_cast<int>(PermissionTypes::END) - 1) {
    for (const auto& [varname, rights] : variable_permissions) {
      auto regions = rights[static_cast<int>(permission)];
      auto perm_in_regions = rights[static_cast<int>(permission)]
                             & ~rights[static_cast<int>(permission) + 1];
      if (perm_in_regions != Regions::Nowhere)
        result.emplace(varname, perm_in_regions);
    }
  } else {
    for (const auto& [varname, rights] : variable_permissions) {
      auto regions = rights[static_cast<int>(permission)];
      if (regions != Regions::Nowhere)
        result.emplace(varname, regions);
    }
  }
  return result;
}

std::string Permissions::regionNames(const Regions regions) {
  std::vector<std::string> regions_present(fundamental_regions.size());
  for (auto & [region, name] : fundamental_regions) {
    if ((regions & region) == region) {
      regions_present.push_back(name);
    }
  }
  return fmt::format("{}", fmt::join(regions_present, ", "));
}

auto fmt::formatter<Permissions>::format(const Permissions& p, format_context& ctx) const
    -> format_context::iterator {
  return formatter<std::map<std::string, Permissions::AccessRights>>::format(
      p.variable_permissions, ctx);
}

std::ostream& operator<<(std::ostream& os, const Permissions& permissions) {
  os << fmt::format("{}", permissions);
  return os;
}

std::istream& operator>>(std::istream& is, Permissions& permissions) {
  std::string tmp, object;
  std::getline(is, tmp, '{');
  if (is.eof()) {
    throw BoutException("Error parsing Permissions data; no opening bracket.");
  }
  std::getline(is, object, '}');
  if (is.eof()) {
    throw BoutException("Error parsing Permissions data; no closing bracket.");
  }
  permissions.variable_permissions.clear();

  // FIXME: This will just skip over any malformed content, without an error or warning
  std::vector<std::string> rights_patterns;
  for (int i = 0; i < static_cast<int>(PermissionTypes::END); i++) {
    rights_patterns.push_back("\\s*(\\d+)\\s*");
  }
  const std::regex re(
      fmt::format("\\s*\"([^\"]+)\"\\s*:\\s*\\[{}\\]", fmt::join(rights_patterns, ",")));
  auto items_begin = std::sregex_iterator(object.begin(), object.end(), re);
  auto items_end = std::sregex_iterator();

  for (std::sregex_iterator i = items_begin; i != items_end; ++i) {
    std::smatch item = *i;
    Permissions::AccessRights rights;
    for (int i = 0; i < static_cast<int>(PermissionTypes::END); i++) {
      std::string val = item[2 + i];
      rights[i] = static_cast<Regions>(std::stoi(val));
    }
    permissions.setAccess(item[1], rights);
  }

  return is;
}
