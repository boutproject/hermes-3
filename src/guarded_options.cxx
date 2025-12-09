#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <tuple>

#include <bout/assert.hxx>
#include <bout/boutexception.hxx>
#include <bout/options.hxx>

#include "../include/guarded_options.hxx"
#include "../include/permissions.hxx"

/// Check whether an option is set, without creating any parent
/// Options objects in the process.
bool isSetRecursive(Options& opt, const std::string& varname) {
  const size_t colon = varname.find(":");
  if (colon == std::string::npos) {
    return opt.isSet(varname);
  }
  const std::string fragment = varname.substr(0, colon);
  if (not opt.isSection(fragment)) {
    return false;
  }
  return isSetRecursive(opt[fragment], varname.substr(colon + 1));
}

GuardedOptions::GuardedOptions(Options* options, Permissions* permissions)
    : options(options), permissions(permissions),
      unread_variables(std::make_shared<std::map<std::string, Regions>>()),
      unwritten_variables(std::make_shared<std::map<std::string, Regions>>()) {
  if (options == nullptr) {
    throw BoutException("Can not construct GuardedOptions with null options pointer.");
  }
  if (permissions == nullptr) {
    throw BoutException(
        "Can not construct GuardedOptions with null permissions pointer.");
  }
#if CHECKLEVEL >= 1
  *unread_variables = permissions->getVariablesWithPermission(PermissionTypes::Read);
  // Only add variables with permission ReadIfSet to
  // unread_variables if they are already present in the options
  // object
  for (auto& [varname, region] :
       permissions->getVariablesWithPermission(PermissionTypes::ReadIfSet)) {
    if (isSetRecursive(*options, varname)) {
      unread_variables->insert({varname, region});
    }
  }
  *unwritten_variables =
      permissions->getVariablesWithPermission(PermissionTypes::Write, false);
#endif
}

std::map<std::string, GuardedOptions> GuardedOptions::getChildren() {
  std::map<std::string, GuardedOptions> result;
  for (const auto& [varname, _] : options->getChildren()) {
    result.insert({varname, (*this)[varname]});
  }
  return result;
}

void updateAccessRecords(std::map<std::string, Regions>& records, const std::string& name,
                         Regions region) {
  if (records.count(name) > 0) {
    const Regions new_region = records[name] & ~region;
    if (new_region == Regions::Nowhere) {
      records.erase(name);
    } else {
      records[name] = new_region;
    }
  }
}

const Options& GuardedOptions::get(Regions region) const {
#if CHECKLEVEL >= 1
  const std::string name = options->str();
  auto [permission, varname] = permissions->getHighestPermission(name, region);
  if (permission >= PermissionTypes::ReadIfSet) {
    if (permission == PermissionTypes::ReadIfSet && !options->isSet()) {
      throw BoutException(
          "Only have permission to read {} if it is already set, which it is not.", name);
    }
    updateAccessRecords(*unread_variables, varname, region);
    return *options;
  }
  throw BoutException("Do not have read permission for {}.", name);
#else
  return *options;
#endif
}

Options& GuardedOptions::getWritable(Regions region) {
#if CHECKLEVEL >= 1
  const std::string name = options->str();
  auto [access, varname] = permissions->canAccess(name, PermissionTypes::Write, region);
  if (access) {
    updateAccessRecords(*unwritten_variables, varname, region);
    return *options;
  }
  throw BoutException("Do not have write permission for {}.", options->str());
#else
  return *options;
#endif
}

bool GuardedOptions::operator==(const GuardedOptions& other) const {
  return std::tie(options, permissions, unread_variables, unwritten_variables)
         == std::tie(other.options, other.permissions, other.unread_variables,
                     other.unwritten_variables);
}

bool GuardedOptions::operator!=(const GuardedOptions& other) const {
  return !(*this == other);
}
