#include <fmt/core.h>

#include "../include/guarded_options.hxx"

/// Check whether an option is set, without creating any parent
/// Options objects in the process.
bool isSetRecursive(Options& opt, std::string varname) {
  size_t colon = varname.find(":");
  if (colon == std::string::npos) {
    return opt.isSet(varname);
  }
  std::string fragment = varname.substr(0, colon);
  if (not opt.isSection(fragment)) {
    return false;
  }
  return isSetRecursive(opt[fragment], varname.substr(colon + 1));
}

GuardedOptions::GuardedOptions(Options* options, Permissions* permissions)
    : options(options), permissions(permissions),
      unread_variables(std::make_shared<std::map<std::string, Permissions::Regions>>()),
      unwritten_variables(
          std::make_shared<std::map<std::string, Permissions::Regions>>()) {
#if CHECKLEVEL >= 1
  if (permissions != nullptr) {
    *unread_variables = permissions->getVariablesWithPermission(Permissions::Read);
    // Only add variables with permission ReadIfSet to
    // unread_variables if they are already present in the options
    // object
    if (options != nullptr) {
      for (auto& [varname, region] :
           permissions->getVariablesWithPermission(Permissions::ReadIfSet)) {
        if (isSetRecursive(*options, varname)) {
          unread_variables->insert({varname, region});
        }
      }
    }
    *unwritten_variables =
        permissions->getVariablesWithPermission(Permissions::Write, false);
  }
#endif
}

GuardedOptions GuardedOptions::operator[](const std::string& name) {
  if (options == nullptr)
    throw BoutException(
        "Trying to access GuardedOptions when underlying options are nullptr.");
  return GuardedOptions(&(*options)[name], permissions, unread_variables,
                        unwritten_variables);
}

const GuardedOptions GuardedOptions::operator[](const std::string& name) const {
  if (options == nullptr)
    throw BoutException(
        "Trying to access GuardedOptions when underlying options are nullptr.");
  return GuardedOptions(&(*options)[name], permissions, unread_variables,
                        unwritten_variables);
}

std::map<std::string, GuardedOptions> GuardedOptions::getChildren() {
  std::map<std::string, GuardedOptions> result;
  for (const auto& [varname, _] : options->getChildren()) {
    result.insert({varname, (*this)[varname]});
  }
  return result;
}

void updateAccessRecords(std::map<std::string, Permissions::Regions>& records,
                         const std::string& name, Permissions::Regions region) {
  if (records.count(name) > 0) {
    Permissions::Regions new_region =
        static_cast<Permissions::Regions>(records[name] & ~region);
    if (new_region == Permissions::Nowhere) {
      records.erase(name);
    } else {
      records[name] = new_region;
    }
  }
}

const Options& GuardedOptions::get(Permissions::Regions region) const {
  if (options == nullptr)
    throw BoutException(
        "Trying to access GuardedOptions when underlying options are nullptr.");
#if CHECKLEVEL >= 1
  std::string name = options->str();
  if (permissions != nullptr) {
    auto [permission, varname] = permissions->getHighestPermission(name, region);
    if (permission >= Permissions::ReadIfSet) {
      if (permission == Permissions::ReadIfSet && !options->isSet()) {
        throw BoutException(
            "Only have permission to read {} if it is already set, which it is not.",
            name);
      }
      updateAccessRecords(*unread_variables, varname, region);
      return *options;
    }
  }
  throw BoutException("Do not have read permission for {}.", name);
#else
  return *options;
#endif
}

Options& GuardedOptions::getWritable(Permissions::Regions region) {
  if (options == nullptr)
    throw BoutException(
        "Trying to access GuardedOptions when underlying options are nullptr.");
#if CHECKLEVEL >= 1
  std::string name = options->str();
  if (permissions != nullptr) {
    auto [access, varname] = permissions->canAccess(name, Permissions::Write, region);
    if (access) {
      updateAccessRecords(*unwritten_variables, varname, region);
      return *options;
    }
  }
  throw BoutException("Do not have write permission for {}.", options->str());
#else
  return *options;
#endif
}

std::map<std::string, Permissions::Regions> GuardedOptions::unreadItems() const {
#if CHECKLEVEL >= 1
  return *unread_variables;
#else
  throw BoutException(
      "Reading of items in GuardedOptions is not tracked when CHECKLEVEL < 1");
#endif
}

std::map<std::string, Permissions::Regions> GuardedOptions::unwrittenItems() const {
#if CHECKLEVEL >= 1
  return *unwritten_variables;
#else
  throw BoutException(
      "Reading of items in GuardedOptions is not tracked when CHECKLEVEL < 1");
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
