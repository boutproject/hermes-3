#include <fmt/core.h>

#include "../include/guarded_options.hxx"

GuardedOptions::GuardedOptions(Options* options, Permissions* permissions)
    : options(options), permissions(permissions),
      unread_variables(std::make_shared<std::map<std::string, Permissions::Regions>>()),
      unwritten_variables(
          std::make_shared<std::map<std::string, Permissions::Regions>>()) {
  if (permissions != nullptr) {
    *unread_variables = permissions->getVariablesWithPermission(Permissions::Read);
    // Only add variables with permission ReadIfSet to
    // unread_variables if they are already present in the options
    // object
    if (options != nullptr) {
      for (auto& [varname, region] :
           permissions->getVariablesWithPermission(Permissions::ReadIfSet)) {
        // find the last component of the full varname path
        size_t last_colon = varname.find_last_of(":");
        Options& parent = (last_colon != std::string::npos)
                              ? (*options)[varname.substr(0, last_colon)]
                              : *options;
        if (parent.isSet((last_colon != std::string::npos)
                             ? varname.substr(last_colon + 1)
                             : varname)) {
          unread_variables->insert({varname, region});
        }
      }
    }
    *unwritten_variables =
        permissions->getVariablesWithPermission(Permissions::Write, false);
  }
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
}

Options& GuardedOptions::getWritable(Permissions::Regions region) {
  if (options == nullptr)
    throw BoutException(
        "Trying to access GuardedOptions when underlying options are nullptr.");
  std::string name = options->str();
  if (permissions != nullptr) {
    auto [access, varname] = permissions->canAccess(name, Permissions::Write, region);
    if (access) {
      updateAccessRecords(*unwritten_variables, varname, region);
      return *options;
    }
  }
  throw BoutException("Do not have read permission for {}.", options->str());
}

std::map<std::string, Permissions::Regions> GuardedOptions::unreadItems() const {
  return *unread_variables;
}

std::map<std::string, Permissions::Regions> GuardedOptions::unwrittenItems() const {
  return *unwritten_variables;
}
