#include <fmt/core.h>

#include "../include/guarded_options.hxx"

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
    auto [access, varname] = permissions->canAccess(name, Permissions::Read, region);
    if (access) {
      updateAccessRecords(*unread_variables, varname, region);
      return *options;
    }
  }
  throw BoutException(fmt::format("Do not have read permission for {}.", name));
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
  throw BoutException(fmt::format("Do not have read permission for {}.", options->str()));
}

std::map<std::string, Permissions::Regions> GuardedOptions::unreadItems() const {
  return *unread_variables;
}

std::map<std::string, Permissions::Regions> GuardedOptions::unwrittenItems() const {
  return *unwritten_variables;
}
