#pragma once

#include <memory>

#include <bout/options.hxx>

#include "permissions.hxx"

/// A wrapper class around BOUT++ Options objects. It uses access
/// right data, stored using a Permissions object, to control reading
/// from and writing to the underlying data.
class GuardedOptions {
public:
  GuardedOptions() = default;
  // GuardedOptions(GuardedOptions &other) = default;
  // GuardedOptions(GuardedOptions &&other) = default;
  
  /// Create a guarded options object which applies the specified
  /// permissions to the underlying options object. Note that the
  /// variable names used in the Permissions object must always be the
  /// full names, relative to the highest-level of the Options
  /// hierarchy.
  GuardedOptions(Options* options, Permissions* permissions)
      : options(options), permissions(permissions),
        unread_variables(std::make_shared<std::map<std::string, Permissions::Regions>>()),
        unwritten_variables(
            std::make_shared<std::map<std::string, Permissions::Regions>>()) {
    if (permissions != nullptr) {
      *unread_variables = permissions->getVariablesWithPermission(Permissions::Read);
      *unwritten_variables =
          permissions->getVariablesWithPermission(Permissions::Write, false);
    }
  }

  /// Get a subsection or value. The result will also be wrapped in a
  /// GuardedOptions object, with the same permissions as this one.
  GuardedOptions operator[](const std::string& name);
  GuardedOptions operator[](const char* name) { return (*this)[std::string(name)]; }

  const GuardedOptions operator[](const std::string& name) const;
  const GuardedOptions operator[](const char* name) const { return (*this)[std::string(name)]; }

  std::map<std::string, GuardedOptions> getChildren();
  bool isSection(const std::string& name) const { return options->isSection(name); }
  bool isSection(const char* name) const { return (*this).isSection(std::string(name)); }
  bool isSection() const { return options->isSection(); }
  bool isSet(const std::string& name) const { return options->isSet(name); }
  bool isSet(const char* name) const { return (*this).isSet((std::string(name))); }
  bool isSet() const { return options->isSet(); }

  /// Get read-only access to the underlying Options object. Throws
  /// BoutException if there is not read-permission for this object.
  const Options& get(Permissions::Regions region = Permissions::AllRegions) const;
  /// Get read-write access to the underlying Options object. Throws
  /// BoutException if there is not write-permission for this object.
  Options& getWritable(Permissions::Regions region = Permissions::AllRegions);

  /// Returns a list of variables with read-only permission but which
  /// have not been accessed using the `get()` method.
  std::map<std::string, Permissions::Regions> unreadItems() const;

  /// Returns a list of variables with read-write permission but which
  /// have not been accessed using the `getWritable()` method.
  std::map<std::string, Permissions::Regions> unwrittenItems() const;

private:
  Options* options{nullptr};
  Permissions* permissions{nullptr};
  mutable std::shared_ptr<std::map<std::string, Permissions::Regions>> unread_variables,
      unwritten_variables;

  GuardedOptions(
      Options* options, Permissions* permissions,
      std::shared_ptr<std::map<std::string, Permissions::Regions>> unread_vars,
      std::shared_ptr<std::map<std::string, Permissions::Regions>> unwritten_vars)
      : options(options), permissions(permissions), unread_variables(unread_vars),
        unwritten_variables(unwritten_vars) {}
};
