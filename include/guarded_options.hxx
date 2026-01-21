#pragma once

#include <map>
#include <memory>
#include <string>
#include <utility>

#include <bout/assert.hxx>
#include <bout/options.hxx>

#include "permissions.hxx"

/// A wrapper class around BOUT++ Options objects. It uses access
/// right data, stored using a Permissions object, to control reading
/// from and writing to the underlying data.
class GuardedOptions {
public:
  GuardedOptions() = delete;

  /// Create a guarded options object which applies the specified
  /// permissions to the underlying options object. Note that the
  /// variable names used in the Permissions object must always be the
  /// full names, relative to the highest-level of the Options
  /// hierarchy.
  GuardedOptions(Options* options, Permissions* permissions);

  /// Get a subsection or value. The result will also be wrapped in a
  /// GuardedOptions object, with the same permissions as this one.
  GuardedOptions operator[](const std::string& name) const {
    return GuardedOptions(&(*options)[name], permissions, unread_variables,
                          unwritten_variables);
  }
  GuardedOptions operator[](const char* name) const { return (*this)[std::string(name)]; }

  std::map<std::string, GuardedOptions> getChildren();
  bool isSection(const std::string& name) const { return options->isSection(name); }
  bool isSection(const char* name) const { return (*this).isSection(std::string(name)); }
  bool isSection() const { return options->isSection(); }
  bool isSet(const std::string& name) const { return options->isSet(name); }
  bool isSet(const char* name) const { return (*this).isSet((std::string(name))); }
  bool isSet() const { return options->isSet(); }
  std::string name() const { return options->name(); }

  /// Get read-only access to the underlying Options object. Throws
  /// BoutException if there is not read-permission for this object.
  const Options& get(Regions region = Regions::All) const;
  /// Get read-write access to the underlying Options object. Throws
  /// BoutException if there is not write-permission for this object.
  Options& getWritable(Regions region = Regions::All);

  /// Returns a list of variables with read-only permission but which
  /// have not been accessed using the `get()` method.
  std::map<std::string, Regions> unreadItems() const {
#if CHECKLEVEL >= 1
    return *unread_variables;
#else
    throw BoutException(
        "Reading of items in GuardedOptions is not tracked when CHECKLEVEL < 1");
#endif
  }

  /// Returns a list of variables with read-write permission but which
  /// have not been accessed using the `getWritable()` method.
  std::map<std::string, Regions> unwrittenItems() const {
#if CHECKLEVEL >= 1
    return *unwritten_variables;
#else
    throw BoutException(
        "Reading of items in GuardedOptions is not tracked when CHECKLEVEL < 1");
#endif
  }

  bool operator==(const GuardedOptions& other) const;
  bool operator!=(const GuardedOptions& other) const;

  PermissionTypes getHighestPermission(Regions region = Regions::All) const {
    return permissions->getHighestPermission(options->str(), region).first;
  }

private:
  Options* options;
  Permissions* permissions;
  mutable std::shared_ptr<std::map<std::string, Regions>> unread_variables,
      unwritten_variables;

  GuardedOptions(Options* options, Permissions* permissions,
                 std::shared_ptr<std::map<std::string, Regions>> unread_vars,
                 std::shared_ptr<std::map<std::string, Regions>> unwritten_vars)
      : options(std::move(options)), permissions(std::move(permissions)),
        unread_variables(std::move(unread_vars)),
        unwritten_variables(std::move(unwritten_vars)) {}
};
