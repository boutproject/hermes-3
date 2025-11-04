#pragma once
#include <array>
#include <initializer_list>
#include <map>
#include <string>
#include <vector>

/// Class to store information on whether particular variables an be
/// read from and/or written to. These permissions can apply on
/// particular regions of the domain.
class Permissions {
public:
  /// Ways in which someone is allowed to access the variable, with
  /// increasing levels of rights. "ReadIfSet" indicates that
  /// variables should only be read if already set. "Final" refers to
  /// the last time anyone is allowed to write to the variable. These
  /// two concepts need to be captured to decide the order in which to
  /// execute components.
  enum PermissionTypes { None = -1, ReadIfSet, Read, Write, Final, PERMISSION_TYPES_END };

  /// The regions of the domain to which a particular permission
  /// apply. These are designed to be used as bit-flags.
  enum Regions {
    Nowhere = 0,
    Interior = 1 << 0,
    Boundaries = 1 << 1,
    AllRegions = Interior | Boundaries
  };

  const static std::map<Regions, std::string> fundamental_regions;

  /// Data type for storing the regions of a variable which have a
  /// particular level of permission. Some examples can be seen below:
  ///
  ///     AccessRights only_read_if_set = { AllRegions, Nowhere, Nowhere, Nowhere },
  ///                  read_only = { Nowhere, AllRegions, Nowhere, Nowhere },
  ///                  write_boundaries = { Nowhere, Nowhere, Boundaries, Nowhere },
  ///                  read_and_write_everywhere = { Nowhere, AllReginos, AllRegions,
  ///                  Nowhere
  ///                  }, final_write_boundaries_read_interior = { Interior, Nowhere,
  ///                      Boundaries };
  ///
  using AccessRights = std::array<Regions, PERMISSION_TYPES_END>;

  Permissions() = default;
  Permissions(Permissions& other) = default;
  Permissions(Permissions&& other) = default;

  /// Create permission from an initialiser list. Each item in the
  /// initialiser list should be a pair made up of the name of a
  /// variable stored in an Options object (with a colon separating
  /// section names) and an array describing which regions of the
  /// domain each level of permission is applied to. The order of the
  /// permissions in the array is: read, write, final write. Note that
  /// higher permissions are taken to imply all lower permissions. See
  /// the examples below.
  ///
  ///     Permissions example({
  ///         // Permission to read charge only if it has been set
  ///         {"species:he:charge", {Permissions::AllRegions, Permissions::Nowhere,
  ///         Permissions::Nowhere,, Permissions::Nowhere}},
  ///         // Read permission for atomic mass
  ///         {"species:he:AA", {Permissions::Nowhere, Permissions::AllRegions,
  ///         Permissions::Nowhere, Permissions::Nowhere}},
  ///         // Read permissions for density
  ///         {"species:he:density", {Permissions::Nowhere, Permissions::AllRegions,
  ///         Permissions::Nowhere, Permissions::Nowhere}},
  ///         // Read and write permissions for pressure in the interior region
  ///         {"species:he:pressure", {Permissions::Nowhere, Permissions::Nowhere,
  ///         Permissions::Interior, Permissions::Nowhere}},
  ///         // Set the final value for collision frequency
  ///         {"species:he:collision_frequency", {Permissions::Nowhere,
  ///         Permissions::Nowhere, Permissions::Nowhere, Permissions::AllRegions}}
  ///     });
  ///
  /// If a variable is not included in the initialiser list then it is
  /// assumed there are no access rights. If a section name appears in
  /// the list then those permissions apply to all children of that
  /// section. The section name must end in a colon (e.g.,
  /// "species:he:"). If an additional entry is present for somethign
  /// located in that section, then the more specific entry takes
  /// precidences.
  ///
  /// Placeholders can be used in variable names by surrounding a
  /// label with curly braces. Multiple values can then be substituted
  /// for this placeholder using the Permissions::substitute
  /// method. For example, if you wanted to read the collision
  /// frequency for every species you would write:
  ///
  ///     Permissions example2({
  ///         {"species:{name}:collision_frequency", {Permissions::Nowhere,
  ///         Permissions::AllRegions, Permissions::Nowhere, Permissions::Nowhere}}
  ///     });
  ///     example2.substitute("name", {"he+", "d+", "e", "d", "he"});
  ///
  Permissions(std::initializer_list<std::pair<std::string, AccessRights>> data);

  /// Set the level of access for the various regions of the
  /// variable. This uses the same logic as the constructor. For
  /// example, to indicate that the density of helium is readable
  /// everywhere but only writeable in the interior, you would use
  ///
  ///     permissions.setAccess("species:he:density",
  ///                           {Permissions::Nowhere, Permissions::AllRegions,
  ///                           Permissions::Interior, Permissions::Nowhere})
  /// 0
  /// or, equivalently,
  ///
  ///     permissions.setAccess("species:he:density",
  ///                           {Permissions::Nowhere, Permissions::Boundary,
  ///                           Permissions::Interior, Permissions::Nowhere});
  ///
  /// As in the constructor, if the variable name is just a section in
  /// an Options object then the permissions apply to all children of
  /// that section. Placeholder names can also be used.
  void setAccess(const std::string& variable, const AccessRights& rights);

  void setAccess(const std::pair<std::string, AccessRights>& info) {
    setAccess(info.first, info.second);
  }

  /// Replace a placeholder in the names of variables stored in this
  /// object. This is useful if you need to access the same variable
  /// for multiple species. For example, the following code gives
  /// permission to read the density and write the collision frequency
  /// for every species.
  ///
  ///     Permissions example({
  ///         {"species:{name}:density", {Permissions::Nowhere, Permissions::AllRegions,
  ///         Permissions::Nowhere, Permissions::Nowhere}},
  ///         {"species:{name}:collision_frequency", {Permissions::Nowhere,
  ///         Permissions::Nowhere, Permissions::AllRegions, Permissions::Nowhere}},
  ///     });
  ///     example.substitute("name", {"d", "d+", "t", "t+", "he", "he+", "c", "c+", "e"});
  ///
  void substitute(const std::string& label,
                  const std::vector<std::string>& substitutions);

  /// Check whether users are allowed to access this variable to the
  /// given permission level, in the given region. The second item
  /// returned indicates the name of the variable or section from
  /// which the access rights are derived. If there is no matching
  /// section then it will be an empty string.
  std::pair<bool, std::string> canAccess(const std::string& variable,
                                         PermissionTypes permission = Read,
                                         Regions region = AllRegions) const;

  /// Get the highest permission level with which the given variable
  /// can be accessed in the given region. The second item
  /// returned indicates the name of the variable or section from
  /// which the access rights are derived. If there is no matching
  /// section then it will be an empty string.
  std::pair<PermissionTypes, std::string>
  getHighestPermission(const std::string& variable, Regions region = AllRegions) const;

  /// Get a set of variables and regions for which there is the
  /// specified level of permission to access. If ``highestOnly`` is
  /// true then it will only include variables/regions for which this
  /// is the highest permission.
  ///
  ///     Permissions example({"test", {Permissions::Nowhere, Permissions::AllRegions,
  ///     Permissions::AllRegions, Permissions:Nowhere}});
  ///     // Print variables which can be read
  ///     for (const auto [varname, region] :
  ///     example.getVariablesWithPermission(Permissions::Read, false))
  ///         std::cout << "Variable name: " << varname << ", Region ID: " << region <<
  ///         "\n";
  ///     // Print variables which can only be read (not written)
  ///     for (const auto [varname, region] :
  ///     example.getVariablesWithPermission(Permissions::Read, true))
  ///         std::cout << "Variable name: " << varname << ", Region ID: " << region <<
  ///         "\n";
  ///
  /// The above code would write a line of output from the first
  /// for-loop, but not the second.
  std::map<std::string, Regions>
  getVariablesWithPermission(PermissionTypes permission, bool highestOnly = true) const;

  /// Return a string version of the region names
  static std::string regionNames(const Regions regions);

private:
  /// Returns the access rights for the most specific entry in this
  /// object which matches the variable name. If there are no matching
  /// entries then the result will indicate no access rights. The
  /// string indicates the name of the variable from which the access
  /// rights were derived. It will be empty if there are no matching
  /// entries.
  std::pair<std::string, AccessRights> bestMatchRights(const std::string& variable) const;

  /// Return a set of access rights where the lower permissions have
  /// been updated so that they reflect higher permissions (e.g., read
  /// permission will be set in all cases where write permission was
  /// set).
  static AccessRights applyLowerPermissions(const AccessRights& rights);

  std::map<std::string, AccessRights> variable_permissions;
};

/// Convenience function to return an object expressing that the
/// variable should have ReadIfSet permissions in the specified regions.
std::pair<std::string, Permissions::AccessRights>
readIfSet(std::string varname, Permissions::Regions region = Permissions::AllRegions);

/// Convenience function to return an object expressing that the
/// variable should have Read permissions in the specified regions.
std::pair<std::string, Permissions::AccessRights>
readOnly(std::string varname, Permissions::Regions region = Permissions::AllRegions);

/// Convenience function to return an object expressing that the
/// variable should have Write permissions in the specified regions.
std::pair<std::string, Permissions::AccessRights>
readWrite(std::string varname, Permissions::Regions region = Permissions::AllRegions);

/// Convenience function to return an object expressing that the
/// variable should have Final permissions in the specified regions.
std::pair<std::string, Permissions::AccessRights>
writeFinal(std::string varname, Permissions::Regions region = Permissions::AllRegions);

/// Convenience function to return an object expressing that the
/// variable should have Write permissions on the boundaries. It will
/// have Read permissions in the interior, as this is normally
/// required to set the boundaries correctly.
std::pair<std::string, Permissions::AccessRights> writeBoundary(std::string varname);
