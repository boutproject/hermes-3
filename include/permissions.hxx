#pragma once
#include <array>
#include <cstddef>
#include <initializer_list>
#include <istream>
#include <map>
#include <ostream>
#include <regex>
#include <string>
#include <utility>
#include <vector>

#include <fmt/base.h>
#include <fmt/format.h>
#include <fmt/ranges.h>

/// Ways in which someone is allowed to access the variable, with
/// increasing levels of rights. "ReadIfSet" indicates that
/// variables should only be read if already set. "Final" refers to
/// the last time anyone is allowed to write to the variable. These
/// two concepts need to be captured to decide the order in which to
/// execute components.
enum class PermissionTypes { None = -1, ReadIfSet, Read, Write, Final, END };

/// \defgroup RegionsGroup
/// @{

/// The regions of the domain to which a particular permission
/// apply. These are designed to be used as bit-flags.
enum class Regions {
  Nowhere = 0,
  Interior = 1 << 0,
  Boundaries = 1 << 1,
  All = Interior | Boundaries
};

constexpr Regions operator&(Regions a, Regions b) {
  return static_cast<Regions>(static_cast<int>(a) & static_cast<int>(b));
}

constexpr Regions operator|(Regions a, Regions b) {
  return static_cast<Regions>(static_cast<int>(a) | static_cast<int>(b));
}

constexpr Regions operator~(Regions a) {
  return static_cast<Regions>(~static_cast<int>(a));
}

/// @}

inline auto format_as(Regions r) { return fmt::underlying(r); }

/// Class to store information on whether particular variables an be
/// read from and/or written to. These permissions can apply on
/// particular regions of the domain.
class Permissions {
public:
  Permissions() = default;
  Permissions(Permissions& x) = default;
  Permissions(Permissions&& x) = default;
  Permissions& operator=(Permissions& x) = default;
  Permissions& operator=(Permissions&& x) = default;

  const static std::map<Regions, std::string> fundamental_regions;

  /// Data type for storing the regions of a variable which have a
  /// particular level of permission. Some examples can be seen below:
  ///
  ///     AccessRights only_read_if_set = { Regions::All, Regions::Nowhere,
  ///                      Regions::Nowhere, Regions::Nowhere };
  ///     AccessRights read_only = { Regions::Nowhere, Regions::All, Regions::Nowhere,
  ///                      Regions::Nowhere };
  ///     AccessRights write_boundaries = { Regions::Nowhere, Regions::Nowhere,
  ///                      Regions::Boundaries, Regions::Nowhere };
  ///     AccessRights read_and_write_everywhere = { Regions::Nowhere,
  ///                      Regions::AllRegions, Regions::All, Regions::Nowhere };
  ///     AccessRights final_write_boundaries_read_interior = { Regions::Nowhere,
  ///                      Regions::Interior, Regions::Nowhere, Regions::Boundaries };
  ///
  using AccessRights = std::array<Regions, static_cast<size_t>(PermissionTypes::END)>;

  /// Data used to specify what access rights apply to the named variable.
  struct VarRights {
    std::string name;
    AccessRights rights;
  };

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
  ///         {"species:he:charge", {Regions::All, Regions::Nowhere,
  ///         Regions::Nowhere, Regions::Nowhere}},
  ///         // Read permission for atomic mass
  ///         {"species:he:AA", {Regions::Nowhere, Regions::All,
  ///         Regions::Nowhere, Regions::Nowhere}},
  ///         // Read permissions for density
  ///         {"species:he:density", {Regions::Nowhere, Regions::All,
  ///         Regions::Nowhere, Regions::Nowhere}},
  ///         // Read and write permissions for pressure in the interior region
  ///         {"species:he:pressure", {Regions::Nowhere, Regions::Nowhere,
  ///         Regions::Interior, Regions::Nowhere}},
  ///         // Set the final value for collision frequency
  ///         {"species:he:collision_frequency", {Regions::Nowhere,
  ///         Regions::Nowhere, Regions::Nowhere, Regions::All}}
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
  ///         {"species:{name}:collision_frequency", {Regions::Nowhere,
  ///         Regions::All, Regions::Nowhere, Regions::Nowhere}}
  ///     });
  ///     example2.substitute("name", {"he+", "d+", "e", "d", "he"});
  ///
  Permissions(std::initializer_list<VarRights> data);

  /// Set the level of access for the various regions of the
  /// variable. This uses the same logic as the constructor. For
  /// example, to indicate that the density of helium is readable
  /// everywhere but only writeable in the interior, you would use
  ///
  ///     permissions.setAccess("species:he:density",
  ///                           {Regions::Nowhere, Regions::All,
  ///                           Regions::Interior, Regions::Nowhere})
  /// 
  /// or, equivalently,
  ///
  ///     permissions.setAccess("species:he:density",
  ///                           {Regions::Nowhere, Permissions::Boundary,
  ///                           Regions::Interior, Regions::Nowhere});
  ///
  /// As in the constructor, if the variable name is just a section in
  /// an Options object then the permissions apply to all children of
  /// that section. Placeholder names can also be used.
  void setAccess(const std::string& variable, const AccessRights& rights);

  void setAccess(const VarRights& info) { setAccess(info.name, info.rights); }

  /// Replace a placeholder in the names of variables stored in this
  /// object. This is useful if you need to access the same variable
  /// for multiple species. For example, the following code gives
  /// permission to read the density and write the collision frequency
  /// for every species.
  ///
  ///     Permissions example({
  ///         {"species:{name}:density", {Regions::Nowhere, Regions::All,
  ///         Regions::Nowhere, Regions::Nowhere}},
  ///         {"species:{name}:collision_frequency", {Regions::Nowhere,
  ///         Regions::Nowhere, Regions::All, Regions::Nowhere}},
  ///     });
  ///     example.substitute("name", {"d", "d+", "t", "t+", "he", "he+", "c", "c+", "e"});
  ///
  /// Note that variable names which already have a permission set
  /// will not be overwritten.
  ///
  void substitute(const std::string& label,
                  const std::vector<std::string>& substitutions);

  /// Check there are no remaining placeholders that have not been
  /// substituted in any variables names. If there are, then throw an
  /// exception.
  void checkNoRemainingSubstitutions() const;

  /// Check whether users are allowed to access this variable to the
  /// given permission level, in the given region. The second item
  /// returned indicates the name of the variable or section from
  /// which the access rights are derived. If there is no matching
  /// section then it will be an empty string.
  std::pair<bool, std::string>
  canAccess(const std::string& variable,
            PermissionTypes permission = PermissionTypes::Read,
            Regions region = Regions::All) const;

  /// Get the highest permission level with which the given variable
  /// can be accessed in the given region. The second item returned
  /// indicates the name of the variable or section from which the
  /// access rights are derived. If there is no matching section then
  /// it will return `PermissionTypes::None` and an empty string.
  std::pair<PermissionTypes, std::string>
  getHighestPermission(const std::string& variable, Regions region = Regions::All) const;

  /// Get a set of variables and regions for which there is the
  /// specified level of permission to access. If ``highestOnly`` is
  /// true then it will only include variables/regions for which this
  /// is the highest permission.
  ///
  ///     Permissions example({"test", {Regions::Nowhere, Regions::All,
  ///                                   Regions::All, Permissions:Nowhere}});
  ///     // Print variables which can be read
  ///     for (const auto [varname, region] :
  ///     example.getVariablesWithPermission(PermissionTypes::Read), false))
  ///         std::cout << "Variable name: " << varname << ", Region ID: " << region <<
  ///         "\n";
  ///     // Print variables which can only be read (not written)
  ///     for (const auto [varname, region] :
  ///     example.getVariablesWithPermission(PermissionTypes::Read), true))
  ///         std::cout << "Variable name: " << varname << ", Region ID: " << region <<
  ///         "\n";
  ///
  /// The above code would write a line of output from the first
  /// for-loop, but not the second.
  std::map<std::string, Regions>
  getVariablesWithPermission(PermissionTypes permission, bool highestOnly = true) const;

  /// Return a string version of the region names
  static std::string regionNames(const Regions regions);

  // Allow to be streamed, so can be stored in Options.
  friend fmt::formatter<Permissions>;
  friend std::ostream& operator<<(std::ostream& os, const Permissions& permissions);
  friend std::istream& operator>>(std::istream& is, Permissions& permissions);

  /// Returns the access rights for the most specific entry in this
  /// object which matches the variable name. If there are no matching
  /// entries then the result will indicate no access rights. The
  /// string indicates the name of the variable from which the access
  /// rights were derived. It will be empty if there are no matching
  /// entries.
  VarRights bestMatchRights(const std::string& variable) const;

private:
  std::map<std::string, AccessRights> variable_permissions;

  static const std::regex LABEL_RE;
};

/// Format `Permissions` to string. Format string specification is the
/// same as when formatting a `std::map<std::string, Permissions::AccessRights>`.
/// See https://fmt.dev/12.0/syntax/#range-format-specifications.
template <>
struct fmt::formatter<Permissions>
    : formatter<std::map<std::string, Permissions::AccessRights>> {
  auto format(const Permissions& p, format_context& ctx) const
      -> format_context::iterator;
};

/// \defgroup PermissionFactories
/// @{

/// Convenience function to return an object expressing that the
/// variable should have ReadIfSet permissions in the specified regions.
inline Permissions::VarRights readIfSet(std::string varname,
                                        Regions region = Regions::All) {
  return {std::move(varname),
          {region, Regions::Nowhere, Regions::Nowhere, Regions::Nowhere}};
}

/// Convenience function to return an object expressing that the
/// variable should have Read permissions in the specified regions.
inline Permissions::VarRights readOnly(std::string varname,
                                       Regions region = Regions::All) {
  return {std::move(varname),
          {Regions::Nowhere, region, Regions::Nowhere, Regions::Nowhere}};
}

/// Convenience function to return an object expressing that the
/// variable should have Write permissions in the specified regions.
inline Permissions::VarRights readWrite(std::string varname,
                                        Regions region = Regions::All) {
  return {std::move(varname),
          {Regions::Nowhere, Regions::Nowhere, region, Regions::Nowhere}};
}

/// Convenience function to return an object expressing that the
/// variable should have Final permissions in the specified regions.
inline Permissions::VarRights writeFinal(std::string varname,
                                         Regions region = Regions::All) {
  return {std::move(varname),
          {Regions::Nowhere, Regions::Nowhere, Regions::Nowhere, region}};
}

/// Convenience function to return an object expressing that the
/// variable should have Final write permissions on the boundaries. It
/// will have Read permissions in the interior, as this is normally
/// required to set the boundaries correctly.
inline Permissions::VarRights writeBoundary(std::string varname) {
  return {std::move(varname),
          {Regions::Nowhere, Regions::Interior, Regions::Nowhere, Regions::Boundaries}};
}

/// Convenience function to return an object expressing that, if the
/// interior has been set, then the variable should have Final write
/// permissions on the boundaries and read permissions for the
/// interior.
///
/// FIXME: Currently these permissiosn are not expressed properly, due
/// to limitations in how the permission system. The boundary will
/// have write permission regardless of whether or not the interior is
/// set.
inline Permissions::VarRights writeBoundaryIfSet(std::string varname) {
  return {std::move(varname),
          {Regions::Interior, Regions::Nowhere, Regions::Nowhere, Regions::Boundaries}};
}

/// Convenience function to return an object expressing that the
/// interior has read-if-set permissions and the boundary has write
/// permissions. This differs from what `writeBoundaryIfSet`
/// is supposed to do because this will write the boundary
/// unconditionally, regardless of whether the interior is set.
inline Permissions::VarRights writeBoundaryReadInteriorIfSet(std::string varname) {
  return {std::move(varname),
          {Regions::Interior, Regions::Nowhere, Regions::Nowhere, Regions::Boundaries}};
}
/// @}

// FIXME: Ideally there would be some way to express write permissions only if set
// FIXME: Ideally we could express to write a boundary only if the interior is set

/// Write a string-representation of the permissions to an IO
/// stream. This is useful for allowing permission data to be stored
/// in Options objects.
std::ostream& operator<<(std::ostream& os, const Permissions& permissions);

/// Read a string representation of permissions from an IO
/// stream. This is useful for allowing permissions data to be
/// retrieved from Options objects. Note that there is undefined
/// behaviour if the input is corrupted; an exception may be thrown or
/// the permissions that are read may be incomplete.
std::istream& operator>>(std::istream& is, Permissions& permissions);

std::string toString(const Permissions& value);
