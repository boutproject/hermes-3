#include <map>
#include <sstream>
#include <string>
#include <utility>

#include <bout/boutexception.hxx>
#include "gtest/gtest.h"

#include "../include/permissions.hxx"

auto make_access = std::make_pair<bool, std::string>;
auto make_permission = std::make_pair<PermissionTypes, std::string>;

TEST(PermissionsTests, TestCanAccess) {
  const Permissions example({
      readIfSet("species:he:charge"),
      readOnly("species:he:density"),
      // Read and write permissions for pressure in the interior region
      {"species:he:pressure",
       {Regions::Nowhere, Regions::Nowhere, Regions::Interior, Regions::Nowhere}},
      // Set the final value for collision frequency
      writeFinal("species:he:collision_frequency"),
      // Only allow reading of boundary velocity
      {"species:he:velocity",
       {Regions::Nowhere, Regions::Boundaries, Regions::Nowhere, Regions::Nowhere}},
      readOnly("species:d"),
      {"species:d:pressure",
       {Regions::Nowhere, Regions::Nowhere, Regions::Interior, Regions::Nowhere}},
      {"species:d:collision_frequencies",
       {Regions::Nowhere, Regions::Nowhere, Regions::Boundaries, Regions::Nowhere}},
  });

  auto no_access = make_access(false, "");

  // Check whether we have read permission for the variables across the entire domain
  EXPECT_EQ(example.canAccess("species:he:charge"), no_access);
  EXPECT_EQ(example.canAccess("species:he:density"),
            make_access(true, "species:he:density"));
  EXPECT_EQ(example.canAccess("species:he:pressure"), no_access);
  EXPECT_EQ(example.canAccess("species:he:collision_frequency"),
            make_access(true, "species:he:collision_frequency"));
  EXPECT_EQ(example.canAccess("species:he:velocity"), no_access);
  EXPECT_EQ(example.canAccess("unset"), no_access);

  // Check whether we have ReadIfSet permissions
  EXPECT_EQ(example.canAccess("species:he:charge", PermissionTypes::ReadIfSet),
            make_access(true, "species:he:charge"));
  EXPECT_EQ(example.canAccess("species:he:density", PermissionTypes::ReadIfSet),
            make_access(true, "species:he:density"));
  EXPECT_EQ(example.canAccess("unset", PermissionTypes::ReadIfSet), no_access);

  // Check whether we have write permission for the variables across the entire domain
  EXPECT_EQ(example.canAccess("species:he:charge", PermissionTypes::Write), no_access);
  EXPECT_EQ(example.canAccess("species:he:density", PermissionTypes::Write), no_access);
  EXPECT_EQ(example.canAccess("species:he:pressure", PermissionTypes::Write), no_access);
  EXPECT_EQ(example.canAccess("species:he:collision_frequency", PermissionTypes::Write),
            make_access(true, "species:he:collision_frequency"));
  EXPECT_EQ(example.canAccess("species:he:velocity", PermissionTypes::Write), no_access);
  EXPECT_EQ(example.canAccess("unset", PermissionTypes::Write), no_access);

  // Check whether we have read permission at the boundaries
  EXPECT_EQ(
      example.canAccess("species:he:charge", PermissionTypes::Read, Regions::Boundaries),
      no_access);
  EXPECT_EQ(
      example.canAccess("species:he:density", PermissionTypes::Read, Regions::Boundaries),
      make_access(true, "species:he:density"));
  EXPECT_EQ(example.canAccess("species:he:pressure", PermissionTypes::Read,
                              Regions::Boundaries),
            no_access);
  EXPECT_EQ(example.canAccess("species:he:collision_frequency", PermissionTypes::Read,
                              Regions::Boundaries),
            make_access(true, "species:he:collision_frequency"));
  EXPECT_EQ(example.canAccess("species:he:velocity", PermissionTypes::Read,
                              Regions::Boundaries),
            make_access(true, "species:he:velocity"));
  EXPECT_EQ(example.canAccess("unset", PermissionTypes::Read, Regions::Boundaries),
            no_access);

  // Check permissions set for whole sections
  EXPECT_EQ(example.canAccess("species:d:pressure"), no_access);
  EXPECT_EQ(example.canAccess("species:d:pressure", PermissionTypes::Write), no_access);
  EXPECT_EQ(
      example.canAccess("species:d:pressure", PermissionTypes::Write, Regions::Interior),
      make_access(true, "species:d:pressure"));
  EXPECT_EQ(example.canAccess("species:d:velocity"), make_access(true, "species:d"));
  EXPECT_EQ(example.canAccess("species:d:velocity", PermissionTypes::Write), no_access);
  EXPECT_EQ(
      example.canAccess("species:d:velocity", PermissionTypes::Read, Regions::Boundaries),
      make_access(true, "species:d"));
  EXPECT_EQ(example.canAccess("species:d:collision_frequencies:d_d_coll"), no_access);
  EXPECT_EQ(example.canAccess("species:d:collision_frequencies:d_d_coll",
                              PermissionTypes::Write),
            no_access);
  EXPECT_EQ(example.canAccess("species:d:collision_frequencies:d_d_coll",
                              PermissionTypes::Write, Regions::Boundaries),
            make_access(true, "species:d:collision_frequencies"));

  // Check permissions for a species that might be mistaken for one of
  // the sections we've given permissions for
  EXPECT_EQ(example.canAccess("species:d+"), no_access);
  EXPECT_EQ(example.canAccess("species:d+", PermissionTypes::Write), no_access);
  EXPECT_EQ(example.canAccess("species:d+", PermissionTypes::Read, Regions::Interior),
            no_access);
  EXPECT_EQ(example.canAccess("species:d+", PermissionTypes::Read, Regions::Boundaries),
            no_access);
}

TEST(PermissionsTests, TestGetHighestPermission) {
  const Permissions example({
      {"species:he:charge",
       {Regions::All, Regions::Nowhere, Regions::Nowhere, Regions::Nowhere}},
      {"species:he:density",
       {Regions::Nowhere, Regions::All, Regions::Nowhere, Regions::Boundaries}},
      // Read and write permissions for pressure in the interior region
      {"species:he:pressure",
       {Regions::Nowhere, Regions::Boundaries, Regions::Interior, Regions::Nowhere}},
      // Set the final value for collision frequency
      {"species:he:collision_frequency",
       {Regions::Nowhere, Regions::Interior, Regions::Nowhere, Regions::All}},
      // Only allow reading of boundary velocity
      {"species:he:velocity",
       {Regions::Nowhere, Regions::Boundaries, Regions::Nowhere, Regions::Nowhere}},
      {"species:d", {Regions::Nowhere, Regions::All, Regions::Nowhere, Regions::Nowhere}},
      {"species:d:pressure",
       {Regions::Nowhere, Regions::Nowhere, Regions::Interior, Regions::Nowhere}},
      {"species:d:collision_frequencies",
       {Regions::Nowhere, Regions::Nowhere, Regions::Boundaries, Regions::Nowhere}},
  });

  auto no_permission = make_permission(PermissionTypes::None, "");

  // Get the highest permission that covers the entire domain
  EXPECT_EQ(example.getHighestPermission("species:he:charge"),
            make_permission(PermissionTypes::ReadIfSet, "species:he:charge"));
  EXPECT_EQ(example.getHighestPermission("species:he:density"),
            make_permission(PermissionTypes::Read, "species:he:density"));
  EXPECT_EQ(example.getHighestPermission("species:he:pressure"),
            make_permission(PermissionTypes::Read, "species:he:pressure"));
  EXPECT_EQ(example.getHighestPermission("species:he:collision_frequency"),
            make_permission(PermissionTypes::Final, "species:he:collision_frequency"));
  EXPECT_EQ(example.getHighestPermission("species:he:velocity"),
            make_permission(PermissionTypes::None, "species:he:velocity"));
  EXPECT_EQ(example.getHighestPermission("species:d:pressure"),
            make_permission(PermissionTypes::None, "species:d:pressure"));
  EXPECT_EQ(example.getHighestPermission("species:d:velocity"),
            make_permission(PermissionTypes::Read, "species:d"));
  EXPECT_EQ(example.getHighestPermission("species:d:collision_frequencies:d_d_coll"),
            make_permission(PermissionTypes::None, "species:d:collision_frequencies"));
  EXPECT_EQ(example.getHighestPermission("unset"), no_permission);

  // Get the highest permission on the boundaries
  EXPECT_EQ(example.getHighestPermission("species:he:charge", Regions::Boundaries),
            make_permission(PermissionTypes::ReadIfSet, "species:he:charge"));
  EXPECT_EQ(example.getHighestPermission("species:he:density", Regions::Boundaries),
            make_permission(PermissionTypes::Final, "species:he:density"));
  EXPECT_EQ(example.getHighestPermission("species:he:pressure", Regions::Boundaries),
            make_permission(PermissionTypes::Read, "species:he:pressure"));
  EXPECT_EQ(
      example.getHighestPermission("species:he:collision_frequency", Regions::Boundaries),
      make_permission(PermissionTypes::Final, "species:he:collision_frequency"));
  EXPECT_EQ(example.getHighestPermission("species:he:velocity", Regions::Boundaries),
            make_permission(PermissionTypes::Read, "species:he:velocity"));
  EXPECT_EQ(example.getHighestPermission("species:d:pressure", Regions::Boundaries),
            make_permission(PermissionTypes::None, "species:d:pressure"));
  EXPECT_EQ(example.getHighestPermission("species:d:velocity", Regions::Boundaries),
            make_permission(PermissionTypes::Read, "species:d"));
  EXPECT_EQ(example.getHighestPermission("species:d:collision_frequencies:d_d_coll",
                                         Regions::Boundaries),
            make_permission(PermissionTypes::Write, "species:d:collision_frequencies"));
  EXPECT_EQ(example.getHighestPermission("unset", Regions::Boundaries), no_permission);

  // Get the highest permission on the interior
  EXPECT_EQ(example.getHighestPermission("species:he:charge", Regions::Interior),
            make_permission(PermissionTypes::ReadIfSet, "species:he:charge"));
  EXPECT_EQ(example.getHighestPermission("species:he:density", Regions::Interior),
            make_permission(PermissionTypes::Read, "species:he:density"));
  EXPECT_EQ(example.getHighestPermission("species:he:pressure", Regions::Interior),
            make_permission(PermissionTypes::Write, "species:he:pressure"));
  EXPECT_EQ(
      example.getHighestPermission("species:he:collision_frequency", Regions::Interior),
      make_permission(PermissionTypes::Final, "species:he:collision_frequency"));
  EXPECT_EQ(example.getHighestPermission("species:he:velocity", Regions::Interior),
            make_permission(PermissionTypes::None, "species:he:velocity"));
  EXPECT_EQ(example.getHighestPermission("species:d:pressure", Regions::Interior),
            make_permission(PermissionTypes::Write, "species:d:pressure"));
  EXPECT_EQ(example.getHighestPermission("species:d:velocity", Regions::Interior),
            make_permission(PermissionTypes::Read, "species:d"));
  EXPECT_EQ(example.getHighestPermission("species:d:collision_frequencies:d_d_coll",
                                         Regions::Interior),
            make_permission(PermissionTypes::None, "species:d:collision_frequencies"));
  EXPECT_EQ(example.getHighestPermission("unset", Regions::Interior), no_permission);

  // Check the permission for the "Nowhere" region is always "None"
  EXPECT_EQ(example.getHighestPermission("species:he:charge", Regions::Nowhere),
            no_permission);
  EXPECT_EQ(example.getHighestPermission("species:he:density", Regions::Nowhere),
            no_permission);
  EXPECT_EQ(example.getHighestPermission("species:he:pressure", Regions::Nowhere),
            no_permission);
  EXPECT_EQ(
      example.getHighestPermission("species:he:collision_frequency", Regions::Nowhere),
      no_permission);
  EXPECT_EQ(example.getHighestPermission("species:he:velocity", Regions::Nowhere),
            no_permission);
  EXPECT_EQ(example.getHighestPermission("species:d:pressure", Regions::Nowhere),
            no_permission);
  EXPECT_EQ(example.getHighestPermission("species:d:velocity", Regions::Nowhere),
            no_permission);
  EXPECT_EQ(example.getHighestPermission("species:d:collision_frequencies:d_d_coll",
                                         Regions::Nowhere),
            no_permission);
  EXPECT_EQ(example.getHighestPermission("unset", Regions::Nowhere), no_permission);

  // Check permissions for a species that might be mistaken for one of
  // the sections we've given permissions for
  EXPECT_EQ(example.getHighestPermission("species:d+"), no_permission);
  EXPECT_EQ(example.getHighestPermission("species:d+", Regions::Interior), no_permission);
  EXPECT_EQ(example.getHighestPermission("species:d+", Regions::Boundaries),
            no_permission);
}

TEST(PermissionsTests, TestSetAccess) {
  Permissions example({
      {"species:he:density",
       {Regions::Nowhere, Regions::All, Regions::Nowhere, Regions::Nowhere}},
      // Read and write permissions for pressure in the interior region
      {"species:he:pressure",
       {Regions::Nowhere, Regions::Nowhere, Regions::Interior, Regions::Nowhere}},
  });

  EXPECT_EQ(example.getHighestPermission("species:he:density"),
            make_permission(PermissionTypes::Read, "species:he:density"));
  example.setAccess("species:he:density", {Regions::Nowhere, Regions::Nowhere,
                                           Regions::Boundaries, Regions::Nowhere});
  EXPECT_EQ(example.getHighestPermission("species:he:density"),
            make_permission(PermissionTypes::None, "species:he:density"));
  EXPECT_EQ(example.getHighestPermission("species:he:density", Regions::Boundaries),
            make_permission(PermissionTypes::Write, "species:he:density"));

  EXPECT_EQ(example.getHighestPermission("species:he:pressure", Regions::Interior),
            make_permission(PermissionTypes::Write, "species:he:pressure"));
  EXPECT_EQ(example.getHighestPermission("species:he:pressure", Regions::All),
            make_permission(PermissionTypes::None, "species:he:pressure"));
  example.setAccess("species:he:pressure",
                    {Regions::Nowhere, Regions::All, Regions::Nowhere, Regions::Nowhere});
  EXPECT_EQ(example.getHighestPermission("species:he:pressure", Regions::Interior),
            make_permission(PermissionTypes::Read, "species:he:pressure"));
  EXPECT_EQ(example.getHighestPermission("species:he:pressure"),
            make_permission(PermissionTypes::Read, "species:he:pressure"));

  EXPECT_EQ(example.getHighestPermission("unset", Regions::Interior),
            make_permission(PermissionTypes::None, ""));
  EXPECT_EQ(example.getHighestPermission("unset", Regions::Boundaries),
            make_permission(PermissionTypes::None, ""));
  example.setAccess("unset", {Regions::Nowhere, Regions::Interior, Regions::Nowhere,
                              Regions::Boundaries});
  EXPECT_EQ(example.getHighestPermission("unset", Regions::All),
            make_permission(PermissionTypes::Read, "unset"));
  EXPECT_EQ(example.getHighestPermission("unset", Regions::Boundaries),
            make_permission(PermissionTypes::Final, "unset"));
}

TEST(PermissionsTests, TestGetVariablesWithPermissions) {
  const Permissions example(
      {{"species:he:density",
        {Regions::Nowhere, Regions::All, Regions::Nowhere, Regions::Boundaries}},
       // Read and write permissions for pressure in the interior region
       {"species:he:pressure",
        {Regions::Nowhere, Regions::Boundaries, Regions::Interior, Regions::Nowhere}},
       // Set the final value for collision frequency
       {"species:he:collision_frequency",
        {Regions::Nowhere, Regions::Interior, Regions::Nowhere, Regions::All}},
       // Only allow reading of boundary velocity
       {"species:he:velocity",
        {Regions::Nowhere, Regions::Boundaries, Regions::Nowhere, Regions::Nowhere}}});

  auto read_only = example.getVariablesWithPermission(PermissionTypes::Read);
  EXPECT_EQ(read_only.size(), 3);
  EXPECT_EQ(read_only["species:he:density"], Regions::Interior);
  EXPECT_EQ(read_only["species:he:pressure"], Regions::Boundaries);
  EXPECT_EQ(read_only["species:he:velocity"], Regions::Boundaries);

  auto readable = example.getVariablesWithPermission(PermissionTypes::Read, false);
  EXPECT_EQ(readable.size(), 4);
  EXPECT_EQ(readable["species:he:density"], Regions::All);
  EXPECT_EQ(readable["species:he:pressure"], Regions::All);
  EXPECT_EQ(readable["species:he:collision_frequency"], Regions::All);
  EXPECT_EQ(readable["species:he:velocity"], Regions::Boundaries);

  auto write_nonfinal = example.getVariablesWithPermission(PermissionTypes::Write, true);
  EXPECT_EQ(write_nonfinal.size(), 1);
  EXPECT_EQ(write_nonfinal["species:he:pressure"], Regions::Interior);

  auto writable = example.getVariablesWithPermission(PermissionTypes::Write, false);
  EXPECT_EQ(writable.size(), 3);
  EXPECT_EQ(writable["species:he:density"], Regions::Boundaries);
  EXPECT_EQ(writable["species:he:pressure"], Regions::Interior);
  EXPECT_EQ(writable["species:he:collision_frequency"], Regions::All);

  auto final_write = example.getVariablesWithPermission(PermissionTypes::Final);
  EXPECT_EQ(final_write.size(), 2);
  EXPECT_EQ(final_write["species:he:density"], Regions::Boundaries);
  EXPECT_EQ(final_write["species:he:collision_frequency"], Regions::All);

  EXPECT_THROW(example.getVariablesWithPermission(PermissionTypes::None, true),
               BoutException);
  EXPECT_THROW(example.getVariablesWithPermission(PermissionTypes::None, false),
               BoutException);
}

TEST(PermissionsTests, TestSubstitute) {
  Permissions example(
      {{"species:{s1}:collision_frequencies:{s1}_{s2}_coll",
        {Regions::Nowhere, Regions::All, Regions::Nowhere, Regions::Nowhere}},
       readIfSet("d")});
  Permissions example2 = example;

  example.setAccess(
      "{var}", {Regions::Nowhere, Regions::Nowhere, Regions::Interior, Regions::Nowhere});

  example.substitute(
      {{"s1", {"e", "d+"}}, {"s2", {"e", "d+"}}, {"var", {"a", "b", "c", "d"}}});

  auto readable = example.getVariablesWithPermission(PermissionTypes::Read);
  EXPECT_EQ(readable.size(), 4);
  EXPECT_EQ(readable["species:e:collision_frequencies:e_e_coll"], Regions::All);
  EXPECT_EQ(readable["species:e:collision_frequencies:e_d+_coll"], Regions::All);
  EXPECT_EQ(readable["species:d+:collision_frequencies:d+_e_coll"], Regions::All);
  EXPECT_EQ(readable["species:d+:collision_frequencies:d+_d+_coll"], Regions::All);
  auto writable = example.getVariablesWithPermission(PermissionTypes::Write);
  EXPECT_EQ(writable.size(), 3);
  EXPECT_EQ(writable["a"], Regions::Interior);
  EXPECT_EQ(writable["b"], Regions::Interior);
  EXPECT_EQ(writable["c"], Regions::Interior);
  EXPECT_EQ(example.getHighestPermission("d"),
            make_permission(PermissionTypes::ReadIfSet, "d"));
  EXPECT_EQ(example.getVariablesWithPermission(PermissionTypes::ReadIfSet, false).size(),
            8);
  // Trying to perform substitutions that are not present should work fine.
  example.substitute("something", {"one", "two", "three", "four"});
  EXPECT_EQ(example.getVariablesWithPermission(PermissionTypes::ReadIfSet, false).size(),
            8);

  // Trying to perform incomplete substitutions will throw an exception;
  EXPECT_THROW(example2.substitute("s1", {"e", "d+"}), BoutException);
}

TEST(PermissionsTests, TestIO) {
  const Permissions empty({});
  const Permissions single({readOnly("test")});
  const Permissions multiple(
      {readIfSet("a", Regions::Interior), writeBoundary("b"), readWrite("c:d")});
  Permissions new_perm;

  std::stringstream ss1;
  std::stringstream ss2;
  std::stringstream ss3;

  ss1 << empty;
  ss1 >> new_perm;
  EXPECT_EQ(new_perm.getVariablesWithPermission(PermissionTypes::ReadIfSet, false).size(),
            0);

  ss2 << single;
  ss2 >> new_perm;
  EXPECT_EQ(new_perm.getVariablesWithPermission(PermissionTypes::ReadIfSet, false).size(),
            1);
  std::map<std::string, Regions> read_only =
      new_perm.getVariablesWithPermission(PermissionTypes::Read);
  EXPECT_EQ(read_only.size(), 1);
  EXPECT_EQ(read_only["test"], Regions::All);

  ss3 << multiple;
  ss3 >> new_perm;
  EXPECT_EQ(new_perm.getVariablesWithPermission(PermissionTypes::ReadIfSet, false).size(),
            3);
  std::map<std::string, Regions> read_if_set =
      new_perm.getVariablesWithPermission(PermissionTypes::ReadIfSet);
  EXPECT_EQ(read_if_set.size(), 1);
  EXPECT_EQ(read_if_set["a"], Regions::Interior);
  std::map<std::string, Regions> read =
      new_perm.getVariablesWithPermission(PermissionTypes::Read);
  EXPECT_EQ(read.size(), 1);
  EXPECT_EQ(read["b"], Regions::Interior);
  std::map<std::string, Regions> read_write =
      new_perm.getVariablesWithPermission(PermissionTypes::Write);
  EXPECT_EQ(read_write.size(), 1);
  EXPECT_EQ(read_write["c:d"], Regions::All);
  std::map<std::string, Regions> write_final =
      new_perm.getVariablesWithPermission(PermissionTypes::Final);
  EXPECT_EQ(write_final.size(), 1);
  EXPECT_EQ(write_final["b"], Regions::Boundaries);
}
