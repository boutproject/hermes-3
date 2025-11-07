#include "gtest/gtest.h"

#include "../include/permissions.hxx"

auto make_access = std::make_pair<bool, std::string>;
auto make_permission = std::make_pair<Permissions::PermissionTypes, std::string>;

TEST(PermissionsTests, TestCanAccess) {
  Permissions example({
      readIfSet("species:he:charge"),
      readOnly("species:he:density"),
      // Read and write permissions for pressure in the interior region
      {"species:he:pressure",
       {Permissions::Nowhere, Permissions::Nowhere, Permissions::Interior,
        Permissions::Nowhere}},
      // Set the final value for collision frequency
      writeFinal("species:he:collision_frequency"),
      // Only allow reading of boundary velocity
      {"species:he:velocity",
       {Permissions::Nowhere, Permissions::Boundaries, Permissions::Nowhere,
        Permissions::Nowhere}},
      readOnly("species:d"),
      {"species:d:pressure",
       {Permissions::Nowhere, Permissions::Nowhere, Permissions::Interior,
        Permissions::Nowhere}},
      {"species:d:collision_frequencies",
       {Permissions::Nowhere, Permissions::Nowhere, Permissions::Boundaries,
        Permissions::Nowhere}},
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
  EXPECT_EQ(example.canAccess("species:he:charge", Permissions::ReadIfSet),
            make_access(true, "species:he:charge"));
  EXPECT_EQ(example.canAccess("species:he:density", Permissions::ReadIfSet),
            make_access(true, "species:he:density"));
  EXPECT_EQ(example.canAccess("unset", Permissions::ReadIfSet), no_access);

  // Check whether we have write permission for the variables across the entire domain
  EXPECT_EQ(example.canAccess("species:he:charge", Permissions::Write), no_access);
  EXPECT_EQ(example.canAccess("species:he:density", Permissions::Write), no_access);
  EXPECT_EQ(example.canAccess("species:he:pressure", Permissions::Write), no_access);
  EXPECT_EQ(example.canAccess("species:he:collision_frequency", Permissions::Write),
            make_access(true, "species:he:collision_frequency"));
  EXPECT_EQ(example.canAccess("species:he:velocity", Permissions::Write), no_access);
  EXPECT_EQ(example.canAccess("unset", Permissions::Write), no_access);

  // Check whether we have read permission at the boundaries
  EXPECT_EQ(
      example.canAccess("species:he:charge", Permissions::Read, Permissions::Boundaries),
      no_access);
  EXPECT_EQ(
      example.canAccess("species:he:density", Permissions::Read, Permissions::Boundaries),
      make_access(true, "species:he:density"));
  EXPECT_EQ(example.canAccess("species:he:pressure", Permissions::Read,
                              Permissions::Boundaries),
            no_access);
  EXPECT_EQ(example.canAccess("species:he:collision_frequency", Permissions::Read,
                              Permissions::Boundaries),
            make_access(true, "species:he:collision_frequency"));
  EXPECT_EQ(example.canAccess("species:he:velocity", Permissions::Read,
                              Permissions::Boundaries),
            make_access(true, "species:he:velocity"));
  EXPECT_EQ(example.canAccess("unset", Permissions::Read, Permissions::Boundaries),
            no_access);

  // Check permissions set for whole sections
  EXPECT_EQ(example.canAccess("species:d:pressure"), no_access);
  EXPECT_EQ(example.canAccess("species:d:pressure", Permissions::Write), no_access);
  EXPECT_EQ(
      example.canAccess("species:d:pressure", Permissions::Write, Permissions::Interior),
      make_access(true, "species:d:pressure"));
  EXPECT_EQ(example.canAccess("species:d:velocity"), make_access(true, "species:d"));
  EXPECT_EQ(example.canAccess("species:d:velocity", Permissions::Write), no_access);
  EXPECT_EQ(
      example.canAccess("species:d:velocity", Permissions::Read, Permissions::Boundaries),
      make_access(true, "species:d"));
  EXPECT_EQ(example.canAccess("species:d:collision_frequencies:d_d_coll"), no_access);
  EXPECT_EQ(
      example.canAccess("species:d:collision_frequencies:d_d_coll", Permissions::Write),
      no_access);
  EXPECT_EQ(example.canAccess("species:d:collision_frequencies:d_d_coll",
                              Permissions::Write, Permissions::Boundaries),
            make_access(true, "species:d:collision_frequencies"));

  // Check permissions for a species that might be mistaken for one of
  // the sections we've given permissions for
  EXPECT_EQ(example.canAccess("species:d+"), no_access);
  EXPECT_EQ(example.canAccess("species:d+", Permissions::Write), no_access);
  EXPECT_EQ(example.canAccess("species:d+", Permissions::Read, Permissions::Interior),
            no_access);
  EXPECT_EQ(example.canAccess("species:d+", Permissions::Read, Permissions::Boundaries),
            no_access);
}

TEST(PermissionsTests, TestGetHighestPermission) {
  Permissions example({
      {"species:he:charge",
       {Permissions::AllRegions, Permissions::Nowhere, Permissions::Nowhere,
        Permissions::Nowhere}},
      {"species:he:density",
       {Permissions::Nowhere, Permissions::AllRegions, Permissions::Nowhere,
        Permissions::Boundaries}},
      // Read and write permissions for pressure in the interior region
      {"species:he:pressure",
       {Permissions::Nowhere, Permissions::Boundaries, Permissions::Interior,
        Permissions::Nowhere}},
      // Set the final value for collision frequency
      {"species:he:collision_frequency",
       {Permissions::Nowhere, Permissions::Interior, Permissions::Nowhere,
        Permissions::AllRegions}},
      // Only allow reading of boundary velocity
      {"species:he:velocity",
       {Permissions::Nowhere, Permissions::Boundaries, Permissions::Nowhere,
        Permissions::Nowhere}},
      {"species:d",
       {Permissions::Nowhere, Permissions::AllRegions, Permissions::Nowhere,
        Permissions::Nowhere}},
      {"species:d:pressure",
       {Permissions::Nowhere, Permissions::Nowhere, Permissions::Interior,
        Permissions::Nowhere}},
      {"species:d:collision_frequencies",
       {Permissions::Nowhere, Permissions::Nowhere, Permissions::Boundaries,
        Permissions::Nowhere}},
  });

  auto no_permission = make_permission(Permissions::None, "");

  // Get the highest permission that covers the entire domain
  EXPECT_EQ(example.getHighestPermission("species:he:charge"),
            make_permission(Permissions::ReadIfSet, "species:he:charge"));
  EXPECT_EQ(example.getHighestPermission("species:he:density"),
            make_permission(Permissions::Read, "species:he:density"));
  EXPECT_EQ(example.getHighestPermission("species:he:pressure"),
            make_permission(Permissions::Read, "species:he:pressure"));
  EXPECT_EQ(example.getHighestPermission("species:he:collision_frequency"),
            make_permission(Permissions::Final, "species:he:collision_frequency"));
  EXPECT_EQ(example.getHighestPermission("species:he:velocity"),
            make_permission(Permissions::None, "species:he:velocity"));
  EXPECT_EQ(example.getHighestPermission("species:d:pressure"),
            make_permission(Permissions::None, "species:d:pressure"));
  EXPECT_EQ(example.getHighestPermission("species:d:velocity"),
            make_permission(Permissions::Read, "species:d"));
  EXPECT_EQ(example.getHighestPermission("species:d:collision_frequencies:d_d_coll"),
            make_permission(Permissions::None, "species:d:collision_frequencies"));
  EXPECT_EQ(example.getHighestPermission("unset"), no_permission);

  // Get the highest permission on the boundaries
  EXPECT_EQ(example.getHighestPermission("species:he:charge", Permissions::Boundaries),
            make_permission(Permissions::ReadIfSet, "species:he:charge"));
  EXPECT_EQ(example.getHighestPermission("species:he:density", Permissions::Boundaries),
            make_permission(Permissions::Final, "species:he:density"));
  EXPECT_EQ(example.getHighestPermission("species:he:pressure", Permissions::Boundaries),
            make_permission(Permissions::Read, "species:he:pressure"));
  EXPECT_EQ(example.getHighestPermission("species:he:collision_frequency",
                                         Permissions::Boundaries),
            make_permission(Permissions::Final, "species:he:collision_frequency"));
  EXPECT_EQ(example.getHighestPermission("species:he:velocity", Permissions::Boundaries),
            make_permission(Permissions::Read, "species:he:velocity"));
  EXPECT_EQ(example.getHighestPermission("species:d:pressure", Permissions::Boundaries),
            make_permission(Permissions::None, "species:d:pressure"));
  EXPECT_EQ(example.getHighestPermission("species:d:velocity", Permissions::Boundaries),
            make_permission(Permissions::Read, "species:d"));
  EXPECT_EQ(example.getHighestPermission("species:d:collision_frequencies:d_d_coll",
                                         Permissions::Boundaries),
            make_permission(Permissions::Write, "species:d:collision_frequencies"));
  EXPECT_EQ(example.getHighestPermission("unset", Permissions::Boundaries),
            no_permission);

  // Get the highest permission on the interior
  EXPECT_EQ(example.getHighestPermission("species:he:charge", Permissions::Interior),
            make_permission(Permissions::ReadIfSet, "species:he:charge"));
  EXPECT_EQ(example.getHighestPermission("species:he:density", Permissions::Interior),
            make_permission(Permissions::Read, "species:he:density"));
  EXPECT_EQ(example.getHighestPermission("species:he:pressure", Permissions::Interior),
            make_permission(Permissions::Write, "species:he:pressure"));
  EXPECT_EQ(example.getHighestPermission("species:he:collision_frequency",
                                         Permissions::Interior),
            make_permission(Permissions::Final, "species:he:collision_frequency"));
  EXPECT_EQ(example.getHighestPermission("species:he:velocity", Permissions::Interior),
            make_permission(Permissions::None, "species:he:velocity"));
  EXPECT_EQ(example.getHighestPermission("species:d:pressure", Permissions::Interior),
            make_permission(Permissions::Write, "species:d:pressure"));
  EXPECT_EQ(example.getHighestPermission("species:d:velocity", Permissions::Interior),
            make_permission(Permissions::Read, "species:d"));
  EXPECT_EQ(example.getHighestPermission("species:d:collision_frequencies:d_d_coll",
                                         Permissions::Interior),
            make_permission(Permissions::None, "species:d:collision_frequencies"));
  EXPECT_EQ(example.getHighestPermission("unset", Permissions::Interior), no_permission);

  // Check the permission for the "Nowhere" region is always "None"
  EXPECT_EQ(example.getHighestPermission("species:he:charge", Permissions::Nowhere),
            no_permission);
  EXPECT_EQ(example.getHighestPermission("species:he:density", Permissions::Nowhere),
            no_permission);
  EXPECT_EQ(example.getHighestPermission("species:he:pressure", Permissions::Nowhere),
            no_permission);
  EXPECT_EQ(example.getHighestPermission("species:he:collision_frequency",
                                         Permissions::Nowhere),
            no_permission);
  EXPECT_EQ(example.getHighestPermission("species:he:velocity", Permissions::Nowhere),
            no_permission);
  EXPECT_EQ(example.getHighestPermission("species:d:pressure", Permissions::Nowhere),
            no_permission);
  EXPECT_EQ(example.getHighestPermission("species:d:velocity", Permissions::Nowhere),
            no_permission);
  EXPECT_EQ(example.getHighestPermission("species:d:collision_frequencies:d_d_coll",
                                         Permissions::Nowhere),
            no_permission);
  EXPECT_EQ(example.getHighestPermission("unset", Permissions::Nowhere), no_permission);

  // Check permissions for a species that might be mistaken for one of
  // the sections we've given permissions for
  EXPECT_EQ(example.getHighestPermission("species:d+"), no_permission);
  EXPECT_EQ(example.getHighestPermission("species:d+", Permissions::Interior),
            no_permission);
  EXPECT_EQ(example.getHighestPermission("species:d+", Permissions::Boundaries),
            no_permission);
}

TEST(PermissionsTests, TestSetAccess) {
  Permissions example({
      {"species:he:density",
       {Permissions::Nowhere, Permissions::AllRegions, Permissions::Nowhere,
        Permissions::Nowhere}},
      // Read and write permissions for pressure in the interior region
      {"species:he:pressure",
       {Permissions::Nowhere, Permissions::Nowhere, Permissions::Interior,
        Permissions::Nowhere}},
  });

  EXPECT_EQ(example.getHighestPermission("species:he:density"),
            make_permission(Permissions::Read, "species:he:density"));
  example.setAccess("species:he:density",
                    {Permissions::Nowhere, Permissions::Nowhere, Permissions::Boundaries,
                     Permissions::Nowhere});
  EXPECT_EQ(example.getHighestPermission("species:he:density"),
            make_permission(Permissions::None, "species:he:density"));
  EXPECT_EQ(example.getHighestPermission("species:he:density", Permissions::Boundaries),
            make_permission(Permissions::Write, "species:he:density"));

  EXPECT_EQ(example.getHighestPermission("species:he:pressure", Permissions::Interior),
            make_permission(Permissions::Write, "species:he:pressure"));
  EXPECT_EQ(example.getHighestPermission("species:he:pressure", Permissions::AllRegions),
            make_permission(Permissions::None, "species:he:pressure"));
  example.setAccess("species:he:pressure", {Permissions::Nowhere, Permissions::AllRegions,
                                            Permissions::Nowhere, Permissions::Nowhere});
  EXPECT_EQ(example.getHighestPermission("species:he:pressure", Permissions::Interior),
            make_permission(Permissions::Read, "species:he:pressure"));
  EXPECT_EQ(example.getHighestPermission("species:he:pressure"),
            make_permission(Permissions::Read, "species:he:pressure"));

  EXPECT_EQ(example.getHighestPermission("unset", Permissions::Interior),
            make_permission(Permissions::None, ""));
  EXPECT_EQ(example.getHighestPermission("unset", Permissions::Boundaries),
            make_permission(Permissions::None, ""));
  example.setAccess("unset", {Permissions::Nowhere, Permissions::Interior,
                              Permissions::Nowhere, Permissions::Boundaries});
  EXPECT_EQ(example.getHighestPermission("unset", Permissions::AllRegions),
            make_permission(Permissions::Read, "unset"));
  EXPECT_EQ(example.getHighestPermission("unset", Permissions::Boundaries),
            make_permission(Permissions::Final, "unset"));
}

TEST(PermissionsTests, TestGetVariablesWithPermissions) {
  Permissions example({{"species:he:density",
                        {Permissions::Nowhere, Permissions::AllRegions,
                         Permissions::Nowhere, Permissions::Boundaries}},
                       // Read and write permissions for pressure in the interior region
                       {"species:he:pressure",
                        {Permissions::Nowhere, Permissions::Boundaries,
                         Permissions::Interior, Permissions::Nowhere}},
                       // Set the final value for collision frequency
                       {"species:he:collision_frequency",
                        {Permissions::Nowhere, Permissions::Interior,
                         Permissions::Nowhere, Permissions::AllRegions}},
                       // Only allow reading of boundary velocity
                       {"species:he:velocity",
                        {Permissions::Nowhere, Permissions::Boundaries,
                         Permissions::Nowhere, Permissions::Nowhere}}});

  auto read_only = example.getVariablesWithPermission(Permissions::Read);
  EXPECT_EQ(read_only.size(), 3);
  EXPECT_EQ(read_only["species:he:density"], Permissions::Interior);
  EXPECT_EQ(read_only["species:he:pressure"], Permissions::Boundaries);
  EXPECT_EQ(read_only["species:he:velocity"], Permissions::Boundaries);

  auto readable = example.getVariablesWithPermission(Permissions::Read, false);
  EXPECT_EQ(readable.size(), 4);
  EXPECT_EQ(readable["species:he:density"], Permissions::AllRegions);
  EXPECT_EQ(readable["species:he:pressure"], Permissions::AllRegions);
  EXPECT_EQ(readable["species:he:collision_frequency"], Permissions::AllRegions);
  EXPECT_EQ(readable["species:he:velocity"], Permissions::Boundaries);

  auto write_nonfinal = example.getVariablesWithPermission(Permissions::Write, true);
  EXPECT_EQ(write_nonfinal.size(), 1);
  EXPECT_EQ(write_nonfinal["species:he:pressure"], Permissions::Interior);

  auto writable = example.getVariablesWithPermission(Permissions::Write, false);
  EXPECT_EQ(writable.size(), 3);
  EXPECT_EQ(writable["species:he:density"], Permissions::Boundaries);
  EXPECT_EQ(writable["species:he:pressure"], Permissions::Interior);
  EXPECT_EQ(writable["species:he:collision_frequency"], Permissions::AllRegions);

  auto final_write = example.getVariablesWithPermission(Permissions::Final);
  EXPECT_EQ(final_write.size(), 2);
  EXPECT_EQ(final_write["species:he:density"], Permissions::Boundaries);
  EXPECT_EQ(final_write["species:he:collision_frequency"], Permissions::AllRegions);
}

TEST(PermissionsTests, TestSubstitute) {
  Permissions example({{"species:{s1}:collision_frequencies:{s1}_{s2}_coll",
                        {Permissions::Nowhere, Permissions::AllRegions,
                         Permissions::Nowhere, Permissions::Nowhere}},
                       readIfSet("d")});

  example.setAccess("{var}", {Permissions::Nowhere, Permissions::Nowhere,
                              Permissions::Interior, Permissions::Nowhere});

  example.substitute("s1", {"e", "d+"});
  example.substitute("s2", {"e", "d+"});

  auto readable = example.getVariablesWithPermission(Permissions::Read, false);
  EXPECT_EQ(readable.size(), 5);
  EXPECT_EQ(readable["{var}"], Permissions::Interior);
  EXPECT_EQ(readable["species:e:collision_frequencies:e_e_coll"],
            Permissions::AllRegions);
  EXPECT_EQ(readable["species:e:collision_frequencies:e_d+_coll"],
            Permissions::AllRegions);
  EXPECT_EQ(readable["species:d+:collision_frequencies:d+_e_coll"],
            Permissions::AllRegions);
  EXPECT_EQ(readable["species:d+:collision_frequencies:d+_d+_coll"],
            Permissions::AllRegions);

  example.substitute("var", {"a", "b", "c", "d"});

  auto writable = example.getVariablesWithPermission(Permissions::Write);
  EXPECT_EQ(writable.size(), 3);
  EXPECT_EQ(writable["a"], Permissions::Interior);
  EXPECT_EQ(writable["b"], Permissions::Interior);
  EXPECT_EQ(writable["c"], Permissions::Interior);

  EXPECT_EQ(example.getHighestPermission("d"),
            make_permission(Permissions::ReadIfSet, "d"));
}
