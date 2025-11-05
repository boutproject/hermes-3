#include "../include/guarded_options.hxx"
#include "gtest/gtest.h"

class GuardedOptionsTests : public testing::Test {
protected:
  GuardedOptionsTests()
      : permissions(
          {{"species:he:density",
            {Permissions::AllRegions, Permissions::Nowhere, Permissions::Nowhere}},
           {"species:he:pressure",
            {Permissions::Nowhere, Permissions::Interior, Permissions::Nowhere}},
           {"species:he:collision_frequency",
            {Permissions::Nowhere, Permissions::Nowhere, Permissions::AllRegions}},
           {"species:he:velocity",
            {Permissions::Boundaries, Permissions::Nowhere, Permissions::Nowhere}},
           {"species:d",
            {Permissions::AllRegions, Permissions::Nowhere, Permissions::Nowhere}},
           {"species:d:pressure",
            {Permissions::Nowhere, Permissions::Interior, Permissions::Nowhere}},
           {"species:d:collision_frequencies",
            {Permissions::Nowhere, Permissions::Boundaries, Permissions::Nowhere}}}),
        opts({{"species",
               {{"he",
                 {{"temperature", 0}, {"density", 1}, {"pressure", 2}, {"velocity", 4}}},
                {"d",
                 {{"pressure", 5},
                  {"velocity", 6},
                  {"collision_frequencies", {{"d_d_coll", 7}, {"d_he_coll", 8}}}}}}}}),
        guarded_opts(&opts, &permissions){};

  Permissions permissions;
  Options opts;
  GuardedOptions guarded_opts;
};

TEST_F(GuardedOptionsTests, TestGet) {
  EXPECT_EQ(guarded_opts["species:he:density"].get(), 1);
  EXPECT_EQ(guarded_opts["species:he:density"].get(Permissions::Boundaries), 1);
  EXPECT_EQ(guarded_opts["species:he:pressure"].get(Permissions::Interior), 2);
  EXPECT_FALSE(guarded_opts["species:he:collision_frequency"].get().isSet());
  EXPECT_EQ(guarded_opts["species"]["he"]["velocity"].get(Permissions::Boundaries), 4);
  EXPECT_EQ(guarded_opts["species"]["d"]["pressure"].get(Permissions::Interior), 5);
  EXPECT_EQ(guarded_opts["species"]["d"]["velocity"].get(Permissions::AllRegions), 6);
  EXPECT_EQ(guarded_opts["species:d:collision_frequencies"]["d_d_coll"].get(
                Permissions::Boundaries),
            7);
  EXPECT_EQ(guarded_opts["species:d:collision_frequencies"].get(
                Permissions::Boundaries)["d_he_coll"],
            8);
  EXPECT_FALSE(guarded_opts["species:d:collision_frequencies:d_t+_cx"]
                   .get(Permissions::Boundaries)
                   .isSet());
}

TEST_F(GuardedOptionsTests, TestGetException) {
  EXPECT_THROW(guarded_opts["species"]["he"]["temperature"].get(), BoutException);
  EXPECT_THROW(guarded_opts["species:he:pressure"].get(), BoutException);
  EXPECT_THROW(guarded_opts["species"]["he"]["velocity"].get(Permissions::Interior),
               BoutException);
  EXPECT_THROW(guarded_opts["species:d:collision_frequencies"].get(Permissions::Interior),
               BoutException);
  EXPECT_THROW(guarded_opts["species:d:pressure"].get(Permissions::Boundaries),
               BoutException);
  EXPECT_THROW(guarded_opts["no_permission"].get(), BoutException);
  EXPECT_THROW(guarded_opts["species:d+:velocity"].get(), BoutException);
}

TEST_F(GuardedOptionsTests, TestGetWritable) {
  auto& he_pressure =
      guarded_opts["species:he:pressure"].getWritable(Permissions::Interior);
  EXPECT_EQ(he_pressure, 2);
  EXPECT_EQ(opts["species"]["he"]["pressure"], 2);
  he_pressure.force(10);
  EXPECT_EQ(he_pressure, 10);
  EXPECT_EQ(opts["species"]["he"]["pressure"], 10);

  auto& he_freq = guarded_opts["species"]["he"]["collision_frequency"].getWritable();
  EXPECT_FALSE(he_freq.isSet());
  he_freq = 11;
  EXPECT_EQ(opts["species"]["he"]["collision_frequency"], 11);

  auto& d_pressure =
      guarded_opts["species"]["d"]["pressure"].getWritable(Permissions::Interior);
  EXPECT_EQ(opts["species"]["d"]["pressure"], 5);
  d_pressure.force(12);
  EXPECT_EQ(opts["species"]["d"]["pressure"], 12);

  auto& d_d_coll = guarded_opts["species:d:collision_frequencies:d_d_coll"].getWritable(
      Permissions::Boundaries);
  EXPECT_EQ(d_d_coll, 7);
  d_d_coll.force(13);
  EXPECT_EQ(opts["species:d:collision_frequencies:d_d_coll"], 13);

  auto& d_tp_coll = guarded_opts["species:d:collision_frequencies:d_t+_coll"].getWritable(
      Permissions::Boundaries);
  EXPECT_FALSE(d_tp_coll.isSet());
  d_tp_coll = 14;
  EXPECT_EQ(opts["species:d:collision_frequencies:d_t+_coll"], 14);
}

TEST_F(GuardedOptionsTests, TestGetWritableException) {
  EXPECT_THROW(guarded_opts["species"]["he"]["temperature"].getWritable(), BoutException);
  EXPECT_THROW(guarded_opts["unset"].getWritable(), BoutException);
  EXPECT_THROW(guarded_opts["species:he:density"].getWritable(), BoutException);
  EXPECT_THROW(guarded_opts["species:he:density"].getWritable(Permissions::Interior),
               BoutException);
  EXPECT_THROW(guarded_opts["species:he:density"].getWritable(Permissions::Boundaries),
               BoutException);
  EXPECT_THROW(guarded_opts["species:he:pressure"].getWritable(), BoutException);
  EXPECT_THROW(guarded_opts["species:he:pressure"].getWritable(Permissions::Boundaries),
               BoutException);
  EXPECT_THROW(guarded_opts["species"]["d"]["velocity"].getWritable(), BoutException);
  EXPECT_THROW(
      guarded_opts["species"]["d"]["pressure"].getWritable(Permissions::Boundaries),
      BoutException);
  EXPECT_THROW(guarded_opts["species:d:collision_frequencies:unset"].getWritable(),
               BoutException);
  EXPECT_THROW(guarded_opts["species:d:collision_frequencies:unset"].getWritable(
                   Permissions::Interior),
               BoutException);
  EXPECT_THROW(
      guarded_opts["species"]["d"]["pressure_suffix"].getWritable(Permissions::Interior),
      BoutException);
}

TEST_F(GuardedOptionsTests, TestUnreadItems) {
  std::map<std::string, Permissions::Regions>
      expected1 = {{"species:he:density", Permissions::AllRegions},
                   {"species:he:velocity", Permissions::Boundaries},
                   {"species:d", Permissions::AllRegions}},
      expected2 = {{"species:he:density", Permissions::Boundaries},
                   {"species:he:velocity", Permissions::Boundaries},
                   {"species:d", Permissions::AllRegions}},
      expected3 = {{"species:d", Permissions::AllRegions}}, expected4;

  EXPECT_EQ(guarded_opts.unreadItems(), expected1);

  guarded_opts["species:he:density"].get(Permissions::Interior);
  guarded_opts["species:d:pressure"].get(Permissions::Interior);
  guarded_opts["species:d:collision_frequencies:d_d_coll"].getWritable(
      Permissions::Boundaries);
  EXPECT_EQ(guarded_opts.unreadItems(), expected2);

  guarded_opts["species"]["he"]["density"].get();
  guarded_opts["species:he:velocity"].get(Permissions::Boundaries);
  EXPECT_EQ(guarded_opts.unreadItems(), expected3);

  guarded_opts["species:d:velocity"].get();
  EXPECT_EQ(guarded_opts.unreadItems(), expected4);
}

TEST_F(GuardedOptionsTests, TestUnwrittenItems) {
  std::map<std::string, Permissions::Regions>
      expected1 = {{"species:he:pressure", Permissions::Interior},
                   {"species:he:collision_frequency", Permissions::AllRegions},
                   {"species:d:pressure", Permissions::Interior},
                   {"species:d:collision_frequencies", Permissions::Boundaries}},
      expected2 = {{"species:he:collision_frequency", Permissions::AllRegions},
                   {"species:d:pressure", Permissions::Interior}},
      expected3;

  EXPECT_EQ(guarded_opts.unwrittenItems(), expected1);

  guarded_opts["species:he:pressure"].get(Permissions::Interior);
  guarded_opts["species:he:collision_frequency"].get();
  EXPECT_EQ(guarded_opts.unwrittenItems(), expected1);

  guarded_opts["species"]["he"]["pressure"].getWritable(Permissions::Interior);
  guarded_opts["species:d:collision_frequencies:d_d_coll"].getWritable(
      Permissions::Boundaries);
  EXPECT_EQ(guarded_opts.unwrittenItems(), expected2);

  guarded_opts["species:he:collision_frequency"].getWritable();
  guarded_opts["species:d:pressure"].getWritable(Permissions::Interior);
  EXPECT_EQ(guarded_opts.unwrittenItems(), expected3);
}

TEST_F(GuardedOptionsTests, TestNullOptions) {
  GuardedOptions null_opts(nullptr, &permissions);
  EXPECT_THROW(null_opts["species:he:collision_frequency"], BoutException);
  EXPECT_THROW(null_opts.get(), BoutException);
  EXPECT_THROW(null_opts.getWritable(), BoutException);
}

TEST_F(GuardedOptionsTests, TestNullPermissions) {
  GuardedOptions null_perms(&opts, nullptr);
  GuardedOptions sub_opt = null_perms["species:he:pressure"];
  EXPECT_THROW(sub_opt.get(Permissions::Interior), BoutException);
  EXPECT_THROW(sub_opt.getWritable(Permissions::Interior), BoutException);
}

TEST_F(GuardedOptionsTests, TestGetChildren) {
  std::map<std::string, GuardedOptions> guarded_children = guarded_opts["species"].getChildren();
  EXPECT_EQ(guarded_children.size(), 2);
  EXPECT_EQ(guarded_children.count("he"), 1);
  EXPECT_EQ(guarded_children.count("d"), 1);
  EXPECT_EQ(&(guarded_children["d"].get()), &(opts["species"]["d"]));
  // We do not have access to the whole "he" section
  EXPECT_THROW(guarded_children["he"].get(), BoutException);
}

TEST_F(GuardedOptionsTests, TestIsThisSection) {
  EXPECT_TRUE(guarded_opts.isSection());
  EXPECT_TRUE(guarded_opts["species:d:collision_frequencies"].isSection());
  EXPECT_FALSE(guarded_opts["species:he:temperature"].isSection());
  EXPECT_TRUE(guarded_opts["species"]["he"]["collision_frequency"].isSection());
}

TEST_F(GuardedOptionsTests, TestIsChildSection) {
  EXPECT_TRUE(guarded_opts.isSection(""));
  // Unforunately Operator::isSection does not support full paths
  EXPECT_FALSE(guarded_opts.isSection("species:he"));
  EXPECT_FALSE(guarded_opts["species"]["he"].isSection("temperature"));
  EXPECT_FALSE(guarded_opts["species:he"].isSection("pressure"));
  EXPECT_TRUE(guarded_opts["species:d"].isSection("collision_frequencies"));
}

TEST_F(GuardedOptionsTests, TestIsThisSet) {
  EXPECT_TRUE(guarded_opts["species"]["he"]["temperature"].isSet());
  EXPECT_TRUE(guarded_opts["species:he:pressure"].isSet());
  EXPECT_TRUE(guarded_opts["species"]["he"]["velocity"].isSet());
  EXPECT_FALSE(guarded_opts["species"]["he"]["collision_frequency"].isSet());
  EXPECT_FALSE(guarded_opts["unset"].isSet());
}

TEST_F(GuardedOptionsTests, TestIsChildSet) {
  EXPECT_FALSE(guarded_opts["species:he:pressure"].isSet("test"));
  EXPECT_FALSE(guarded_opts["species:he"].isSet("collision_frequency"));
  EXPECT_TRUE(guarded_opts["species"]["he"].isSet("temperature"));
  // Unforunately Operator::isSet does not support full paths
  EXPECT_FALSE(guarded_opts.isSet("species:he:temperature"));
  EXPECT_TRUE(guarded_opts["species"]["d"]["collision_frequencies"].isSet("d_d_coll"));
}
