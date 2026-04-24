#include <map>
#include <string>

#include <bout/boutexception.hxx>
#include <bout/options.hxx>
#include "gtest/gtest.h"

#include "../include/guarded_options.hxx"
#include "../include/permissions.hxx"

class GuardedOptionsTests : public testing::Test {
protected:
  GuardedOptionsTests()
      : permissions(
          {readIfSet("species:he:AA"),
           readIfSet("species:he:charge"),
           readOnly("species:he:density"),
           {"species:he:pressure",
            {Regions::Nowhere, Regions::Nowhere, Regions::Interior, Regions::Nowhere}},
           writeFinal("species:he:collision_frequency"),
           {"species:he:velocity",
            {Regions::Nowhere, Regions::Boundaries, Regions::Nowhere, Regions::Nowhere}},
           readOnly("species:d"),
           {"species:d:pressure",
            {Regions::Nowhere, Regions::Nowhere, Regions::Interior, Regions::Nowhere}},
           {"species:d:collision_frequencies",
            {Regions::Nowhere, Regions::Nowhere, Regions::Boundaries, Regions::Nowhere}},
           readIfSet("fields:phi"),
           readOnly("unused:option")}),
        opts({{"species",
               {{"he",
                 {{"charge", 0},
                  {"temperature", 0},
                  {"density", 1},
                  {"pressure", 2},
                  {"velocity", 4}}},
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
  EXPECT_EQ(guarded_opts["species:he:charge"].get(), 0);
  EXPECT_EQ(guarded_opts["species:he:density"].get(), 1);
  EXPECT_EQ(guarded_opts["species:he:density"].get(Regions::Boundaries), 1);
  EXPECT_EQ(guarded_opts["species:he:pressure"].get(Regions::Interior), 2);
  EXPECT_FALSE(guarded_opts["species:he:collision_frequency"].get().isSet());
  EXPECT_EQ(guarded_opts["species"]["he"]["velocity"].get(Regions::Boundaries), 4);
  EXPECT_EQ(guarded_opts["species"]["d"]["pressure"].get(Regions::Interior), 5);
  EXPECT_EQ(guarded_opts["species"]["d"]["velocity"].get(Regions::All), 6);
  EXPECT_EQ(guarded_opts["species:d:collision_frequencies"]["d_d_coll"].get(
                Regions::Boundaries),
            7);
  EXPECT_EQ(guarded_opts["species:d:collision_frequencies"].get(
                Regions::Boundaries)["d_he_coll"],
            8);
  EXPECT_FALSE(guarded_opts["species:d:collision_frequencies:d_t+_cx"]
                   .get(Regions::Boundaries)
                   .isSet());
}

#if CHECKLEVEL >= 1
TEST_F(GuardedOptionsTests, TestGetException) {
  EXPECT_THROW(guarded_opts["species:he:AA"].get(), BoutException);
  EXPECT_THROW(guarded_opts["species"]["he"]["temperature"].get(), BoutException);
  EXPECT_THROW(guarded_opts["species:he:pressure"].get(), BoutException);
  EXPECT_THROW(guarded_opts["species"]["he"]["velocity"].get(Regions::Interior),
               BoutException);
  EXPECT_THROW(guarded_opts["species:d:collision_frequencies"].get(Regions::Interior),
               BoutException);
  EXPECT_THROW(guarded_opts["species:d:pressure"].get(Regions::Boundaries),
               BoutException);
  EXPECT_THROW(guarded_opts["no_permission"].get(), BoutException);
  EXPECT_THROW(guarded_opts["species:d+:velocity"].get(), BoutException);
}
#endif

TEST_F(GuardedOptionsTests, TestGetWritable) {
  auto& he_pressure = guarded_opts["species:he:pressure"].getWritable(Regions::Interior);
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
      guarded_opts["species"]["d"]["pressure"].getWritable(Regions::Interior);
  EXPECT_EQ(opts["species"]["d"]["pressure"], 5);
  d_pressure.force(12);
  EXPECT_EQ(opts["species"]["d"]["pressure"], 12);

  auto& d_d_coll = guarded_opts["species:d:collision_frequencies:d_d_coll"].getWritable(
      Regions::Boundaries);
  EXPECT_EQ(d_d_coll, 7);
  d_d_coll.force(13);
  EXPECT_EQ(opts["species:d:collision_frequencies:d_d_coll"], 13);

  auto& d_tp_coll = guarded_opts["species:d:collision_frequencies:d_t+_coll"].getWritable(
      Regions::Boundaries);
  EXPECT_FALSE(d_tp_coll.isSet());
  d_tp_coll = 14;
  EXPECT_EQ(opts["species:d:collision_frequencies:d_t+_coll"], 14);
}

#if CHECKLEVEL >= 1
TEST_F(GuardedOptionsTests, TestGetWritableException) {
  EXPECT_THROW(guarded_opts["species"]["he"]["temperature"].getWritable(), BoutException);
  EXPECT_THROW(guarded_opts["unset"].getWritable(), BoutException);
  EXPECT_THROW(guarded_opts["species:he:density"].getWritable(), BoutException);
  EXPECT_THROW(guarded_opts["species:he:density"].getWritable(Regions::Interior),
               BoutException);
  EXPECT_THROW(guarded_opts["species:he:density"].getWritable(Regions::Boundaries),
               BoutException);
  EXPECT_THROW(guarded_opts["species:he:pressure"].getWritable(), BoutException);
  EXPECT_THROW(guarded_opts["species:he:pressure"].getWritable(Regions::Boundaries),
               BoutException);
  EXPECT_THROW(guarded_opts["species"]["d"]["velocity"].getWritable(), BoutException);
  EXPECT_THROW(guarded_opts["species"]["d"]["pressure"].getWritable(Regions::Boundaries),
               BoutException);
  EXPECT_THROW(guarded_opts["species:d:collision_frequencies:unset"].getWritable(),
               BoutException);
  EXPECT_THROW(guarded_opts["species:d:collision_frequencies:unset"].getWritable(
                   Regions::Interior),
               BoutException);
  EXPECT_THROW(
      guarded_opts["species"]["d"]["pressure_suffix"].getWritable(Regions::Interior),
      BoutException);
}
#endif

TEST_F(GuardedOptionsTests, TestUnreadItems) {
  const std::map<std::string, Regions> expected1 = {
      {"species:he:charge", Regions::All},
      {"species:he:density", Regions::All},
      {"species:he:velocity", Regions::Boundaries},
      {"species:d", Regions::All},
      {"unused:option", Regions::All}};
  const std::map<std::string, Regions> expected2 = {
      {"species:he:charge", Regions::All},
      {"species:he:density", Regions::Boundaries},
      {"species:he:velocity", Regions::Boundaries},
      {"species:d", Regions::All},
      {"unused:option", Regions::All}};
  const std::map<std::string, Regions> expected3 = {{"species:d", Regions::All}};
  const std::map<std::string, Regions> expected4;

#if CHECKLEVEL >= 999
  EXPECT_EQ(guarded_opts.unreadItems(), expected1);
#else
  EXPECT_THROW(guarded_opts.unreadItems(), BoutException);
#endif

  guarded_opts["species:he:density"].get(Regions::Interior);
  guarded_opts["species:d:pressure"].get(Regions::Interior);
  guarded_opts["species:d:collision_frequencies:d_d_coll"].getWritable(
      Regions::Boundaries);
#if CHECKLEVEL >= 999
  EXPECT_EQ(guarded_opts.unreadItems(), expected2);
#else
  EXPECT_THROW(guarded_opts.unreadItems(), BoutException);
#endif

  guarded_opts["species"]["he"]["charge"].get();
  guarded_opts["species"]["he"]["density"].get();
  guarded_opts["species:he:velocity"].get(Regions::Boundaries);
  EXPECT_FALSE(guarded_opts["unused"]["option"].get().isSet());
#if CHECKLEVEL >= 999
  EXPECT_EQ(guarded_opts.unreadItems(), expected3);
#else
  EXPECT_THROW(guarded_opts.unreadItems(), BoutException);
#endif

  guarded_opts["species:d:velocity"].get();
#if CHECKLEVEL >= 999
  EXPECT_EQ(guarded_opts.unreadItems(), expected4);
#else
  EXPECT_THROW(guarded_opts.unreadItems(), BoutException);
#endif
}

TEST_F(GuardedOptionsTests, TestUnwrittenItems) {
  const std::map<std::string, Regions> expected1 = {
      {"species:he:pressure", Regions::Interior},
      {"species:he:collision_frequency", Regions::All},
      {"species:d:pressure", Regions::Interior},
      {"species:d:collision_frequencies", Regions::Boundaries}};
  const std::map<std::string, Regions> expected2 = {
      {"species:he:collision_frequency", Regions::All},
      {"species:d:pressure", Regions::Interior}};
  const std::map<std::string, Regions> expected3;

#if CHECKLEVEL >= 999
  EXPECT_EQ(guarded_opts.unwrittenItems(), expected1);
#else
  EXPECT_THROW(guarded_opts.unwrittenItems(), BoutException);
#endif

  guarded_opts["species:he:pressure"].get(Regions::Interior);
  guarded_opts["species:he:collision_frequency"].get();
#if CHECKLEVEL >= 999
  EXPECT_EQ(guarded_opts.unwrittenItems(), expected1);
#else
  EXPECT_THROW(guarded_opts.unwrittenItems(), BoutException);
#endif

  guarded_opts["species"]["he"]["pressure"].getWritable(Regions::Interior);
  guarded_opts["species:d:collision_frequencies:d_d_coll"].getWritable(
      Regions::Boundaries);
#if CHECKLEVEL >= 999
  EXPECT_EQ(guarded_opts.unwrittenItems(), expected2);
#else
  EXPECT_THROW(guarded_opts.unwrittenItems(), BoutException);
#endif

  guarded_opts["species:he:collision_frequency"].getWritable();
  guarded_opts["species:d:pressure"].getWritable(Regions::Interior);
#if CHECKLEVEL >= 999
  EXPECT_EQ(guarded_opts.unwrittenItems(), expected3);
#else
  EXPECT_THROW(guarded_opts.unwrittenItems(), BoutException);
#endif
}

TEST_F(GuardedOptionsTests, TestNullOptions) {
  EXPECT_THROW(GuardedOptions(nullptr, &permissions), BoutException);
}

TEST_F(GuardedOptionsTests, TestNullPermissions) {
  EXPECT_THROW(GuardedOptions(&opts, nullptr), BoutException);
}

TEST_F(GuardedOptionsTests, TestGetChildren) {
  std::map<std::string, GuardedOptions> guarded_children = guarded_opts["species"].getChildren();
  EXPECT_EQ(guarded_children.size(), 2);
  EXPECT_EQ(guarded_children.count("he"), 1);
  EXPECT_EQ(guarded_children.count("d"), 1);
  EXPECT_EQ(&(guarded_children.at("d").get()), &(opts["species"]["d"]));
#if CHECKLEVEL >= 1
  // We do not have access to the whole "he" section
  EXPECT_THROW(guarded_children.at("he").get(), BoutException);
#endif
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

TEST_F(GuardedOptionsTests, TestUnsetNotInChildren) {
  // This is a test for a bug that was found in the initial implementation
  auto children = guarded_opts.getChildren();
  EXPECT_EQ(children.count("fields"), 0);
  EXPECT_EQ(children.count("unused"), 0);
}

TEST_F(GuardedOptionsTests, TestEqual) {
  EXPECT_EQ(guarded_opts, guarded_opts);
  EXPECT_EQ(guarded_opts["species:he:pressure"],
            guarded_opts["species"]["he"]["pressure"]);
  EXPECT_NE(guarded_opts["species:he:pressure"], guarded_opts["species:he:velocity"]);
  EXPECT_NE(guarded_opts, GuardedOptions(&opts, &permissions));
}
