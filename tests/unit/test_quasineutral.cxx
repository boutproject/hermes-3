
#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/quasineutral.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using QuasineutralTest = FakeMeshFixture;

TEST_F(QuasineutralTest, CreateComponent) {
  Options options = {{"test", {{"charge", 1.0}, {"AA", 1.0}}}};

  Quasineutral component("test", options, nullptr);
}

TEST_F(QuasineutralTest, ChargeNotZero) {
  Options options = {{"test", {{"charge", 0.0}, {"AA", 1.0}}}};

  ASSERT_THROW(Quasineutral component("test", options, nullptr), BoutException);
}

TEST_F(QuasineutralTest, TransformWithoutSpecies) {
  Options options = {{"test", {{"charge", -1.0}, {"AA", 2.0}}}};

  Quasineutral component("test", options, nullptr);

  Options state;
  component.transform(state);
  ASSERT_TRUE(state["species"]["test"].isSet("density"));
  EXPECT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["test"]["density"]), 0.0));

  ASSERT_TRUE(state["species"]["test"].isSet("charge"));
  EXPECT_DOUBLE_EQ(get<BoutReal>(state["species"]["test"]["charge"]), -1.0);

  ASSERT_TRUE(state["species"]["test"].isSet("AA"));
  EXPECT_DOUBLE_EQ(get<BoutReal>(state["species"]["test"]["AA"]), 2.0);
}

TEST_F(QuasineutralTest, TransformOneSpecies) {
  Options options = {{"test", {{"charge", -1.0}, {"AA", 2.0}}}};

  Quasineutral component("test", options, nullptr);
  component.declareAllSpecies({{"i"}, {}});

  Options state{{"species", {{"i", {{"charge", 1.5}, {"density", 2.0}}}}}};
  component.transform(state);
  ASSERT_TRUE(state["species"]["test"].isSet("density"));
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["test"]["density"]), 1.5 * 2.0));

  ASSERT_TRUE(state["species"]["test"].isSet("charge"));
  EXPECT_DOUBLE_EQ(get<BoutReal>(state["species"]["test"]["charge"]), -1.0);

  ASSERT_TRUE(state["species"]["test"].isSet("AA"));
  EXPECT_DOUBLE_EQ(get<BoutReal>(state["species"]["test"]["AA"]), 2.0);
}

TEST_F(QuasineutralTest, TransformTwoSpecies) {
  Options options = {{"test", {{"charge", 1.5}, {"AA", 2.0}}}};

  Quasineutral component("test", options, nullptr);
  component.declareAllSpecies({"e", "i"});

  Options state{{"species",
                 {{"e", {{"charge", -1.0}, {"density", 2.5}}},
                  {"i", {{"charge", 1.0}, {"density", 2.0}}}}}};
  component.transform(state);
  ASSERT_TRUE(state["species"]["test"].isSet("density"));
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["test"]["density"]),
                           -(-1.0 * 2.5 + 1.0 * 2.0) / 1.5));

  ASSERT_TRUE(state["species"]["test"].isSet("charge"));
  EXPECT_DOUBLE_EQ(get<BoutReal>(state["species"]["test"]["charge"]), 1.5);

  ASSERT_TRUE(state["species"]["test"].isSet("AA"));
  EXPECT_DOUBLE_EQ(get<BoutReal>(state["species"]["test"]["AA"]), 2.0);
}
