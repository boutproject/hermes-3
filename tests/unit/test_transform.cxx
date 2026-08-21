#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/transform.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using TransformTest = FakeMeshFixture;

TEST_F(TransformTest, SingleTransform) {
  Options options{{"test", {{"transforms", "a = b"}}}};
  Transform component("test", options, nullptr);

  Options state{{"b", 2.}};
  component.transform(state);
  ASSERT_EQ(get<float>(state["a"]), 2.);
  ASSERT_EQ(get<float>(state["b"]), 2.);
}

TEST_F(TransformTest, SingleTransformOverwrite) {
  Options options{{"test", {{"transforms", "a = b"}}}};
  Transform component("test", options, nullptr);

  Options state{{"a", 1.}, {"b", 2.}};
  ASSERT_NE(state["a"], state["b"]);
  component.transform(state);
  ASSERT_EQ(get<float>(state["a"]), 2.);
  ASSERT_EQ(get<float>(state["b"]), 2.);
}

TEST_F(TransformTest, MultipleTransforms) {
  Options options{{"test", {{"transforms", "a =   b , c=d"}}}};
  Transform component("test", options, nullptr);

  Options state{{"b", 2.}, {"d", 1}};
  component.transform(state);
  ASSERT_EQ(get<float>(state["a"]), 2.);
  ASSERT_EQ(get<int>(state["c"]), 1);
}

TEST_F(TransformTest, BadTransformConfigs) {
  Options options{{"test", {{"transforms", "foo bar"}}}};
  EXPECT_THROW(Transform("test", options, nullptr), BoutException);
}
