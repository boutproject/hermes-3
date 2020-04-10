#include "gtest/gtest.h"

#include "../include/component_scheduler.hxx"

namespace {
struct TestComponent : public Component {
  TestComponent(const std::string&, Options&, const MeshMap&) {}
  void transform(Options &state) { state["answer"] = 42; }
};
  

struct TestMultiply : public Component {
  TestMultiply(const std::string&, Options&, const MeshMap&) {}
  void transform(Options &state) { state["answer"].force(state["answer"].as<int>() * 2); }
};

RegisterComponent<TestComponent> registertestcomponent("testcomponent");
RegisterComponent<TestMultiply> registertestcomponent2("multiply");
} // namespace

TEST(SchedulerTest, OneComponent) {
  MeshMap meshes;

  Options options;
  options["components"] = "testcomponent";
  auto scheduler = ComponentScheduler::create(options, meshes);

  EXPECT_FALSE(options.isSet("answer"));
  scheduler->transform(options);
  ASSERT_TRUE(options.isSet("answer"));
  ASSERT_TRUE(options["answer"] == 42);
}

TEST(SchedulerTest, TwoComponents) {
  MeshMap meshes;

  Options options;
  options["components"] = "testcomponent, multiply";
  auto scheduler = ComponentScheduler::create(options, meshes);

  EXPECT_FALSE(options.isSet("answer"));
  scheduler->transform(options);
  ASSERT_TRUE(options.isSet("answer"));
  ASSERT_TRUE(options["answer"] == 42 * 2);
}
