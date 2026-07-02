#include "gtest/gtest.h"

#include "../../include/component_scheduler.hxx"

namespace {
struct TestComponent : public NamedComponent<TestComponent> {
  TestComponent(std::string name, Options&, Solver*)
      : NamedComponent(name, {readWrite("answer")}) {}

  static constexpr auto type = "testcomponent";

private:
  void transform_impl(GuardedOptions& state) override {
    state["answer"].getWritable() = 42;
  }
};

struct TestMultiply : public NamedComponent<TestMultiply> {
  TestMultiply(std::string name, Options&, Solver*)
      : NamedComponent(name, {writeFinal("answer")}) {}

  static constexpr auto type = "multiply";

private:
  void transform_impl(GuardedOptions& state) override {
    // Note: Using set<>() and get<>() for quicker access, avoiding printing
    //       getNonFinal needs to be used because we set the value afterwards
    set(state["answer"], getNonFinal<int>(state["answer"]) * 2);
  }
};

RegisterComponent<TestComponent> registertestcomponent;
RegisterComponent<TestMultiply> registertestcomponent2;
} // namespace

TEST(SchedulerTest, OneComponent) {
  Options options;
  options["components"] = "testcomponent";
  auto scheduler = ComponentScheduler::create(options, options, nullptr);

  EXPECT_FALSE(options.isSet("answer"));
  scheduler->transform(options);
  ASSERT_TRUE(options.isSet("answer"));
  ASSERT_TRUE(options["answer"] == 42);
}

TEST(SchedulerTest, TwoComponents) {
  Options options;
  options["components"] = "testcomponent, multiply";
  auto scheduler = ComponentScheduler::create(options, options, nullptr);

  EXPECT_FALSE(options.isSet("answer"));
  scheduler->transform(options);
  ASSERT_TRUE(options.isSet("answer"));
  ASSERT_TRUE(options["answer"] == 42 * 2);
}

TEST(SchedulerTest, SubComponents) {
  Options options;
  options["components"] = "species";
  options["species"]["type"] = "testcomponent, multiply";

  auto scheduler = ComponentScheduler::create(options, options, nullptr);

  EXPECT_FALSE(options.isSet("answer"));
  scheduler->transform(options);
  ASSERT_TRUE(options.isSet("answer"));
  ASSERT_TRUE(options["answer"] == 42 * 2);
}
