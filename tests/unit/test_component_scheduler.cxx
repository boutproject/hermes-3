#include "gtest/gtest.h"

#include "../../include/component_scheduler.hxx"

namespace {
struct TestComponent : public Component {
  TestComponent(const std::string&, Options&, Solver*)
      : Component({readWrite("answer")}) {}

private:
  void transform_impl(GuardedOptions& state) override {
    state["answer"].getWritable() = 42;
  }
};
  

struct TestMultiply : public Component {
  TestMultiply(const std::string&, Options&, Solver*)
      : Component({writeFinal("answer")}) {}

private:
  void transform_impl(GuardedOptions& state) override {
    // Note: Using set<>() and get<>() for quicker access, avoiding printing
    //       getNonFinal needs to be used because we set the value afterwards
    set(state["answer"],
        getNonFinal<int>(state["answer"]) * 2);
  }
};

RegisterComponent<TestComponent> registertestcomponent("testcomponent");
RegisterComponent<TestMultiply> registertestcomponent2("multiply");
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

