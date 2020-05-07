
#include "gtest/gtest.h"

#include "../../include/component.hxx"

#include <algorithm> // std::any_of

namespace {
struct TestComponent : public Component {
  TestComponent(const std::string&, Options&, Solver *) {}
  void transform(Options &state) { state["answer"] = 42; }
};

RegisterComponent<TestComponent> registertestcomponent("testcomponent");
} // namespace

TEST(ComponentTest, InAvailableList) {
  // Check that the test component is in the list of available components
  auto available = ComponentFactory::getInstance().listAvailable();

  ASSERT_TRUE(std::any_of(
      available.begin(), available.end(),
      [](const std::string &str) { return str == "testcomponent"; }));
}

TEST(ComponentTest, CanCreate) {
  Options options;
  options["mycomponent"]["type"] = "testcomponent"; // This chooses what type of component
  
  auto component = Component::create("mycomponent", options, nullptr);

  EXPECT_FALSE(options.isSet("answer"));
  
  component->transform(options);

  ASSERT_TRUE(options.isSet("answer"));
  ASSERT_TRUE(options["answer"] == 42);
}

TEST(ComponentTest, GetThrowsNoValue) {
  Options option;

  // No value throws
  ASSERT_THROW(get<int>(option), BoutException);
  
  // Compatible value doesn't throw
  option = 42;
  ASSERT_TRUE(option == 42);
}

TEST(ComponentTest, GetThrowsIncompatibleValue) {
  Options option;
  
  option = "hello";
  // Invalid value throws
  ASSERT_THROW(get<int>(option), BoutException);
}

