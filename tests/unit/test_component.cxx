
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
  auto component = Component::create("testcomponent", "species", options, nullptr);

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

TEST(ComponentTest, SetInteger) {
  Options option;

  set<int>(option, 3);

  ASSERT_EQ(getNonFinal<int>(option), 3);
}

TEST(ComponentTest, SetAfterGetThrows) {
  Options option;

  option = 42;

  ASSERT_EQ(get<int>(option), 42);

  // Setting after get should fail
  ASSERT_THROW(set<int>(option, 3), BoutException);
}

TEST(ComponentTest, SetAfterGetNonFinal) {
  Options option;

  option = 42;

  ASSERT_EQ(getNonFinal<int>(option), 42);

  set<int>(option, 3); // Doesn't throw

  ASSERT_EQ(getNonFinal<int>(option), 3);
}

TEST(ComponentTest, SetBoundaryAfterGetThrows) {
  Options option;

  option = 42;

  ASSERT_EQ(get<int>(option), 42);

  // Setting after get should fail because get indicates an assumption
  // that all values are final including boundary cells.
  ASSERT_THROW(setBoundary<int>(option, 3), BoutException);
}

TEST(ComponentTest, SetBoundaryAfterGetNoBoundary) {
  Options option;

  option = 42;

  ASSERT_EQ(getNoBoundary<int>(option), 42);

  setBoundary<int>(option, 3); // ok because boundary not assumed final

  ASSERT_EQ(getNonFinal<int>(option), 3);
}
