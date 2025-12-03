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

struct OrderChecker : public Component {
  OrderChecker(const std::string& name, Options& alloptions, Solver*)
      : Component(getPermissions(name, alloptions)), name(name) {}
  static Permissions getPermissions(const std::string& name, Options& alloptions) {
    return alloptions[name]["permissions"].as<Permissions>();
  }
  static void resetOrderInfo() { execution_order.clear(); }

  std::string name;
  static std::vector<std::string> execution_order;

private:
  void transform_impl(GuardedOptions&) override { execution_order.push_back(name); }
};

std::vector<std::string> OrderChecker::execution_order;

RegisterComponent<TestComponent> registertestcomponent("testcomponent");
RegisterComponent<TestMultiply> registertestcomponent2("multiply");
RegisterComponent<OrderChecker> registercomponentorderchecker("orderchecker");
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

using Parameter = std::pair<Options, std::vector<std::string>>;

class ComponentOrderTest : public testing::TestWithParam<Parameter> {
  void SetUp() override { OrderChecker::resetOrderInfo(); }
};

TEST_P(ComponentOrderTest, Sorted) {
  Options options = GetParam().first.copy();
  auto scheduler = ComponentScheduler::create(options, options, nullptr);
  scheduler->transform(options);
  EXPECT_EQ(OrderChecker::execution_order, GetParam().second);
}

INSTANTIATE_TEST_SUITE_P(
    TopologicalSort, ComponentOrderTest,
    testing::Values(
        Parameter({{"components", ""}}, {}),
        Parameter(
            {{"components", "a"},
             {"a", {{"type", "orderchecker"}, {"permissions", toString(Permissions())}}}},
            {"a"}),
        Parameter({{"components", "a"},
                   {"a",
                    {{"type", "orderchecker"},
                     {"permissions",
                      toString(Permissions({readWrite("1"), readWrite("2")}))}}}},
                  {"a"}),
        Parameter(
            {{"components", "a,b"},
             {"a",
              {{"type", "orderchecker"},
               {"permissions", toString(Permissions({readWrite("1"), readWrite("2")}))}}},
             {"b",
              {{"type", "orderchecker"},
               {"permissions", toString(Permissions({readOnly("1"), readOnly("2")}))}}}},
            {"a", "b"}),
        Parameter({{"components", "a,b"},
                   {"a",
                    {{"type", "orderchecker"},
                     {"permissions", toString(Permissions({readWrite("1"), readWrite("2"),
                                                           readOnly("time")}))}}},
                   {"b",
                    {{"type", "orderchecker"},
                     {"permissions", toString(Permissions({readOnly("1"), readOnly("2"),
                                                           readIfSet("linear"),
                                                           readOnly("units:eV")}))}}}},
                  {"a", "b"}),
        Parameter({{"components", "b,a"},
                   {"b",
                    {{"type", "orderchecker"},
                     {"permissions",
                      toString(Permissions({readOnly("1"), readOnly("2")}))}}},
                   {"a",
                    {{"type", "orderchecker"},
                     {"permissions",
                      toString(Permissions({readWrite("1"), readWrite("2")}))}}}},
                  {"a", "b"}),
        Parameter(
            {{"components", "b,a,c"},
             {"b",
              {{"type", "orderchecker"},
               {"permissions", toString(Permissions({readOnly("1"), readOnly("2")}))}}},
             {"a",
              {{"type", "orderchecker"},
               {"permissions", toString(Permissions({readWrite("1"), readWrite("2")}))}}},
             {"c",
              {{"type", "orderchecker"},
               {"permissions", toString(Permissions({readWrite("2"), readOnly("1")}))}}}},
            {"a", "c", "b"}),
        Parameter({{"components", "a,b"},
                   {"a",
                    {{"type", "orderchecker"},
                     {"permissions", toString(Permissions({readIfSet("1")}))}}},
                   {"b",
                    {{"type", "orderchecker"},
                     {"permissions", toString(Permissions({readWrite("1")}))}}}},
                  {"b", "a"}),
        Parameter({{"components", "a,b,c"},
                   {"a",
                    {{"type", "orderchecker"},
                     {"permissions",
                      toString(Permissions({writeFinal("1"), readIfSet("3")}))}}},
                   {"b",
                    {{"type", "orderchecker"},
                     {"permissions",
                      toString(Permissions({readWrite("1"), readOnly("2")}))}}},
                   {"c",
                    {{"type", "orderchecker"},
                     {"permissions", toString(Permissions({readWrite("1"), readWrite("2"),
                                                           readIfSet("3")}))}}}},
                  {"c", "b", "a"}),
        Parameter({{"components", "a,b,c"},
                   {"a",
                    {{"type", "orderchecker"},
                     {"permissions",
                      toString(Permissions({writeFinal("1:1_1"), readWrite("1:1_2")}))}}},
                   {"b",
                    {{"type", "orderchecker"},
                     {"permissions", toString(Permissions({readWrite("1")}))}}},
                   {"c",
                    {{"type", "orderchecker"},
                     {"permissions", toString(Permissions({readOnly("1:1_1")}))}}}},
                  {"b", "a", "c"}),
        Parameter(
            {{"components", "a,b,c"},
             {"a",
              {{"type", "orderchecker"},
               {"permissions", toString(Permissions({writeFinal("1"), readOnly("2")}))}}},
             {"b",
              {{"type", "orderchecker"},
               {"permissions",
                toString(Permissions({readWrite("2:2_1"), writeFinal("2:2_2")}))}}},
             {"c",
              {{"type", "orderchecker"},
               {"permissions", toString(Permissions({readWrite("1"), readOnly("2:2_1"),
                                                     readIfSet("3")}))}}}},
            {"b", "c", "a"}),
        Parameter({{"components", "a,b,c"},
                   {"a",
                    {{"type", "orderchecker"},
                     {"permissions",
                      toString(Permissions({readOnly("1"), readWrite("1:1_1")}))}}},
                   {"b",
                    {{"type", "orderchecker"},
                     {"permissions",
                      toString(Permissions({readWrite("1:1_1"), readWrite("1:1_2")}))}}},
                   {"c",
                    {{"type", "orderchecker"},
                     {"permissions", toString(Permissions({readOnly("1:1_2")}))}}}},
                  {"b", "a", "c"}),
        Parameter({{"components", "a,b,c"},
                   {"a",
                    {{"type", "orderchecker"},
                     {"permissions", toString(Permissions({writeBoundaryFinal("1")}))}}},
                   {"b",
                    {{"type", "orderchecker"},
                     {"permissions", toString(Permissions({readOnly("1")}))}}},
                   {"c",
                    {{"type", "orderchecker"},
                     {"permissions",
                      toString(Permissions({readWrite("1", Regions::Interior)}))}}}},
                  {"c", "a", "b"}),
        Parameter({{"components", "b,a,c"},
                   {"a",
                    {{"type", "orderchecker"},
                     {"permissions",
                      toString(Permissions({readIfSet("1", Regions::Interior),
                                            readWrite("2", Regions::All)}))}}},
                   {"b",
                    {{"type", "orderchecker"},
                     {"permissions",
                      toString(Permissions({writeFinal("1", Regions::Boundaries),
                                            readIfSet("2", Regions::Interior)}))}}},
                   {"c",
                    {{"type", "orderchecker"},
                     {"permissions",
                      toString(Permissions({readOnly("1", Regions::Boundaries),
                                            readOnly("2", Regions::Boundaries)}))}}}},
                  {"a", "b", "c"}),
        Parameter({{"components", "a,b"},
                   {"a",
                    {{"type", "orderchecker"},
                     {"permissions", toString(Permissions({
                                         readOnly("1"),
                                     }))}}},
                   {"b",
                    {{"type", "orderchecker"},
                     {"permissions", toString(Permissions({writeFinal("1:1_1"),
                                                           readIfSet("1:1_2")}))}}}},
                  {"b", "a"})));

class InvalidComponentOrderTest : public testing::TestWithParam<Options> {
  void SetUp() override { OrderChecker::resetOrderInfo(); }
};

TEST_P(InvalidComponentOrderTest, BadDAG) {
  Options options = GetParam().copy();
  EXPECT_THROW(ComponentScheduler::create(options, options, nullptr), BoutException);
}

INSTANTIATE_TEST_SUITE_P(
    InvalidTopologicalSort, InvalidComponentOrderTest,
    testing::Values(
        // Unsatisfiable dependency
        Options({{"components", "a,b"},
                 {"a",
                  {{"type", "orderchecker"},
                   {"permissions",
                    toString(Permissions({readOnly("1"), readWrite("2")}))}}},
                 {"b",
                  {{"type", "orderchecker"},
                   {"permissions", toString(Permissions({readWrite("3")}))}}}}),
        // Multiple final writes
        Options({{"components", "a,b"},
                 {"a",
                  {{"type", "orderchecker"},
                   {"permissions", toString(Permissions({writeFinal("1")}))}}},
                 {"b",
                  {{"type", "orderchecker"},
                   {"permissions", toString(Permissions({writeFinal("1")}))}}}}),
        // Circular dependency
        Options({{"components", "a,b"},
                 {"a",
                  {{"type", "orderchecker"},
                   {"permissions",
                    toString(Permissions({readOnly("1"), readWrite("2")}))}}},
                 {"b",
                  {{"type", "orderchecker"},
                   {"permissions",
                    toString(Permissions({readWrite("2"), readOnly("1")}))}}}}),
        // Circular dependency from readIfSet
        Options({{"components", "a,b"},
                 {"a",
                  {{"type", "orderchecker"},
                   {"permissions",
                    toString(Permissions({readIfSet("1"), readWrite("2")}))}}},
                 {"b",
                  {{"type", "orderchecker"},
                   {"permissions",
                    toString(Permissions({readOnly("2"), readWrite("1")}))}}}}),
        // Unsatisfiable dependency due to only setting one region
        Options({{"components", "a,b"},
                 {"a",
                  {{"type", "orderchecker"},
                   {"permissions", toString(Permissions({readOnly("1")}))}}},
                 {"b",
                  {{"type", "orderchecker"},
                   {"permissions",
                    toString(Permissions({readWrite("1", Regions::Interior)}))}}}}),
        // Circular dependency on only one region
        Options({{"components", "a,b"},
                 {"a",
                  {{"type", "orderchecker"},
                   {"permissions", toString(Permissions({readOnly("1", Regions::Interior),
                                                         readWrite("2")}))}}},
                 {"b",
                  {{"type", "orderchecker"},
                   {"permissions",
                    toString(Permissions({readWrite("1"), readOnly("2")}))}}}})));
