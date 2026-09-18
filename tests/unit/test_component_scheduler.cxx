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

// A component which, when executed, appends its name to a static
// pulbic member called `execution_order`. It also derives its read
// and write permissions at run-time from the `permissions`
// optino. Taken together, this allows a user to set up different
// dependency chains between OrderChecker component objects and
// confirm that they get executed in the appropriate order.
struct OrderChecker : public NamedComponent<OrderChecker> {
  OrderChecker(const std::string& name, Options& alloptions, Solver*)
      : NamedComponent(name, getPermissions(name, alloptions)) {}
  static Permissions getPermissions(const std::string& name, Options& alloptions) {
    if (alloptions[name].isSet("permissions")) {
      return alloptions[name]["permissions"].as<Permissions>();
    }
    return {};
  }
  static void resetOrderInfo() { execution_order.clear(); }

  static std::vector<std::string> execution_order;

  static constexpr auto type = "orderchecker";

private:
  void transform_impl(GuardedOptions&) override {
    execution_order.push_back(objectName());
  }
};

std::vector<std::string> OrderChecker::execution_order;

RegisterComponent<TestComponent> registertestcomponent;
RegisterComponent<TestMultiply> registertestcomponent2;
RegisterComponent<OrderChecker> registercomponentorderchecker;
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

// Options describing a set of components and the expected order these
// components should be executed.
using Parameter = std::pair<Options, std::vector<std::string>>;

class ComponentOrderTest : public testing::TestWithParam<Parameter> {
  void SetUp() override { OrderChecker::resetOrderInfo(); }
};

// Create a series of components and confirm they are executed in the expected order.
TEST_P(ComponentOrderTest, Sorted) {
  Options options = GetParam().first.copy();
  auto scheduler = ComponentScheduler::create(options, options, nullptr);
  scheduler->transform(options);
  EXPECT_EQ(OrderChecker::execution_order, GetParam().second);
}

INSTANTIATE_TEST_SUITE_P(
    TopologicalSort, ComponentOrderTest,
    testing::Values(
        // Empty set of components
        Parameter({{"components", ""}}, {}),
        // Single component with no permissions
        Parameter(
            {{"components", "a"},
             {"a", {{"type", "orderchecker"}, {"permissions", toString(Permissions())}}}},
            {"a"}),
        // Single component with only read permissions
        Parameter({{"components", "a"},
                   {"a",
                    {{"type", "orderchecker"},
                     {"permissions",
                      toString(Permissions({readWrite("1"), readWrite("2")}))}}}},
                  {"a"}),
        // Two components, one of which reads variables written by the
        // other. The components are listed in the correct order for
        // execution alrady.
        Parameter(
            {{"components", "a,b"},
             {"a",
              {{"type", "orderchecker"},
               {"permissions", toString(Permissions({readWrite("1"), readWrite("2")}))}}},
             {"b",
              {{"type", "orderchecker"},
               {"permissions", toString(Permissions({readOnly("1"), readOnly("2")}))}}}},
            {"a", "b"}),
        // Two components, one of which reads variables written by the
        // other. The second component will try to read an additional
        // variable if it is set elsewhere (but it is not) and one of
        // the default variables set by ComponentScheduler itself. The
        // components are listed in the correct order for execution
        // alrady.
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
        // Two components, one of which reads variables written by the
        // other. The components are not listed in the correct order
        // execution.
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
        // Three components, not listed in the correct order for execution.
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
        // Two components, one of which reads a variable only if
        // another one sets it (which the other one does in this
        // case).
        Parameter({{"components", "a,b"},
                   {"a",
                    {{"type", "orderchecker"},
                     {"permissions", toString(Permissions({readIfSet("1")}))}}},
                   {"b",
                    {{"type", "orderchecker"},
                     {"permissions", toString(Permissions({readWrite("1")}))}}}},
                  {"b", "a"}),
        // Three components, all of which write a particular variable,
        // but one of which must write it last. A second variable
        // determines in what order the remaining two components must
        // be executed. Two of them will try to read a variable if it
        // has been set (which it has not).
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
        // Three components, one of which has permission to write an
        // entire section of the variables. The other two components
        // read and write only particular variables within that
        // section.
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
        // Three components, one of which reads an entire section of
        // variables, while the others read and write only particula
        // variables, either in that section or elsewhere.
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
        // Three components, one of which reads an entire section of
        // variables, but has write permission for a particular
        // variable within that section. The other components can read
        // or write only particular variables.
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
        // The interior and boundary of a variable are written by
        // different components, before being read by a third.
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
        // Only the boundaries are set for one variable, while the
        // entire domain is set for another. The boundaries of both of
        // these are are read by a third component. The first
        // components also read the interior of other variables, if
        // they have been set. At first glance this might appear to
        // cause a circular dependency, but it does not.
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
        // One component reads an entire section, in which one
        // variable is set by a second component. The second component
        // would also read another variable in this section if it were
        // set elsewhere, but it is not.
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

// Check that ComponentScheduler throws an exception when trying to
// sort components with dependency conflicts.
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
