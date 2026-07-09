#include <fmt/base.h>

#include "gtest/gtest.h"

#include "../../include/component.hxx"
#include "../../include/reaction.hxx"
#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include <bout/field_factory.hxx> // For generating functions

#include <algorithm> // std::any_of

namespace {
struct TestComponent : public NamedComponent<TestComponent> {
  TestComponent(const std::string name, Options&, Solver*)
      : NamedComponent(name, {readWrite("answer")}) {}

  static constexpr auto type = "testcomponent";

private:
  void transform_impl(GuardedOptions& state) override {
    state["answer"].getWritable() = 42;
  }
};

RegisterComponent<TestComponent> registertestcomponent;
} // namespace

TEST(ComponentTest, InAvailableList) {
  // Check that the test component is in the list of available components
  auto available = ComponentFactory::getInstance().listAvailable();

  ASSERT_TRUE(std::any_of(available.begin(), available.end(),
                          [](const std::string& str) { return str == "testcomponent"; }));
}

TEST(ComponentTest, CanCreate) {
  Options options;
  auto component = Component::create("testcomponent", "species", options, nullptr);

  EXPECT_FALSE(options.isSet("answer"));

  component->transform(options);

  ASSERT_TRUE(options.isSet("answer"));
  ASSERT_TRUE(options["answer"] == 42);
}

TEST(ComponentTest, ObjectName) {
  Options options;
  auto component = Component::create("testcomponent", "some_name", options, nullptr);
  ASSERT_EQ(component->objectName(), "some_name");
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

#if CHECKLEVEL >= 1
TEST(ComponentTest, SetAfterGetThrows) {
  Options option;

  option = 42;

  ASSERT_EQ(get<int>(option), 42);

  // Setting after get should fail
  ASSERT_THROW(set<int>(option, 3), BoutException);
}
#endif

TEST(ComponentTest, SetAfterGetNonFinal) {
  Options option;

  option = 42;

  ASSERT_EQ(getNonFinal<int>(option), 42);

  set<int>(option, 3); // Doesn't throw

  ASSERT_EQ(getNonFinal<int>(option), 3);
}

#if CHECKLEVEL >= 1
TEST(ComponentTest, SetBoundaryAfterGetThrows) {
  Options option;

  option = 42;

  ASSERT_EQ(get<int>(option), 42);

  // Setting after get should fail because get indicates an assumption
  // that all values are final including boundary cells.
  ASSERT_THROW(setBoundary<int>(option, 3), BoutException);
}
#endif

TEST(ComponentTest, SetBoundaryAfterGetNoBoundary) {
  Options option;

  option = 42;

  ASSERT_EQ(getNoBoundary<int>(option), 42);

  setBoundary<int>(option, 3); // ok because boundary not assumed final

  ASSERT_EQ(getNonFinal<int>(option), 3);
}

TEST(ComponentTest, IsSetFinalStaysFalse) {
  Options option;

  ASSERT_EQ(isSetFinal(option["test"]), false);
  // Shouldn't change if called again
  ASSERT_EQ(isSetFinal(option["test"]), false);
}

TEST(ComponentTest, GetAfterIsSetFinal) {
  Options option;
  option["test"] = 1;

  ASSERT_EQ(isSetFinal(option["test"]), true);
  // Can get the value
  ASSERT_EQ(get<int>(option["test"]), 1);
}

#if CHECKLEVEL >= 1
TEST(ComponentTest, SetAfterIsSetFinal) {
  Options option;

  ASSERT_EQ(isSetFinal(option["test"]), false);
  // Can't now set the value
  ASSERT_THROW(set<int>(option["test"], 3), BoutException);
}
#endif

TEST(ComponentTest, Formatting) {
  Options options;
  auto component_samename =
      Component::create("testcomponent", "testcomponent", options, nullptr);
  EXPECT_EQ(fmt::format("{}", *component_samename), "testcomponent");
  EXPECT_EQ(fmt::format("{:~n}", *component_samename), "testcomponent");
  EXPECT_EQ(fmt::format("{:~t}", *component_samename), "testcomponent");
  EXPECT_EQ(fmt::format("{:~n~t}", *component_samename), "");
  EXPECT_EQ(fmt::format("{:T}", *component_samename), "testcomponent (testcomponent)");

  auto component_diffname =
      Component::create("testcomponent", "object_name", options, nullptr);
  EXPECT_EQ(fmt::format("{}", *component_diffname), "object_name (testcomponent)");
  EXPECT_EQ(fmt::format("{:~n}", *component_diffname), "testcomponent");
  EXPECT_EQ(fmt::format("{:~t}", *component_diffname), "object_name");
  EXPECT_EQ(fmt::format("{:~n~t}", *component_diffname), "");
  EXPECT_EQ(fmt::format("{:T}", *component_diffname), "object_name (testcomponent)");

  EXPECT_EQ(fmt::format(fmt::runtime("{:xT}"), *component_diffname),
            "object_name (testcomponent)");
  EXPECT_THROW((void)fmt::format(fmt::runtime("{:~}"), *component_diffname),
               fmt::format_error);
}

using ComponentCreationTest = FakeMeshFixture;
/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

struct ConcreteComponentTests : public FakeMeshFixture,
                                public testing::WithParamInterface<std::string> {
  static const Options base_options;
  static const Options required_params;
  static constexpr auto objname = "object_name";

  std::string typname;
  Options options;
  FakeSolver solver;

  ConcreteComponentTests()
      : FakeMeshFixture(), typname(GetParam()), options(base_options.copy()) {
    Options::root()["mesh:paralleltransform:type"] = "identity";
    static_cast<FakeMesh*>(bout::globals::mesh)
        ->setGridDataSource(new FakeGridDataSource{
            {{"Rxy", FieldFactory::get()->create2D("1 + x", Options::getRoot(),
                                                   bout::globals::mesh)},
             {"Zxy", FieldFactory::get()->create2D("y", Options::getRoot(),
                                                   bout::globals::mesh)},
             {"hthe", 1.0},
             {"Bpxy", 1.0},
             {"Bxy", 1.0},
             {"external_apar", 1.0}}});
    if (required_params.isSection(typname)) {
      options[objname] = required_params[typname].copy();
    }
    options[objname]["type"] = typname;
    options[objname]["charge"] = 1.;
    options[objname]["AA"] = 1.;
  }

  ~ConcreteComponentTests() override { hermes::ReactionBase::reset_instance_counter(); }
};

const Options ConcreteComponentTests::base_options{
    {"units",
     {{"eV", 1.0},
      {"meters", 1.0},
      {"seconds", 1.0},
      {"inv_meters_cubed", 1.0},
      {"Tesla", 1.0}}},
    {"mesh", {{"length", 1}, {"paralleltransform", {{"type", "identity"}}}}}};
const Options ConcreteComponentTests::required_params{
    {"detachment_controller",
     {{"detachment_front_setpoint", 0.},
      {"neutral_species", "d"},
      {"actuator", "particles"},
      {"species_list", ""},
      {"scaling_factors_list", ""}}},
    {"fixed_density", {{"density", 1.}}},
    {"fixed_fraction_ions", {{"fractions", "h+@1."}}},
    {"fixed_temperature", {{"temperature", 1.}}},
    {"fixed_velocity", {{"velocity", 1.}}},
    {"isothermal", {{"temperature", 1.}}},
    {"neutral_parallel_diffusion", {{"dneut", 1.}}},
    {"recycling", {{"species", ""}}},
    {"sheath_closure", {{"connection_length", 0.}}},
    {"temperature_feedback",
     {{"temperature_setpoint", 1.},
      {"species_for_temperature_feedback", ""},
      {"scaling_factors_for_temperature_feedback", ""}}},
    {"transform", {{"transforms", ""}}},
    {"upstream_density_feedback", {{"density_upstream", 1.}}},
    {"zero_current", {{"charge", 1.}}},
};

// Check the result of the Component::typeName() function is the
// same as the type name used to create the component.
TEST_P(ConcreteComponentTests, CheckComponentTypeName) {
  auto component =
      ComponentFactory::getInstance().create(typname, objname, options, &solver);
  EXPECT_EQ(component->typeName(), typname);
}

INSTANTIATE_TEST_SUITE_P(
    AllRegisteredComponents, ConcreteComponentTests,
    testing::ValuesIn(ComponentFactory::getInstance().listAvailable()));
