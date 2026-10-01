#include <fmt/base.h>

#include "gtest/gtest.h"

#include "../../include/component.hxx"
#include "../../include/reaction.hxx"
#include "fake_mesh_fixture.hxx"
#include "fake_solver.hxx"
#include "guarded_options.hxx"
#include "permissions.hxx"
#include <bout/boutexception.hxx>
#include <bout/field_factory.hxx> // For generating functions

#include <algorithm> // std::any_of

struct ComponentTest : public FakeMeshFixture {
  ComponentTest() : FakeMeshFixture() {
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
  }
};

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

TEST_F(ComponentTest, InAvailableList) {
  // Check that the test component is in the list of available components
  auto available = ComponentFactory::getInstance().listAvailable();

  ASSERT_TRUE(std::any_of(available.begin(), available.end(),
                          [](const std::string& str) { return str == "testcomponent"; }));
}

TEST_F(ComponentTest, CanCreate) {
  Options options;
  auto component = Component::create("testcomponent", "species", options, nullptr);

  EXPECT_FALSE(options.isSet("answer"));

  component->transform(options);

  ASSERT_TRUE(options.isSet("answer"));
  ASSERT_TRUE(options["answer"] == 42);
}

TEST_F(ComponentTest, ObjectName) {
  Options options;
  auto component = Component::create("testcomponent", "some_name", options, nullptr);
  ASSERT_EQ(component->objectName(), "some_name");
}

TEST_F(ComponentTest, GetThrowsNoValue) {
  Options option;

  // No value throws
  ASSERT_THROW(get<int>(option), BoutException);

  // Compatible value doesn't throw
  option = 42;
  ASSERT_TRUE(option == 42);
}

TEST_F(ComponentTest, GuardedGet) {
  Options options{{"a", 1}, {"b", 2}, {"c", 3}, {"d", 4}};
  Permissions perms{readIfSet("a"), readWrite("b", Regions::Interior),
                    readOnly("c", Regions::Boundaries), readIfSet("z")};
  GuardedOptions gopts{&options, &perms};

  ASSERT_EQ(get<int>(gopts["a"]), 1);
  // Don't have permission for entire domain
#if CHECKLEVEL >= 1
  ASSERT_THROW(get<int>(gopts["b"]), BoutException);
  ASSERT_THROW(get<int>(gopts["c"]), BoutException);
  ASSERT_THROW(get<int>(gopts["d"]), BoutException);
#endif
  ASSERT_THROW(get<int>(gopts["z"]), BoutException);
}

#if CHECKLEVEL >= 1
TEST_F(ComponentTest, SetNaN) {
  Options option;
  EXPECT_THROW(set(option, Field3D{BoutNaN, bout::globals::mesh}), BoutException);
}
#endif

TEST_F(ComponentTest, GetThrowsIncompatibleValue) {
  Options option;

  option = "hello";
  // Invalid value throws
  ASSERT_THROW(get<int>(option), BoutException);
}

TEST_F(ComponentTest, SetInteger) {
  Options option;

  set<int>(option, 3);

  ASSERT_EQ(getNonFinal<int>(option), 3);
}

TEST_F(ComponentTest, GuardedSetInteger) {
  Options options;
  Permissions perms{readWrite("a"), readWrite("b", Regions::Interior),
                    readWrite("c", Regions::Boundaries)};
  GuardedOptions gopts{&options, &perms};

  set<int>(gopts["a"], 1);
  ASSERT_EQ(get<int>(gopts["a"]), 1);

  // Don't have permission for entire domain
#if CHECKLEVEL >= 1
  ASSERT_THROW(set<int>(gopts["b"], 2), BoutException);
  ASSERT_THROW(set<int>(gopts["c"], 3), BoutException);
  ASSERT_THROW(set<int>(gopts["d"], 4), BoutException);
#endif
}

#if CHECKLEVEL >= 1
TEST_F(ComponentTest, SetAfterGetThrows) {
  Options option;

  option = 42;

  ASSERT_EQ(get<int>(option), 42);

  // Setting after get should fail
  ASSERT_THROW(set<int>(option, 3), BoutException);
}

// Check it is an exception to write a variable after it has been read.
TEST_F(ComponentTest, GuardedSetAfterGetThrows) {
  Options option;
  Permissions perms{readWrite("test")};
  GuardedOptions gopts{&option, &perms};

  option["test"] = 42;

  ASSERT_EQ(get<int>(gopts["test"]), 42);

  // Setting after get should fail
  ASSERT_THROW(set<int>(gopts["test"], 3), BoutException);
}
#endif

TEST_F(ComponentTest, SetAfterGetNonFinal) {
  Options option;

  option = 42;

  ASSERT_EQ(getNonFinal<int>(option), 42);

  set<int>(option, 3); // Doesn't throw

  ASSERT_EQ(getNonFinal<int>(option), 3);
}

// Check it is permitted to set a variable after it has been read using getNonFinal.
TEST_F(ComponentTest, GuardedSetAfterGetNonFinal) {
  Options option;
  Permissions perms{readWrite("a")};

  option["a"] = 42;
  GuardedOptions gopts{&option, &perms};

  ASSERT_EQ(getNonFinal<int>(gopts["a"]), 42);

  set<int>(gopts["a"], 3); // Doesn't throw

  ASSERT_EQ(getNonFinal<int>(gopts["a"]), 3);
}

#if CHECKLEVEL >= 1
TEST_F(ComponentTest, SetBoundaryAfterGetThrows) {
  Options option;

  option = 42;

  ASSERT_EQ(get<int>(option), 42);

  // Setting after get should fail because get indicates an assumption
  // that all values are final including boundary cells.
  ASSERT_THROW(setBoundary<int>(option, 3), BoutException);
}

// Check can not write boundary after it has been read, even if the
// boundary has write-permissions.
TEST_F(ComponentTest, GuardedSetBoundaryAfterGetThrows) {
  Options option;
  Permissions perms{writeBoundaryReadInteriorIfSet(".")};

  option = 42;
  GuardedOptions gopt{&option, &perms};

  ASSERT_EQ(get<int>(option), 42);

  // Setting after get should fail because get indicates an assumption
  // that all values are final including boundary cells.
  ASSERT_THROW(setBoundary<int>(option, 3), BoutException);
}

// Check an exception is thrown if you try to set the interior after
// reading the entire field
TEST_F(ComponentTest, SetNoBoundaryAfterGetThrows) {
  Options option;

  option = 42;

  ASSERT_EQ(get<int>(option), 42);

  // Setting after get should fail because get indicates an assumption
  // that all values are final including interior cells.
  ASSERT_THROW(setNoBoundary<int>(option, 3), BoutException);
}

// Check an exception is thrown if you try to set the interior after
// reading the entire field, even if the interior has read permission
TEST_F(ComponentTest, GuardedSetNoBoundaryAfterGetThrows) {
  Options option;
  Permissions perms{
      {"a",
       {Regions::Nowhere, Regions::Boundaries, Regions::Interior, Regions::Nowhere}}};
  GuardedOptions gopts{&option, &perms};

  option["a"] = 42;

  ASSERT_EQ(get<int>(gopts["a"]), 42);

  // Setting after get should fail because get indicates an assumption
  // that all values are final including interior cells.
  ASSERT_THROW(setNoBoundary<int>(gopts["a"], 3), BoutException);
}
#endif

TEST_F(ComponentTest, SetBoundaryAfterGetNoBoundary) {
  Options option;

  option = 42;

  ASSERT_EQ(getNoBoundary<int>(option), 42);

  setBoundary<int>(option, 3); // ok because boundary not assumed final

  ASSERT_EQ(getNonFinal<int>(option), 3);
}

// Check still allowed to set boundaries (if have permission) even after reading interior
TEST_F(ComponentTest, GuardedSetBoundaryAfterGetNoBoundary) {
  Options options{{"a", 42}, {"b", 43}, {"c", 44}, {"d", 45}};
  Permissions perms{
      readWrite("a"),
      {"b", {Regions::Nowhere, Regions::Interior, Regions::Boundaries, Regions::Nowhere}},
      {"c", {Regions::Nowhere, Regions::Boundaries, Regions::Interior, Regions::Nowhere}},
      readOnly("d")};
  GuardedOptions gopts{&options, &perms};

  ASSERT_EQ(getNoBoundary<int>(gopts["a"]), 42);
  ASSERT_EQ(getNoBoundary<int>(gopts["b"]), 43);
  ASSERT_EQ(getNoBoundary<int>(gopts["c"]), 44);
  ASSERT_EQ(getNoBoundary<int>(gopts["d"]), 45);

  setBoundary<int>(gopts["a"], 3); // ok because boundary not assumed final
  setBoundary<int>(gopts["b"], 4); // ok because boundary not assumed final

#if CHECKLEVEL >= 1
  ASSERT_THROW(setBoundary<int>(gopts["c"], 5),
               BoutException); // Don't have permission to write boundary
  ASSERT_THROW(setBoundary<int>(gopts["d"], 6),
               BoutException); // Don't have permission to write boundary
#endif

  ASSERT_EQ(getNonFinal<int>(gopts["a"]), 3);
  ASSERT_EQ(getNonFinal<int>(gopts["b"]), 4);
}

// Check still allowed to set interior even after having read the boundary
TEST_F(ComponentTest, SetNoBoundaryAfterGetBoundary) {
  Options option;

  option = 42;

  ASSERT_EQ(getBoundary<int>(option), 42);

  setNoBoundary<int>(option, 3); // ok because domain not assumed final

  ASSERT_EQ(getNonFinal<int>(option), 3);
}

// Check still allowed to set interior (if have permission) even after reading boundaries
TEST_F(ComponentTest, GuardedSetNoBoundaryAfterGetBoundary) {
  Options options{{"a", 42}, {"b", 43}, {"c", 44}, {"d", 45}};
  Permissions perms{
      readWrite("a"),
      {"b", {Regions::Nowhere, Regions::Interior, Regions::Boundaries, Regions::Nowhere}},
      {"c", {Regions::Nowhere, Regions::Boundaries, Regions::Interior, Regions::Nowhere}},
      readOnly("d")};
  GuardedOptions gopts{&options, &perms};

  ASSERT_EQ(getBoundary<int>(gopts["a"]), 42);
  ASSERT_EQ(getBoundary<int>(gopts["b"]), 43);
  ASSERT_EQ(getBoundary<int>(gopts["c"]), 44);
  ASSERT_EQ(getBoundary<int>(gopts["d"]), 45);

  setNoBoundary<int>(gopts["a"], 3); // ok because boundary not assumed final
  setNoBoundary<int>(gopts["c"], 4); // ok because boundary not assumed final

#if CHECKLEVEL >= 1
  ASSERT_THROW(setNoBoundary<int>(gopts["b"], 5),
               BoutException); // Don't have permission to write interior
  ASSERT_THROW(setNoBoundary<int>(gopts["d"], 6),
               BoutException); // Don't have permission to write interior
#endif

  ASSERT_EQ(getNonFinal<int>(gopts["a"]), 3);
  ASSERT_EQ(getNonFinal<int>(gopts["c"]), 4);
}

TEST_F(ComponentTest, IsSetFinalStaysFalse) {
  Options option;

  ASSERT_EQ(isSetFinal(option["test"]), false);
  // Shouldn't change if called again
  ASSERT_EQ(isSetFinal(option["test"]), false);
}

// Calling isSetFinal shouldn't mark as variable as being final. It
// can be used regardless of permissions; we often need to check if
// something is set to know whether or not to read it.
TEST_F(ComponentTest, GuardedIsSetFinalStaysFalse) {
  Options option;
  Permissions perms{readOnly("a"), readOnly("b", Regions::Boundaries),
                    readOnly("c", Regions::Interior)};
  GuardedOptions gopts{&option, &perms};

  ASSERT_EQ(isSetFinal(gopts["a"]), false);
  // Shouldn't change if called again
  ASSERT_EQ(isSetFinal(gopts["a"]), false);

  ASSERT_EQ(isSetFinal(gopts["b"]), false);
  // Shouldn't change if called again
  ASSERT_EQ(isSetFinal(gopts["b"]), false);

  ASSERT_EQ(isSetFinal(gopts["c"]), false);
  // Shouldn't change if called again
  ASSERT_EQ(isSetFinal(gopts["c"]), false);

  ASSERT_EQ(isSetFinal(gopts["d"]), false);
  // Shouldn't change if called again
  ASSERT_EQ(isSetFinal(gopts["d"]), false);
}

// Confirm that checking whether the interior has been set does not
// mark the interior as set.
TEST_F(ComponentTest, IsSetFinalNoBoundaryStaysFalse) {
  Options option;

  ASSERT_EQ(isSetFinalNoBoundary(option["test"]), false);
  // Shouldn't change if called again
  ASSERT_EQ(isSetFinalNoBoundary(option["test"]), false);
}

// Confirm that checking whether the interior has been set does not
// mark the interior as set.
TEST_F(ComponentTest, GuardedIsSetFinalNoBoundaryStaysFalse) {
  Options option;
  Permissions perms{readOnly("a"), readOnly("b", Regions::Boundaries),
                    readOnly("c", Regions::Interior)};
  GuardedOptions gopts{&option, &perms};

  ASSERT_EQ(isSetFinalNoBoundary(gopts["a"]), false);
  // Shouldn't change if called again
  ASSERT_EQ(isSetFinalNoBoundary(gopts["a"]), false);

  ASSERT_EQ(isSetFinalNoBoundary(gopts["b"]), false);
  // Shouldn't change if called again
  ASSERT_EQ(isSetFinalNoBoundary(gopts["b"]), false);

  ASSERT_EQ(isSetFinalNoBoundary(gopts["c"]), false);
  // Shouldn't change if called again
  ASSERT_EQ(isSetFinalNoBoundary(gopts["c"]), false);

  ASSERT_EQ(isSetFinalNoBoundary(gopts["d"]), false);
  // Shouldn't change if called again
  ASSERT_EQ(isSetFinalNoBoundary(gopts["d"]), false);
}

// Confirm that checking whether the boundary has been set does not
// mark the boundary as set.
TEST_F(ComponentTest, IsSetFinalBoundaryStaysFalse) {
  Options option;

  ASSERT_EQ(isSetFinalBoundary(option["test"]), false);
  // Shouldn't change if called again
  ASSERT_EQ(isSetFinalBoundary(option["test"]), false);
}

// Confirm that checking whether the boundary has been set does not
// mark the boundary as set.
TEST_F(ComponentTest, GuardedIsSetFinalBoundaryStaysFalse) {
  Options option;
  Permissions perms{readOnly("a"), readOnly("b", Regions::Boundaries),
                    readOnly("c", Regions::Interior)};
  GuardedOptions gopts{&option, &perms};

  ASSERT_EQ(isSetFinalNoBoundary(gopts["a"]), false);
  // Shouldn't change if called again
  ASSERT_EQ(isSetFinalNoBoundary(gopts["a"]), false);

  ASSERT_EQ(isSetFinalNoBoundary(gopts["b"]), false);
  // Shouldn't change if called again
  ASSERT_EQ(isSetFinalNoBoundary(gopts["b"]), false);

  ASSERT_EQ(isSetFinalNoBoundary(gopts["c"]), false);
  // Shouldn't change if called again
  ASSERT_EQ(isSetFinalNoBoundary(gopts["c"]), false);

  ASSERT_EQ(isSetFinalNoBoundary(gopts["d"]), false);
  // Shouldn't change if called again
  ASSERT_EQ(isSetFinalNoBoundary(gopts["d"]), false);
}

// Confirm you can read a value after having confirmed it's been set.
TEST_F(ComponentTest, GetAfterIsSetFinal) {
  Options option;
  option["test"] = 1;

  ASSERT_EQ(isSetFinal(option["test"]), true);
  // Can get the value
  ASSERT_EQ(get<int>(option["test"]), 1);
}

// Confirm you can read a value after having confirmed it's been set
// and that read permissions are respected.
TEST_F(ComponentTest, GuardedGetAfterIsSetFinal) {
  Options option;
  Permissions perms{readOnly("a")};
  GuardedOptions gopts{&option, &perms};
  option["a"] = 1;

  ASSERT_EQ(isSetFinal(gopts["a"]), true);
  // Can get the value
  ASSERT_EQ(get<int>(gopts["a"]), 1);
}

// Confirm you can read an interior value after having confirmed the
// interior has been been set.
TEST_F(ComponentTest, GetAfterIsSetFinalNoBoundary) {
  Options option;
  option["test"] = 1;

  ASSERT_EQ(isSetFinalNoBoundary(option["test"]), true);
  // Can get the value
  ASSERT_EQ(getNoBoundary<int>(option["test"]), 1);
}

// Confirm you can read an interior value after having confirmed the
// interior has been been set.
TEST_F(ComponentTest, GuardedGetAfterIsSetFinalNoBoundary) {
  Options option;
  Permissions perms{readOnly("a"), readOnly("b", Regions::Interior)};
  GuardedOptions gopts{&option, &perms};
  option["a"] = 1;
  option["b"] = 1;

  ASSERT_EQ(isSetFinalNoBoundary(gopts["a"]), true);
  // Can get the value
  ASSERT_EQ(getNoBoundary<int>(gopts["a"]), 1);

  ASSERT_EQ(isSetFinalNoBoundary(gopts["b"]), true);
  // Can get the value
  ASSERT_EQ(getNoBoundary<int>(gopts["b"]), 1);
}

// Confirm you can read a boundary value after having confirmed the
// boundary has been been set.
TEST_F(ComponentTest, GetAfterIsSetFinalBoundary) {
  Options option;
  option["test"] = 1;

  ASSERT_EQ(isSetFinalBoundary(option["test"]), true);
  // Can get the value
  ASSERT_EQ(getBoundary<int>(option["test"]), 1);
}

// Confirm you can read a boundary value after having confirmed the
// boundary has been been set.
TEST_F(ComponentTest, GuardedGetAfterIsSetFinalBoundary) {
  Options option;
  Permissions perms{readOnly("a"), readOnly("b", Regions::Boundaries)};
  GuardedOptions gopts{&option, &perms};
  option["a"] = 1;
  option["b"] = 1;

  ASSERT_EQ(isSetFinalBoundary(gopts["a"]), true);
  // Can get the value
  ASSERT_EQ(getBoundary<int>(gopts["a"]), 1);

  ASSERT_EQ(isSetFinalBoundary(gopts["b"]), true);
  // Can get the value
  ASSERT_EQ(getBoundary<int>(gopts["b"]), 1);
}

#if CHECKLEVEL >= 1
// Confirm you can not set the value on any part of the domain after
// checking whether it has been set.
TEST_F(ComponentTest, SetAfterIsSetFinal) {
  Options option;

  ASSERT_EQ(isSetFinal(option["test"]), false);
  // Can't now set the value
  ASSERT_THROW(set<int>(option["test"], 3), BoutException);
  ASSERT_THROW(setNoBoundary<int>(option["test"], 3), BoutException);
  ASSERT_THROW(setBoundary<int>(option["test"], 3), BoutException);
}

// Confirm you can not set the value on any part of the domain after
// checking whether it has been set.
TEST_F(ComponentTest, GuardedSetAfterIsSetFinal) {
  Options option;
  Permissions perms{readWrite("a"), readOnly("b"), readOnly("c", Regions::Interior),
                    readOnly("d", Regions::Boundaries), readIfSet("e")};
  GuardedOptions gopts{&option, &perms};

  ASSERT_EQ(isSetFinal(gopts["a"]), false);
  // Can't now set the value
  ASSERT_THROW(set<int>(gopts["a"], 3), BoutException);
  ASSERT_THROW(setNoBoundary<int>(gopts["a"], 3), BoutException);
  ASSERT_THROW(setBoundary<int>(gopts["a"], 3), BoutException);

  ASSERT_EQ(isSetFinal(gopts["b"]), false);
  // Can't now set the value
  ASSERT_THROW(set<int>(option["b"], 3), BoutException);
  ASSERT_THROW(setNoBoundary<int>(option["b"], 3), BoutException);
  ASSERT_THROW(setBoundary<int>(option["b"], 3), BoutException);

  ASSERT_EQ(isSetFinal(gopts["c"]), false);
  // Can't now set the value
  ASSERT_THROW(set<int>(option["c"], 3), BoutException);
  ASSERT_THROW(setNoBoundary<int>(option["c"], 3), BoutException);
  ASSERT_THROW(setBoundary<int>(option["c"], 3), BoutException);

  ASSERT_EQ(isSetFinal(gopts["d"]), false);
  // Can't now set the value
  ASSERT_THROW(set<int>(option["d"], 3), BoutException);
  ASSERT_THROW(setNoBoundary<int>(option["d"], 3), BoutException);
  ASSERT_THROW(setBoundary<int>(option["d"], 3), BoutException);

  ASSERT_EQ(isSetFinal(gopts["e"]), false);
  // Can still set the value, as isSetFinal doesn't mark variables it
  // doesn't have permission to read.
  set<int>(option["e"], 3);
  ASSERT_EQ(get<int>(option["e"]), 3);

  ASSERT_EQ(isSetFinal(gopts["f"]), false);
  // Can still set the value, as isSetFinal doesn't mark variables it
  // doesn't have permission to read.
  set<int>(option["f"], 3);
  ASSERT_EQ(get<int>(option["f"]), 3);
}

// Confirm you can not set the value in the interior after checking
// whether it has been set there, but you can still set the boundary.
TEST_F(ComponentTest, SetAfterIsSetFinalNoBoundary) {
  Options option;

  ASSERT_EQ(isSetFinalNoBoundary(option["test"]), false);
  // Can't now set the value in the domain
  ASSERT_THROW(set<int>(option["test"], 3), BoutException);
  ASSERT_THROW(setNoBoundary<int>(option["test"], 3), BoutException);
  // Can still set the value for the bounds
  setBoundary<int>(option["test"], 1);
  ASSERT_EQ(getBoundary<int>(option["test"]), 1);
}

// Confirm you can not set the value in the interior after checking
// whether it has been set there, but you can still set the boundary.
TEST_F(ComponentTest, GuardedSetAfterIsSetFinalNoBoundary) {
  Options option;
  Permissions perms{readWrite("a"), readOnly("b"), readOnly("c", Regions::Interior),
                    readOnly("d", Regions::Boundaries), readIfSet("e")};
  GuardedOptions gopts{&option, &perms};

  ASSERT_EQ(isSetFinalNoBoundary(gopts["a"]), false);
  // Can't now set the value in the domain
  ASSERT_THROW(set<int>(gopts["a"], 3), BoutException);
  ASSERT_THROW(setNoBoundary<int>(gopts["a"], 3), BoutException);
  // Can still set the value for the bounds
  setBoundary<int>(gopts["a"], 1);
  ASSERT_EQ(getBoundary<int>(gopts["a"]), 1);

  ASSERT_EQ(isSetFinalNoBoundary(gopts["b"]), false);
  // Can't now set the value in the domain
  ASSERT_THROW(set<int>(option["b"], 3), BoutException);
  ASSERT_THROW(setNoBoundary<int>(option["b"], 3), BoutException);
  // Can still set the value for the bounds
  setBoundary<int>(option["b"], 1);
  ASSERT_EQ(getBoundary<int>(option["b"]), 1);

  ASSERT_EQ(isSetFinalNoBoundary(gopts["c"]), false);
  // Can't now set the value in the domain
  ASSERT_THROW(set<int>(option["c"], 3), BoutException);
  ASSERT_THROW(setNoBoundary<int>(option["c"], 3), BoutException);
  // Can still set the value for the bounds
  setBoundary<int>(option["c"], 1);
  ASSERT_EQ(getBoundary<int>(option["c"]), 1);

  ASSERT_EQ(isSetFinalNoBoundary(gopts["d"]), false);
  // Can still set the value, as isSetFinalNoBoundary doesn't mark variables it
  // doesn't have permission to read.
  set<int>(option["d"], 1);
  ASSERT_EQ(get<int>(option["d"]), 1);

  ASSERT_EQ(isSetFinalNoBoundary(gopts["e"]), false);
  // Can still set the value, as isSetFinalNoBoundary doesn't mark variables it
  // doesn't have permission to read.
  set<int>(option["e"], 1);
  ASSERT_EQ(get<int>(option["e"]), 1);

  ASSERT_EQ(isSetFinalNoBoundary(gopts["f"]), false);
  // Can still set the value, as isSetFinalNoBoundary doesn't mark variables it
  // doesn't have permission to read.
  set<int>(option["f"], 1);
  ASSERT_EQ(get<int>(option["f"]), 1);
}

// Confirm you can not set the value in the boundary after checking
// whether it has been set there, but you can still set the interior.
TEST_F(ComponentTest, SetAfterIsSetFinalBoundary) {
  Options option;

  ASSERT_EQ(isSetFinalBoundary(option["test"]), false);
  // Can't now set the value in the bounds
  ASSERT_THROW(set<int>(option["test"], 3), BoutException);
  ASSERT_THROW(setBoundary<int>(option["test"], 3), BoutException);
  // Can still set the value for the domain
  setNoBoundary<int>(option["test"], 1);
  ASSERT_EQ(getNoBoundary<int>(option["test"]), 1);
}

// Confirm you can not set the value in the boundary after checking
// whether it has been set there, but you can still set the interior.
TEST_F(ComponentTest, GuardedSetAfterIsSetFinalBoundary) {
  Options option;
  Permissions perms{readWrite("a"), readOnly("b"), readOnly("c", Regions::Interior),
                    readOnly("d", Regions::Boundaries), readIfSet("e")};
  GuardedOptions gopts{&option, &perms};

  ASSERT_EQ(isSetFinalBoundary(gopts["a"]), false);
  // Can't now set the value in the boundary
  ASSERT_THROW(set<int>(gopts["a"], 3), BoutException);
  ASSERT_THROW(setBoundary<int>(gopts["a"], 3), BoutException);
  // Can still set the value for the domain
  setNoBoundary<int>(gopts["a"], 1);
  ASSERT_EQ(getNoBoundary<int>(gopts["a"]), 1);

  ASSERT_EQ(isSetFinalBoundary(gopts["b"]), false);
  // Can't now set the value in the boundary
  ASSERT_THROW(set<int>(option["b"], 3), BoutException);
  ASSERT_THROW(setBoundary<int>(option["b"], 3), BoutException);
  // Can still set the value for the domain
  setNoBoundary<int>(option["b"], 1);
  ASSERT_EQ(getNoBoundary<int>(option["b"]), 1);

  ASSERT_EQ(isSetFinalBoundary(gopts["c"]), false);
  // Can still set the value, as isSetFinalBoundary doesn't mark variables it
  // doesn't have permission to read.
  set<int>(option["c"], 1);
  ASSERT_EQ(get<int>(option["c"]), 1);

  ASSERT_EQ(isSetFinalBoundary(gopts["d"]), false);
  // Can't now set the value in the boundary
  ASSERT_THROW(set<int>(option["d"], 3), BoutException);
  ASSERT_THROW(setBoundary<int>(option["d"], 3), BoutException);
  // Can still set the value for the domain
  setNoBoundary<int>(option["d"], 1);
  ASSERT_EQ(getNoBoundary<int>(option["d"]), 1);

  ASSERT_EQ(isSetFinalBoundary(gopts["e"]), false);
  // Can still set the value, as isSetFinalBoundary doesn't mark variables it
  // doesn't have permission to read.
  set<int>(option["e"], 1);
  ASSERT_EQ(get<int>(option["e"]), 1);

  ASSERT_EQ(isSetFinalBoundary(gopts["f"]), false);
  // Can still set the value, as isSetFinalBoundary doesn't mark variables it
  // doesn't have permission to read.
  set<int>(option["f"], 1);
  ASSERT_EQ(get<int>(option["f"]), 1);
}
#endif

TEST_F(ComponentTest, Formatting) {
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

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

struct ConcreteComponentTests : public ComponentTest,
                                public testing::WithParamInterface<std::string> {
  static const Options base_options;
  static const Options required_params;
  static constexpr auto objname = "object_name";

  std::string typname;
  Options options;
  FakeSolver solver;

  ConcreteComponentTests()
      : ComponentTest(), typname(GetParam()), options(base_options.copy()) {
    Options::root()["mesh:paralleltransform:type"] = "identity";
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
    {"fieldline_geometry",
     {{"lambda_int", "0.01"},
      {"fieldline_radius", "2.0"},
      {"poloidal_magnetic_field", "0.5"},
      {"upstream_toroidal_magnetic_field", "1.0"}}},
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
