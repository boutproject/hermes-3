#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/component.hxx"
#include "../../include/external_apar.hxx"

#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/options.hxx>

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// Reuse the "standard" fixture for FakeMesh
using ExternalAparTest = FakeMeshFixture;

TEST_F(ExternalAparTest, CreateComponentNeedsApar) {
  Options options = {{"units", {{"Tesla", 2.0}, {"meters", 1.0e-3}}}};

  // No external_apar available
  EXPECT_THROW(ExternalApar component("test", options, nullptr), BoutException);
}

TEST_F(ExternalAparTest, CreateComponent) {
  Options options = {{"units", {{"Tesla", 2.0}, {"meters", 1.0e-3}}}};

  dynamic_cast<FakeMesh*>(bout::globals::mesh)
      ->setGridDataSource(new FakeGridDataSource{{{"external_apar", 1.0}}});
  // external_apar available
  ExternalApar component("test", options, nullptr);
}

TEST_F(ExternalAparTest, CreateComponentSetName) {
  Options options = {{"units", {{"Tesla", 2.0}, {"meters", 1.0e-3}}},
                     {"test", {{"apar_name", "AparExt"}}}};

  dynamic_cast<FakeMesh*>(bout::globals::mesh)
      ->setGridDataSource(new FakeGridDataSource{{{"AparExt", 1.0}}});
  // Read AparExt
  ExternalApar component("test", options, nullptr);
}

TEST_F(ExternalAparTest, TransformSet) {
  const BoutReal Bnorm = 2.0;
  const BoutReal Lnorm = 1.0e-3;
  Options options = {{"units", {{"Tesla", Bnorm}, {"meters", Lnorm}}}};

  const BoutReal external_apar = 3.1;
  dynamic_cast<FakeMesh*>(bout::globals::mesh)
      ->setGridDataSource(new FakeGridDataSource{{{"external_apar", external_apar}}});

  ExternalApar component("test", options, nullptr);

  Options state; // No Apar_flutter set
  component.transform(state);

  ASSERT_TRUE(state["fields"].isSet("Apar_flutter"));
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(state["fields"]["Apar_flutter"]),
                           external_apar / (Bnorm * Lnorm)));
}

TEST_F(ExternalAparTest, TransformAdd) {
  const BoutReal Bnorm = 2.0;
  const BoutReal Lnorm = 1.0e-3;
  Options options = {{"units", {{"Tesla", Bnorm}, {"meters", Lnorm}}}};

  const BoutReal external_apar = 3.1;
  dynamic_cast<FakeMesh*>(bout::globals::mesh)
      ->setGridDataSource(new FakeGridDataSource{{{"external_apar", external_apar}}});

  ExternalApar component("test", options, nullptr);

  Options state{{"fields", {{"Apar_flutter", 2.0}}}};
  component.transform(state);

  // Adds to existing Apar_flutter
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(state["fields"]["Apar_flutter"]),
                           2.0 + (external_apar / (Bnorm * Lnorm))));
}

TEST_F(ExternalAparTest, TransformAddScale) {
  const BoutReal Bnorm = 2.0;
  const BoutReal Lnorm = 1.0e-3;
  Options options = {{"units", {{"Tesla", Bnorm}, {"meters", Lnorm}}},
                     {"test", {{"scale", 4.0}}}};

  const BoutReal external_apar = 3.1;
  dynamic_cast<FakeMesh*>(bout::globals::mesh)
      ->setGridDataSource(new FakeGridDataSource{{{"external_apar", external_apar}}});

  ExternalApar component("test", options, nullptr);

  Options state{{"fields", {{"Apar_flutter", 2.0}}}};
  component.transform(state);

  // Adds to existing Apar_flutter
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(state["fields"]["Apar_flutter"]),
                           2.0 + (4.0 * external_apar / (Bnorm * Lnorm))));
}

TEST_F(ExternalAparTest, OutputScaleWithTransform) {
  const BoutReal Bnorm = 2.0;
  const BoutReal Lnorm = 1.0e-3;
  Options options = {{"units", {{"Tesla", Bnorm}, {"meters", Lnorm}}},
                     {"test", {{"scale", 4.0}}}};

  const BoutReal external_apar = 3.1;
  dynamic_cast<FakeMesh*>(bout::globals::mesh)
      ->setGridDataSource(new FakeGridDataSource{{{"external_apar", external_apar}}});

  ExternalApar component("test", options, nullptr);

  // Check that transform doesn't do anything funky to the output
  Options state{{"fields", {{"Apar_flutter", 2.0}}}};
  component.transform(state);

  Options output_state{{"Bnorm", Bnorm}, {"rho_s0", Lnorm}};
  component.outputVars(output_state);

  ASSERT_TRUE(output_state.isSet("external_apar"));

  // Did not add the existing Apar, only external Apar
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(output_state["external_apar"]),
                           (4.0 * external_apar / (Bnorm * Lnorm))));

  ASSERT_TRUE(output_state["external_apar"].attributes.contains("conversion"));
  ASSERT_DOUBLE_EQ(output_state["external_apar"].attributes["conversion"].as<BoutReal>(),
                   Bnorm * Lnorm);
}

TEST_F(ExternalAparTest, OutputScaleNoTransform) {
  const BoutReal Bnorm = 2.0;
  const BoutReal Lnorm = 1.0e-3;
  Options options = {{"units", {{"Tesla", Bnorm}, {"meters", Lnorm}}},
                     {"test", {{"scale", 4.0}}}};

  const BoutReal external_apar = 3.1;
  dynamic_cast<FakeMesh*>(bout::globals::mesh)
      ->setGridDataSource(new FakeGridDataSource{{{"external_apar", external_apar}}});

  ExternalApar component("test", options, nullptr);

  // Check that calling transform isn't necessary

  Options output_state{{"Bnorm", Bnorm}, {"rho_s0", Lnorm}};
  component.outputVars(output_state);

  ASSERT_TRUE(output_state.isSet("external_apar"));

  // Did not add the existing Apar, only external Apar
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(output_state["external_apar"]),
                           (4.0 * external_apar / (Bnorm * Lnorm))));

  ASSERT_TRUE(output_state["external_apar"].attributes.contains("conversion"));
  ASSERT_DOUBLE_EQ(output_state["external_apar"].attributes["conversion"].as<BoutReal>(),
                   Bnorm * Lnorm);
}
