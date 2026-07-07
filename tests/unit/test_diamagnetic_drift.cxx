
#include "gtest/gtest.h"

#include <bout/boutexception.hxx>

#include "fake_mesh.hxx"
#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/diamagnetic_drift.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using DiamagneticDriftTest = FakeMeshFixture;
using WideDiamagneticDriftTest = FakeMeshFixture_tmpl<5, 5, 7>;

namespace {

Options diamagneticOptions() {
  Options options;
  options["units"]["Tesla"] = 1.0;
  options["units"]["meters"] = 1.0;
  return options;
}

void setCurvature(BoutReal bx = 1.0, BoutReal by = 0.0, BoutReal bz = 0.0) {
  mesh->getCoordinates()->Bxy = 1.0;
  static_cast<FakeMesh*>(mesh)->setGridDataSource(
      new FakeGridDataSource{{{"bxcvx", bx}, {"bxcvy", by}, {"bxcvz", bz}}});
}

BoutReal volumeIntegral(const Field3D& field) {
  const auto* coords = mesh->getCoordinates();
  const Field2D cell_volume = coords->dx * coords->dy * coords->dz * coords->J;

  BoutReal result = 0.0;
  for (int jx = mesh->xstart; jx <= mesh->xend; ++jx) {
    for (int jy = mesh->ystart; jy <= mesh->yend; ++jy) {
      for (int jz = mesh->zstart; jz <= mesh->zend; ++jz) {
        result += cell_volume(jx, jy) * field(jx, jy, jz);
      }
    }
  }
  return result;
}

class ScopedGlobalMesh {
public:
  explicit ScopedGlobalMesh(Mesh* new_mesh) : previous_mesh(mesh) { mesh = new_mesh; }
  ~ScopedGlobalMesh() { mesh = previous_mesh; }

private:
  Mesh* previous_mesh;
};

std::shared_ptr<Coordinates> makeUnitCoordinates(Mesh* target_mesh) {
  auto coords = std::make_shared<Coordinates>(
      target_mesh, Field2D{1.0, target_mesh}, Field2D{1.0, target_mesh},
      Field2D{1.0, target_mesh}, Field2D{1.0, target_mesh}, Field2D{1.0, target_mesh},
      Field2D{1.0, target_mesh}, Field2D{1.0, target_mesh}, Field2D{1.0, target_mesh},
      Field2D{0.0, target_mesh}, Field2D{0.0, target_mesh}, Field2D{0.0, target_mesh},
      Field2D{1.0, target_mesh}, Field2D{1.0, target_mesh}, Field2D{1.0, target_mesh},
      Field2D{0.0, target_mesh}, Field2D{0.0, target_mesh}, Field2D{0.0, target_mesh},
      Field2D{0.0, target_mesh}, Field2D{0.0, target_mesh});

  coords->G1 = coords->G2 = coords->G3 = 0.1;
  coords->non_uniform = true;
  coords->d1_dx = coords->d1_dy = 0.2;
  coords->d1_dz = 0.0;
#if BOUT_USE_METRIC_3D
  coords->Bxy.splitParallelSlices();
  coords->Bxy.yup() = coords->Bxy.ydown() = coords->Bxy;
#endif
  coords->setParallelTransform(
      bout::utils::make_unique<ParallelTransformIdentity>(*target_mesh));

  return coords;
}

} // namespace

TEST_F(DiamagneticDriftTest, CreateComponent) {
  Options options = diamagneticOptions();
  setCurvature();
  DiamagneticDrift component("test", options, nullptr);
}

TEST_F(DiamagneticDriftTest, AddDiamagneticSourcesSkipsSpeciesWithoutRequiredInputs) {
  setCurvature();
  Options options = diamagneticOptions();
  options["diamagnetic_drift"]["average_core"] = false;

  DiamagneticDrift component("diamagnetic_drift", options, nullptr);

  Options state{
      {"species",
       {{"missing_temperature", {{"charge", 1.0}, {"density", 1.0}}},
        {"missing_charge", {{"temperature", 1.0}, {"density", 1.0}}},
        {"zero_charge", {{"charge", 0.0}, {"temperature", 1.0}, {"density", 1.0}}}}}};
  Permissions permissions{{readOnly("species"),
                           readWrite("species:missing_temperature:density_source"),
                           readWrite("species:missing_charge:density_source"),
                           readWrite("species:zero_charge:density_source")}};
  GuardedOptions guarded_state{&state, &permissions};

  auto missing_temperature = guarded_state["species"]["missing_temperature"];
  component.addDiamagneticSources(missing_temperature);
  auto missing_charge = guarded_state["species"]["missing_charge"];
  component.addDiamagneticSources(missing_charge);
  auto zero_charge = guarded_state["species"]["zero_charge"];
  component.addDiamagneticSources(zero_charge);

  ASSERT_FALSE(state["species"]["missing_temperature"].isSet("density_source"));
  ASSERT_FALSE(state["species"]["missing_charge"].isSet("density_source"));
  ASSERT_FALSE(state["species"]["zero_charge"].isSet("density_source"));
}

TEST_F(DiamagneticDriftTest, AddDiamagneticSourcesWritesOnlyPresentOutputs) {
  setCurvature();
  Options options = diamagneticOptions();
  options["diamagnetic_drift"]["average_core"] = false;
  options["diamagnetic_drift"]["divergence_form"] = false;

  DiamagneticDrift component("diamagnetic_drift", options, nullptr);

  Options state{
      {"species",
       {{"d+",
         {{"charge", 1.0}, {"temperature", 1.0}, {"density", 1.0}, {"momentum", 1.0}}},
        {"e", {{"charge", -1.0}, {"temperature", 1.0}, {"pressure", 1.0}}}}}};
  Permissions permissions{{readOnly("species"), readWrite("species:d+:density_source"),
                           readWrite("species:d+:momentum_source"),
                           readWrite("species:e:energy_source")}};
  GuardedOptions guarded_state{&state, &permissions};

  auto d_ion = guarded_state["species"]["d+"];
  component.addDiamagneticSources(d_ion);
  auto electron = guarded_state["species"]["e"];
  component.addDiamagneticSources(electron);

  ASSERT_TRUE(state["species"]["d+"].isSet("density_source"));
  ASSERT_TRUE(state["species"]["d+"].isSet("momentum_source"));
  ASSERT_FALSE(state["species"]["d+"].isSet("energy_source"));

  ASSERT_FALSE(state["species"]["e"].isSet("density_source"));
  ASSERT_TRUE(state["species"]["e"].isSet("energy_source"));
  ASSERT_FALSE(state["species"]["e"].isSet("momentum_source"));
}

TEST_F(DiamagneticDriftTest, CoreAverageRequiresAverageCore) {
  setCurvature();
  Options options = diamagneticOptions();
  options["diamagnetic_drift"]["average_core"] = false;

  DiamagneticDrift component("diamagnetic_drift", options, nullptr);

  Field3D field{0.0};
  EXPECT_THROW(component.coreAverage(field), BoutException);
}

TEST_F(DiamagneticDriftTest, CoreAverageAveragesInnermostCoreRing) {
  setCurvature();

  auto* coords = mesh->getCoordinates();
  coords->J = 1.0;
  coords->J(mesh->xstart, mesh->ystart) = 1.0;
  coords->J(mesh->xstart, mesh->ystart + 1) = 2.0;
  coords->J(mesh->xstart, mesh->yend) = 3.0;

  Options options = diamagneticOptions();
  DiamagneticDrift component("diamagnetic_drift", options, nullptr);

  Field3D field{0.0};
  field = 0.0;
  for (int jy = mesh->ystart; jy <= mesh->yend; ++jy) {
    const BoutReal ring_value =
        jy == mesh->ystart ? 1.0 : (jy == mesh->ystart + 1 ? 2.0 : 4.0);
    for (int jz = mesh->zstart; jz <= mesh->zend; ++jz) {
      field(mesh->xstart, jy, jz) = ring_value;
      field(mesh->xstart + 1, jy, jz) = 9.0;
    }
  }

  component.coreAverage(field);

  const BoutReal expected_average = (1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 4.0) / 6.0;
  for (int jy = mesh->ystart; jy <= mesh->yend; ++jy) {
    for (int jz = mesh->zstart; jz <= mesh->zend; ++jz) {
      ASSERT_DOUBLE_EQ(field(mesh->xstart, jy, jz), expected_average);
      ASSERT_DOUBLE_EQ(field(mesh->xstart + 1, jy, jz), 9.0);
    }
  }
}

TEST_F(DiamagneticDriftTest, DivergenceFormIsConservativeWithoutBoundaryFluxes) {
  FakeMesh div_mesh(6, 5, 7, *bout::globals::mpi);
  div_mesh.xstart = 2;
  div_mesh.xend = div_mesh.LocalNx - 3;
  div_mesh.createDefaultRegions();
  div_mesh.setCoordinates(nullptr);
  auto div_coords = makeUnitCoordinates(&div_mesh);
  div_mesh.setCoordinates(div_coords);
  div_mesh.setGridDataSource(new FakeGridDataSource());
  div_mesh.createBoundaryRegions();

  ScopedGlobalMesh use_div_mesh(&div_mesh);

  setCurvature(0.0, 1.0, 0.0);
  Options options = diamagneticOptions();
  options["diamagnetic_drift"]["average_core"] = false;
  options["diamagnetic_drift"]["divergence_form"] = true;
  options["diamagnetic_drift"]["bndry_flux"] = false;

  DiamagneticDrift component("diamagnetic_drift", options, nullptr);

  Field3D quantity = makeField<Field3D>([](Ind3D& i) {
    const bool on_y_boundary = i.y() <= mesh->ystart || i.y() >= mesh->yend;
    return on_y_boundary ? 0.0 : 1.0 + 0.1 * i.z();
  });
  Field3D temperature{1.0};

  Field3D sink = component.calculateDivergenceForm(quantity, temperature, 1.0);

  ASSERT_NEAR(volumeIntegral(sink), 0.0, 1e-12);
}

TEST_F(WideDiamagneticDriftTest, GradientFormIsConservativeWithConstantBoundaryProduct) {
  setCurvature(0.0, 1.0, 0.0);
  Options options = diamagneticOptions();
  options["diamagnetic_drift"]["average_core"] = false;
  options["diamagnetic_drift"]["divergence_form"] = false;

  DiamagneticDrift component("diamagnetic_drift", options, nullptr);

  Field3D quantity = makeField<Field3D>([](Ind3D& i) {
    const bool on_boundary = i.x() <= mesh->xstart || i.x() >= mesh->xend
                             || i.y() <= mesh->ystart || i.y() >= mesh->yend;
    return on_boundary ? 1.0 : 1.5 + 0.1 * i.y() + 0.05 * i.z();
  });
  Field3D temperature{1.0};

  Field3D sink = component.calculateGradientForm(quantity, temperature, 1.0);

  ASSERT_NEAR(volumeIntegral(sink), 0.0, 1e-12);
}
