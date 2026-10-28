
#include "gtest/gtest.h"

#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/coordinates.hxx>

#include "fake_mesh.hxx"
#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/diamagnetic_drift.hxx"

#include <memory>

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
  mesh->getCoordinates()->setBxy(1.0);
  static_cast<FakeMesh*>(mesh)->setGridDataSource(
      new FakeGridDataSource{{{"bxcvx", bx}, {"bxcvy", by}, {"bxcvz", bz}}});
}

BoutReal volumeIntegral(const Field3D& field) {
  const auto* coords = mesh->getCoordinates();
  const Field2D cell_volume = coords->dx() * coords->dy() * coords->dz() * coords->J();

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

  coords->setNon_uniform(true);
  coords->setD1_dx(0.2);
  coords->setD1_dy(0.2);
  coords->setD1_dz(0.0);
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

TEST_F(DiamagneticDriftTest, DefaultOptions) {
  setCurvature();
  Options options = diamagneticOptions();
  DiamagneticDrift component("test", options, nullptr);

  EXPECT_EQ(options["test"]["form"].as<DiamagneticForm>(), DiamagneticForm::divergence);
  EXPECT_EQ(options["test"]["radial_boundary"].as<DiamagneticBoundary>(),
            DiamagneticBoundary::closed_average);
}

TEST_F(DiamagneticDriftTest, AcceptsEveryFormAndRadialBoundary) {
  setCurvature();
  for (const auto form : {DiamagneticForm::gradient, DiamagneticForm::divergence}) {
    for (const auto boundary : {DiamagneticBoundary::open, DiamagneticBoundary::closed,
                                DiamagneticBoundary::closed_average}) {
      SCOPED_TRACE(toString(form) + ", " + toString(boundary));
      Options options = diamagneticOptions();
      options["test"]["form"] = toString(form);
      options["test"]["radial_boundary"] = toString(boundary);

      EXPECT_NO_THROW(DiamagneticDrift("test", options, nullptr));
      EXPECT_EQ(options["test"]["form"].as<DiamagneticForm>(), form);
      EXPECT_EQ(options["test"]["radial_boundary"].as<DiamagneticBoundary>(), boundary);
    }
  }
}

TEST_F(DiamagneticDriftTest, RejectsInvalidOptions) {
  setCurvature();
  for (const auto* key : {"form", "radial_boundary"}) {
    SCOPED_TRACE(key);
    Options options = diamagneticOptions();
    options["test"][key] = "invalid";
    EXPECT_THROW(DiamagneticDrift("test", options, nullptr), BoutException);
  }
}

TEST_F(DiamagneticDriftTest, RemovedDiamagFormExplainsMigration) {
  setCurvature();
  for (const BoutReal value : {0.0, 1.0, 0.5}) {
    SCOPED_TRACE(value);
    Options options = diamagneticOptions();
    options["test"]["diamag_form"] = value;
    try {
      DiamagneticDrift component("test", options, nullptr);
      FAIL() << "diamag_form must be rejected";
    } catch (const BoutException& error) {
      const std::string message = error.what();
      EXPECT_NE(message.find("diamag_form = 1 with form = divergence"),
                std::string::npos);
      EXPECT_NE(message.find("diamag_form = 0 with form = gradient"), std::string::npos);
      EXPECT_NE(message.find("Spatially varying blends were removed"), std::string::npos);
      EXPECT_NE(message.find("not conservative"), std::string::npos);
    }
  }
}

TEST_F(DiamagneticDriftTest, DeprecatedBoundaryFluxMapsToRadialBoundary) {
  setCurvature();
  for (const bool allow_flux : {false, true}) {
    SCOPED_TRACE(allow_flux);
    Options options = diamagneticOptions();
    options["test"]["bndry_flux"] = allow_flux;
    DiamagneticDrift component("test", options, nullptr);

    EXPECT_EQ(options["test"]["radial_boundary"].as<DiamagneticBoundary>(),
              allow_flux ? DiamagneticBoundary::open : DiamagneticBoundary::closed);
    EXPECT_EQ(options["test"]["form"].as<DiamagneticForm>(), DiamagneticForm::divergence);
  }
}

TEST_F(DiamagneticDriftTest, DeprecatedBoundaryFluxConflictsWithRadialBoundary) {
  setCurvature();
  for (const bool allow_flux : {false, true}) {
    for (const auto* boundary : {"open", "closed", "closed_average"}) {
      SCOPED_TRACE(boundary);
      SCOPED_TRACE(allow_flux);
      Options options = diamagneticOptions();
      options["test"]["bndry_flux"] = allow_flux;
      options["test"]["radial_boundary"] = boundary;
      EXPECT_THROW(DiamagneticDrift("test", options, nullptr), BoutException);
    }
  }
}

TEST_F(DiamagneticDriftTest, TransformSkipsSpeciesWithoutRequiredInputs) {
  setCurvature();
  Options options = diamagneticOptions();
  options["diamagnetic_drift"]["radial_boundary"] = "closed";

  DiamagneticDrift component("diamagnetic_drift", options, nullptr);

  Options state{
      {"species",
       {{"missing_temperature", {{"charge", 1.0}, {"density", 1.0}}},
        {"missing_charge", {{"temperature", 1.0}, {"density", 1.0}}},
        {"zero_charge", {{"charge", 0.0}, {"temperature", 1.0}, {"density", 1.0}}}}}};
  component.declareAllSpecies({"missing_temperature", "missing_charge", "zero_charge"});
  component.transform(state);

  ASSERT_FALSE(state["species"]["missing_temperature"].isSet("density_source"));
  ASSERT_FALSE(state["species"]["missing_charge"].isSet("density_source"));
  ASSERT_FALSE(state["species"]["zero_charge"].isSet("density_source"));
}

TEST_F(DiamagneticDriftTest, TransformWritesOnlyPresentOutputs) {
  setCurvature();
  Options options = diamagneticOptions();
  options["diamagnetic_drift"]["radial_boundary"] = "closed";
  options["diamagnetic_drift"]["form"] = "gradient";

  DiamagneticDrift component("diamagnetic_drift", options, nullptr);

  Options state{
      {"species",
       {{"d+",
         {{"charge", 1.0}, {"temperature", 1.0}, {"density", 1.0}, {"momentum", 1.0}}},
        {"e", {{"charge", -1.0}, {"temperature", 1.0}, {"pressure", 1.0}}}}}};
  component.declareAllSpecies({"d+", "e"});
  component.transform(state);

  ASSERT_TRUE(state["species"]["d+"].isSet("density_source"));
  ASSERT_TRUE(state["species"]["d+"].isSet("momentum_source"));
  ASSERT_FALSE(state["species"]["d+"].isSet("energy_source"));

  ASSERT_FALSE(state["species"]["e"].isSet("density_source"));
  ASSERT_TRUE(state["species"]["e"].isSet("energy_source"));
  ASSERT_FALSE(state["species"]["e"].isSet("momentum_source"));

  ASSERT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["d+"]["density_source"]), 0.0,
                           "RGN_NOBNDRY"));
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["d+"]["momentum_source"]), 0.0,
                           "RGN_NOBNDRY"));
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["e"]["energy_source"]), 0.0,
                           "RGN_NOBNDRY"));
}

TEST_F(WideDiamagneticDriftTest, TransformAveragesInnermostCoreRing) {
  setCurvature();

  auto* coords = mesh->getCoordinates();
  Coordinates::FieldMetric J{1.0};
  J(mesh->xstart, mesh->ystart) = 1.0;
  J(mesh->xstart, mesh->ystart + 1) = 2.0;
  J(mesh->xstart, mesh->yend) = 3.0;
  coords->setJ(J);

  Options options = diamagneticOptions();
  options["diamagnetic_drift"]["form"] = "gradient";
  DiamagneticDrift component("diamagnetic_drift", options, nullptr);

  Field3D quantity = makeField<Field3D>([](Ind3D& i) {
    const BoutReal slope =
        i.y() == mesh->ystart ? 1.0 : (i.y() == mesh->ystart + 1 ? 2.0 : 4.0);
    return 0.5 * i.x() * slope;
  });
  Options state;
  auto& species = state["species"]["d+"];
  species["charge"] = 1.0;
  species["temperature"] = 1.0;
  species["density"] = quantity;
  species["pressure"] = quantity;
  species["momentum"] = quantity;

  Options closed_state = state.copy();
  Options closed_options = diamagneticOptions();
  closed_options["diamagnetic_drift"]["form"] = "gradient";
  closed_options["diamagnetic_drift"]["radial_boundary"] = "closed";
  DiamagneticDrift closed_component("diamagnetic_drift", closed_options, nullptr);
  closed_component.declareAllSpecies({"d+"});
  closed_component.transform(closed_state);

  component.declareAllSpecies({"d+"});
  component.transform(state);

  for (const auto* source : {"density_source", "energy_source", "momentum_source"}) {
    SCOPED_TRACE(source);
    const auto averaged = get<Field3D>(species[source]);
    const auto closed = get<Field3D>(closed_state["species"]["d+"][source]);
    EXPECT_NE(closed(mesh->xstart, mesh->ystart, mesh->zstart),
              closed(mesh->xstart, mesh->yend, mesh->zstart));

    // Compare against the unaveraged closed source using the unequal cell volumes.
    BoutReal ring_integral = 0.0;
    BoutReal ring_volume = 0.0;
    const auto cell_volume = coords->cell_volume();
    for (int jy = mesh->ystart; jy <= mesh->yend; ++jy) {
      for (int jz = mesh->zstart; jz <= mesh->zend; ++jz) {
        ring_integral += cell_volume(mesh->xstart, jy, jz) * closed(mesh->xstart, jy, jz);
        ring_volume += cell_volume(mesh->xstart, jy, jz);
      }
    }
    const BoutReal expected_average = ring_integral / ring_volume;
    for (int jy = mesh->ystart; jy <= mesh->yend; ++jy) {
      for (int jz = mesh->zstart; jz <= mesh->zend; ++jz) {
        EXPECT_NEAR(averaged(mesh->xstart, jy, jz), expected_average, 1e-12);
        EXPECT_DOUBLE_EQ(averaged(mesh->xstart + 1, jy, jz),
                         closed(mesh->xstart + 1, jy, jz));
      }
    }
    EXPECT_NEAR(volumeIntegral(averaged), volumeIntegral(closed), 1e-12);
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
  options["diamagnetic_drift"]["radial_boundary"] = "closed";
  options["diamagnetic_drift"]["form"] = "divergence";

  DiamagneticDrift component("diamagnetic_drift", options, nullptr);

  Field3D quantity = makeField<Field3D>([](Ind3D& i) {
    const bool on_y_boundary = i.y() <= mesh->ystart || i.y() >= mesh->yend;
    return on_y_boundary ? 0.0 : 1.0 + 0.1 * i.z();
  });
  Field3D temperature{1.0};

  Options state;
  state["species"]["d+"]["charge"] = 1.0;
  state["species"]["d+"]["temperature"] = temperature;
  state["species"]["d+"]["density"] = quantity;
  component.declareAllSpecies({"d+"});
  component.transform(state);

  const auto source = get<Field3D>(state["species"]["d+"]["density_source"]);
  ASSERT_NEAR(volumeIntegral(source), 0.0, 1e-12);
}

TEST_F(WideDiamagneticDriftTest, GradientFormOpenIsConservativeWhenBoundaryValueUniform) {
  setCurvature(0.0, 1.0, 0.0);
  Options options = diamagneticOptions();
  options["diamagnetic_drift"]["radial_boundary"] = "open";
  options["diamagnetic_drift"]["form"] = "gradient";

  DiamagneticDrift component("diamagnetic_drift", options, nullptr);

  // With open radial boundaries, the gradient form's implied core flux vanishes
  // only when g = quantity * temperature / charge is uniform around the ring in
  // both the x guard cell and the first interior cell.
  Field3D quantity = makeField<Field3D>([](Ind3D& i) {
    const bool on_boundary = i.x() <= mesh->xstart || i.x() >= mesh->xend
                             || i.y() <= mesh->ystart || i.y() >= mesh->yend;
    return on_boundary ? 1.0 : 1.5 + 0.1 * i.y() + 0.05 * i.z();
  });
  Field3D temperature{1.0};

  Options state;
  state["species"]["d+"]["charge"] = 1.0;
  state["species"]["d+"]["temperature"] = temperature;
  state["species"]["d+"]["density"] = quantity;
  component.declareAllSpecies({"d+"});
  component.transform(state);

  const auto source = get<Field3D>(state["species"]["d+"]["density_source"]);
  ASSERT_NEAR(volumeIntegral(source), 0.0, 1e-12);
}
