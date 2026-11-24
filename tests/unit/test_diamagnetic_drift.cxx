#include "gtest/gtest.h"

#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/build_defines.hxx>
#include <bout/constants.hxx>
#include <bout/coordinates.hxx>
#include <bout/field.hxx>
#include <bout/field3d.hxx>
#include <bout/globals.hxx>
#include <bout/mesh.hxx>
#include <bout/paralleltransform.hxx>
#include <bout/region.hxx>
#include <bout/sys/range.hxx>

#include "fake_mesh.hxx"
#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/component.hxx"
#include "../../include/diamagnetic_drift.hxx"

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// The unit tests use the global mesh
using namespace bout::globals;

// Make target queries consistent with periodicY(): core rings have no targets,
// while open field lines have a target at each end. FakeMesh defaults do not
// distinguish these cases in their target iterators.
class DiamagneticFakeMesh : public FakeMesh {
public:
  using FakeMesh::FakeMesh; // Inherit the grid-size and guard-cell constructors.
  // Keep the no-argument overloads visible when overriding firstY(int)/lastY(int).
  using FakeMesh::firstY;
  using FakeMesh::lastY;

  RangeIterator iterateBndryLowerY() const override {
    const int first_target = std::max(xstart, ix_separatrix);
    return first_target <= xend ? RangeIterator(first_target, xend) : RangeIterator{};
  }
  RangeIterator iterateBndryUpperY() const override { return iterateBndryLowerY(); }
  bool hasBndryLowerY() const override { return ix_separatrix <= xend; }
  bool hasBndryUpperY() const override { return hasBndryLowerY(); }
  bool firstY(int jx) const override { return !periodicY(jx); }
  bool lastY(int jx) const override { return !periodicY(jx); }
};

namespace {

Options diamagneticOptions() {
  Options options;
  options["units"]["Tesla"] = 1.0;
  options["units"]["meters"] = 1.0;
  return options;
}

Options diamagneticState(const Field3D& quantity, BoutReal temperature = 1.0,
                         BoutReal charge = 1.0) {
  Options state;
  auto& species = state["species"]["d+"];
  species["charge"] = charge;
  species["temperature"] = temperature;
  species["density"] = quantity;
  species["pressure"] = quantity;
  species["momentum"] = quantity;
  return state;
}

void setCurvature(BoutReal bx = 1.0, BoutReal by = 0.0, BoutReal bz = 0.0) {
  mesh->getCoordinates()->setBxy(1.0);
  static_cast<FakeMesh*>(mesh)->setGridDataSource(
      new FakeGridDataSource{{{"bxcvx", bx}, {"bxcvy", by}, {"bxcvz", bz}}});
}

BoutReal volumeIntegral(const Field3D& field) {
  const auto* coords = mesh->getCoordinates();
  const Coordinates::FieldMetric cell_volume = coords->cell_volume();

  BoutReal result = 0.0;
  BOUT_FOR(i, mesh->getRegion3D("RGN_NOBNDRY")) { result += cell_volume[i] * field[i]; }
  return result;
}

class DivergenceDiamagneticMesh : public DiamagneticFakeMesh {
public:
  DivergenceDiamagneticMesh(int nx, int ny, int nz, MpiWrapper& mpi, int mxg)
      : DiamagneticFakeMesh(nx, ny, nz, mpi, mxg) {
    ix_separatrix = 0;
    xend = LocalNx - mxg - 1;
  }
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
  auto Bxy = coords->Bxy();
  Bxy.splitParallelSlices();
  Bxy.yup() = Bxy.ydown() = coords->Bxy();
  coords->setBxy(Bxy);
#endif
  coords->setParallelTransform(
      bout::utils::make_unique<ParallelTransformIdentity>(*target_mesh));

  return coords;
}

} // namespace

// Aligned gradient/divergence tests start with B = 1 and radial curvature.
// Tests needing another direction or a varying B override this default.
template <class Base>
class CurvatureDiamagneticFixture : public Base {
protected:
  void SetUp() override {
    Base::SetUp();
    setCurvature();
  }
};

// One physical x cell is sufficient for option and species/source-selection tests.
using DiamagneticDriftTest = CurvatureDiamagneticFixture<
    FakeMeshFixture_tmpl<3, 5, 7, false, 1, DiamagneticFakeMesh>>;
// Three physical x cells leave a nonuniform interior between uniform boundary
// cells for the open-gradient conservation test.
using WideDiamagneticDriftTest = CurvatureDiamagneticFixture<
    FakeMeshFixture_tmpl<5, 5, 7, false, 1, DiamagneticFakeMesh>>;
// Four physical x cells allow comparison of the averaged and adjacent rings.
// Two x guards are required by FV::Div_f_v's reconstruction stencil.
using RadialDiamagneticDriftTest = CurvatureDiamagneticFixture<
    FakeMeshFixture_tmpl<8, 5, 7, false, 2, DiamagneticFakeMesh>>;
// The FCI mock cannot convert to/from aligned coordinates, so gradient/divergence
// construction must fail before reading curvature; a minimal mesh suffices.
using UnsupportedDiamagneticDriftTest =
    FakeMeshFixture_tmpl<3, 5, 7, true, 1, DiamagneticFakeMesh>;
// A wider FCI mesh supports nonuniform bracket profiles and core-ring averaging.
// Deliberately omit curvature: the bracket uses B and the metrics instead.
using FciDiamagneticDriftTest =
    FakeMeshFixture_tmpl<8, 5, 7, true, 2, DiamagneticFakeMesh>;

// Open field lines with two targets, unit metrics and poloidal curvature test
// divergence conservation without flux through any boundary.
class DivergenceDiamagneticDriftTest
    : public FakeMeshFixture_tmpl<6, 5, 7, false, 2, DivergenceDiamagneticMesh> {
protected:
  void SetUp() override {
    FakeMeshFixture_tmpl::SetUp();
    setCurvature(0.0, 1.0, 0.0);
  }

  std::shared_ptr<Coordinates> getCoordinates() override {
    return makeUnitCoordinates(&getMesh());
  }
};

// Unspecified options select divergence form and closed_average radial boundaries.
TEST_F(DiamagneticDriftTest, DefaultOptions) {
  Options options = diamagneticOptions();
  DiamagneticDrift component("test", options, nullptr);

  EXPECT_EQ(options["test"]["form"].as<DiamagneticForm>(), DiamagneticForm::divergence);
  EXPECT_EQ(options["test"]["radial_boundary"].as<DiamagneticBoundary>(),
            DiamagneticBoundary::closed_average);
}

// Every valid form/boundary combination constructs successfully and retains its options.
TEST_F(DiamagneticDriftTest, AcceptsEveryFormAndRadialBoundary) {
  for (const auto form : {DiamagneticForm::gradient, DiamagneticForm::divergence,
                          DiamagneticForm::bracket}) {
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

// Unknown form or radial_boundary names are rejected during construction.
TEST_F(DiamagneticDriftTest, RejectsInvalidOptions) {
  for (const auto* key : {"form", "radial_boundary"}) {
    SCOPED_TRACE(key);
    Options options = diamagneticOptions();
    options["test"][key] = "invalid";
    EXPECT_THROW(DiamagneticDrift("test", options, nullptr), BoutException);
  }
}

// The removed diamag_form option reports its replacements and why blends were removed.
TEST_F(DiamagneticDriftTest, RemovedDiamagFormExplainsMigration) {
  for (const BoutReal value : {0.0, 1.0, 0.5}) {
    SCOPED_TRACE(value);
    Options options = diamagneticOptions();
    options["test"]["diamag_form"] = value;
    try {
      const DiamagneticDrift component("test", options, nullptr);
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

// Removed bndry_flux reports its replacements rather than acting as an alias.
TEST_F(DiamagneticDriftTest, RemovedBoundaryFluxExplainsMigration) {
  for (const bool allow_flux : {false, true}) {
    SCOPED_TRACE(allow_flux);
    Options options = diamagneticOptions();
    options["test"]["bndry_flux"] = allow_flux;
    try {
      const DiamagneticDrift component("test", options, nullptr);
      FAIL() << "bndry_flux must be rejected";
    } catch (const BoutException& error) {
      const std::string message = error.what();
      EXPECT_NE(message.find("bndry_flux has been removed"), std::string::npos);
      EXPECT_NE(message.find("radial_boundary = open"), std::string::npos);
      EXPECT_NE(message.find("radial_boundary = closed"), std::string::npos);
    }
  }
}

// FCI meshes default to the bracket form without requiring a curvature vector.
TEST_F(FciDiamagneticDriftTest, DefaultOptionsSelectBracket) {
  ASSERT_TRUE(mesh->isFci());
  Options options = diamagneticOptions();
  const DiamagneticDrift component("test", options, nullptr);
  EXPECT_EQ(options["test"]["form"].as<DiamagneticForm>(), DiamagneticForm::bracket);
  EXPECT_EQ(options["test"]["radial_boundary"].as<DiamagneticBoundary>(),
            DiamagneticBoundary::closed_average);
}

// For log(B) linear in x, the bracket differentiates a periodic z profile.
// Check its sign, metric factor, charge scaling, energy factor and source accumulation.
TEST_F(FciDiamagneticDriftTest, BracketSourcesUseFciMetricsAndCharge) {
  const BoutReal wave_number = TWOPI / mesh->LocalNz;
  const BoutReal slope = 0.2;
  auto* coords = mesh->getCoordinates();
  coords->setBxy(makeField<Coordinates::FieldMetric>(
      [slope](auto& i) { return std::exp(slope * i.x()); }));
  coords->setMetricTensorJB(
      ContravariantMetricTensor(1.0, 1.0 / 9.0, 1.0, 0.0, 0.0, 0.0, mesh),
      CovariantMetricTensor(1.0, 9.0, 1.0, 0.0, 0.0, 0.0, mesh),
      Coordinates::FieldMetric{2.0}, coords->Bxy());
  coords->setDx(0.5);
  coords->setDz(2.0);
  const auto quantity = makeField<Field3D>(
      [wave_number](Ind3D& i) { return 2.0 + std::sin(wave_number * i.z()); });

  for (const BoutReal charge : {1.0, -2.0}) {
    SCOPED_TRACE(charge);
    Options options = diamagneticOptions();
    options["test"]["radial_boundary"] = "open";
    DiamagneticDrift component("test", options, nullptr);
    component.declareAllSpecies({"d+"});
    auto state = diamagneticState(quantity, 4.0, charge);
    for (const auto* output : {"density_source", "energy_source", "momentum_source"}) {
      state["species"]["d+"][output] = 7.0;
    }
    component.transform(state);

    for (const auto* output : {"density_source", "energy_source", "momentum_source"}) {
      SCOPED_TRACE(output);
      const auto source = get<Field3D>(state["species"]["d+"][output]);
      const BoutReal factor = std::string(output) == "energy_source" ? 2.5 : 1.0;
      BOUT_FOR_SERIAL(i, source.getRegion("RGN_NOBNDRY")) {
        // [log(B), g] = -d_x log(B) * d_z g, with C2's sin(k) factor.
        const BoutReal expected =
            7.0
            + 2.0 * factor * slope * (4.0 / charge) * std::sin(wave_number)
                  * std::cos(wave_number * i.z()) * 3.0 / (2.0 * std::exp(slope * i.x()));
        EXPECT_NEAR(source[i], expected, 1e-12);
      }
    }
  }
}

// Explicit bracket selection works with an aligned transform and omits the FCI factor.
TEST_F(RadialDiamagneticDriftTest, ExplicitBracketUsesAlignedMeshNormalisation) {
  const BoutReal wave_number = TWOPI / mesh->LocalNz;
  auto* coords = mesh->getCoordinates();
  coords->setBxy(
      makeField<Coordinates::FieldMetric>([](auto& i) { return std::exp(0.2 * i.x()); }));
  coords->setMetricTensorJB(
      ContravariantMetricTensor(1.0, 1.0 / 9.0, 1.0, 0.0, 0.0, 0.0, mesh),
      CovariantMetricTensor(1.0, 9.0, 1.0, 0.0, 0.0, 0.0, mesh),
      Coordinates::FieldMetric{2.0}, coords->Bxy());
  Options options = diamagneticOptions();
  options["test"]["form"] = "bracket";
  options["test"]["radial_boundary"] = "open";
  DiamagneticDrift component("test", options, nullptr);
  component.declareAllSpecies({"d+"});
  auto state = diamagneticState(makeField<Field3D>(
      [wave_number](Ind3D& i) { return std::sin(wave_number * i.z()); }));
  component.transform(state);
  const auto source = get<Field3D>(state["species"]["d+"]["density_source"]);
  BOUT_FOR_SERIAL(i, source.getRegion("RGN_NOBNDRY")) {
    EXPECT_NEAR(source[i], 0.4 * std::sin(wave_number) * std::cos(wave_number * i.z()),
                1e-12);
  }
}

// Core averaging remains available for the default FCI bracket form and preserves sums.
TEST_F(FciDiamagneticDriftTest, BracketCoreAveragingPreservesRingIntegral) {
  const BoutReal wave_number = TWOPI / mesh->LocalNz;
  auto* coords = mesh->getCoordinates();
  coords->setBxy(
      makeField<Coordinates::FieldMetric>([](auto& i) { return std::exp(0.2 * i.x()); }));
  coords->setJ(
      makeField<Coordinates::FieldMetric>([](auto& i) { return 1.0 + 0.2 * i.y(); }));
  const auto quantity = makeField<Field3D>([wave_number](Ind3D& i) {
    return (1.0 + 0.1 * i.y()) * std::sin(wave_number * i.z());
  });
  Options options = diamagneticOptions();
  Options closed_options = diamagneticOptions();
  closed_options["test"]["radial_boundary"] = "closed";
  DiamagneticDrift averaged("test", options, nullptr);
  DiamagneticDrift closed("test", closed_options, nullptr);
  averaged.declareAllSpecies({"d+"});
  closed.declareAllSpecies({"d+"});
  auto averaged_state = diamagneticState(quantity);
  auto closed_state = diamagneticState(quantity);
  averaged.transform(averaged_state);
  closed.transform(closed_state);

  for (const auto* output : {"density_source", "energy_source", "momentum_source"}) {
    SCOPED_TRACE(output);
    const auto source = get<Field3D>(averaged_state["species"]["d+"][output]);
    const auto unaveraged = get<Field3D>(closed_state["species"]["d+"][output]);
    BoutReal ring_sum = 0.0;
    BoutReal ring_volume = 0.0;
    for (int jy = mesh->ystart; jy <= mesh->yend; ++jy) {
      for (int jz = mesh->zstart; jz <= mesh->zend; ++jz) {
        const BoutReal volume = coords->cell_volume()(mesh->xstart, jy, jz);
        ring_sum += volume * unaveraged(mesh->xstart, jy, jz);
        ring_volume += volume;
      }
    }
    EXPECT_GT(std::abs(unaveraged(mesh->xstart, mesh->ystart, mesh->zstart)), 1e-6);
    BOUT_FOR_SERIAL(i, source.getRegion("RGN_NOBNDRY")) {
      EXPECT_NEAR(source[i],
                  i.x() == mesh->xstart ? ring_sum / ring_volume : unaveraged[i], 1e-12);
    }
    EXPECT_NEAR(volumeIntegral(source), volumeIntegral(unaveraged), 1e-12);
  }
}

// Species without charge or temperature, and neutral species, acquire no density source.
TEST_F(DiamagneticDriftTest, TransformSkipsSpeciesWithoutRequiredInputs) {
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

// Only present density, pressure and momentum inputs produce their corresponding sources;
// these uniform profiles give zero sources, including for negative charge.
TEST_F(DiamagneticDriftTest, TransformWritesOnlyPresentOutputs) {
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

// Both forms replace core-ring sources with their volume-weighted means, preserving
// the volume integrals and leaving the adjacent ring unchanged.
TEST_F(RadialDiamagneticDriftTest, TransformAveragesInnermostCoreRing) {
  auto* coords = mesh->getCoordinates();
  Coordinates::FieldMetric J{1.0};
  for (int jz = mesh->zstart; jz <= mesh->zend; ++jz) {
    J(mesh->xstart, mesh->ystart + 1, jz) = 2.0;
    J(mesh->xstart, mesh->yend, jz) = 3.0;
  }
  coords->setJ(J);

  for (const auto* form : {"gradient", "divergence"}) {
    SCOPED_TRACE(form);
    Options options = diamagneticOptions();
    options["diamagnetic_drift"]["form"] = form;
    DiamagneticDrift component("diamagnetic_drift", options, nullptr);

    const Field3D quantity = makeField<Field3D>([](Ind3D& i) {
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
    closed_options["diamagnetic_drift"]["form"] = form;
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
          ring_integral +=
              cell_volume(mesh->xstart, jy, jz) * closed(mesh->xstart, jy, jz);
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
}

#if BOUT_USE_METRIC_3D
// Core means include z-dependent cell volumes, preserving ring integrals rather than
// using arithmetic means; the adjacent ring is unchanged and closed sources conserve.
TEST_F(RadialDiamagneticDriftTest, CoreAveragePreservesRingIntegralWithZVaryingJacobian) {
  const BoutReal wave_number = TWOPI / mesh->LocalNz;
  const auto J = makeField<Field3D>([wave_number](Ind3D& i) {
    return 1.5 + 0.2 * i.y() + 0.3 * std::cos(wave_number * i.z());
  });
  auto* coords = mesh->getCoordinates();
  coords->setJ(J);
  // Bxy = J makes J*C.x constant, so the discrete curvature is divergence-free.
  coords->setBxy(J);
  const auto quantity = makeField<Field3D>([wave_number](Ind3D& i) {
    return 1.0 + 0.2 * i.x() * i.x() + 0.1 * i.x() * i.y()
           + 0.5 * i.x() * std::sin(wave_number * i.z());
  });
  const auto& cell_volume = coords->cell_volume();

  for (const auto* form : {"gradient", "divergence"}) {
    SCOPED_TRACE(form);
    Options options = diamagneticOptions();
    options["test"]["form"] = form;
    DiamagneticDrift component("test", options, nullptr);
    Options closed_options = diamagneticOptions();
    closed_options["test"]["form"] = form;
    closed_options["test"]["radial_boundary"] = "closed";
    DiamagneticDrift closed_component("test", closed_options, nullptr);
    auto state = diamagneticState(quantity, 1.7, -2.0);
    auto closed_state = state.copy();
    component.declareAllSpecies({"d+"});
    closed_component.declareAllSpecies({"d+"});
    component.transform(state);
    closed_component.transform(closed_state);

    for (const auto* name : {"density_source", "energy_source", "momentum_source"}) {
      SCOPED_TRACE(name);
      const auto averaged = get<Field3D>(state["species"]["d+"][name]);
      const auto closed = get<Field3D>(closed_state["species"]["d+"][name]);
      BoutReal ring_volume = 0.0;
      BoutReal closed_integral = 0.0;
      BoutReal averaged_integral = 0.0;
      BoutReal unweighted_sum = 0.0;
      int ring_cells = 0;
      for (int jy = mesh->ystart; jy <= mesh->yend; ++jy) {
        for (int jz = mesh->zstart; jz <= mesh->zend; ++jz) {
          const BoutReal volume = cell_volume(mesh->xstart, jy, jz);
          ring_volume += volume;
          closed_integral += volume * closed(mesh->xstart, jy, jz);
          averaged_integral += volume * averaged(mesh->xstart, jy, jz);
          unweighted_sum += closed(mesh->xstart, jy, jz);
          ++ring_cells;
        }
      }
      const BoutReal expected = closed_integral / ring_volume;
      // The profiles distinguish volume weighting from an arithmetic mean.
      EXPECT_GT(std::abs(expected - unweighted_sum / ring_cells), 1e-4);
      EXPECT_NEAR(averaged_integral, closed_integral, 1e-12);
      for (int jy = mesh->ystart; jy <= mesh->yend; ++jy) {
        for (int jz = mesh->zstart; jz <= mesh->zend; ++jz) {
          EXPECT_NEAR(averaged(mesh->xstart, jy, jz), expected, 1e-12);
          EXPECT_DOUBLE_EQ(averaged(mesh->xstart + 1, jy, jz),
                           closed(mesh->xstart + 1, jy, jz));
        }
      }
      EXPECT_NEAR(volumeIntegral(averaged), 0.0, 1e-12);
      EXPECT_NEAR(volumeIntegral(closed), 0.0, 1e-12);
    }
  }
}

// The shifted curvature correction retains the full 3D sinty profile in every source.
TEST_F(RadialDiamagneticDriftTest, ShiftedCurvatureReadsThreeDimensionalSinty) {
  const BoutReal wave_number = TWOPI / mesh->LocalNz;
  const auto sinty = makeField<Field3D>([wave_number](Ind3D& i) {
    return 0.2 + 0.05 * i.y() + 0.1 * std::sin(wave_number * i.z());
  });
  static_cast<FakeMesh*>(mesh)->setGridDataSource(new FakeGridDataSource{
      {{"bxcvx", 1.0}, {"bxcvy", 0.0}, {"bxcvz", 0.0}, {"sinty", sinty}}});
  // Select the curvature correction while retaining the identity transform
  // to isolate its dependence on the 3D sinty field.
  Options::root()["mesh"]["paralleltransform"]["type"] = "shifted";
  Options options = diamagneticOptions();
  options["test"]["form"] = "gradient";
  options["test"]["radial_boundary"] = "open";
  DiamagneticDrift component("test", options, nullptr);
  const auto quantity = makeField<Field3D>(
      [wave_number](Ind3D& i) { return std::sin(wave_number * i.z()); });
  auto state = diamagneticState(quantity);
  component.declareAllSpecies({"d+"});
  component.transform(state);

  for (const auto* name : {"density_source", "energy_source", "momentum_source"}) {
    SCOPED_TRACE(name);
    const auto source = get<Field3D>(state["species"]["d+"][name]);
    const BoutReal factor = std::string(name) == "energy_source" ? 2.5 : 1.0;
    BOUT_FOR_SERIAL(i, source.getRegion("RGN_NOBNDRY")) {
      const BoutReal expected = -2.0 * factor * sinty[i] * std::sin(wave_number)
                                * std::cos(wave_number * i.z());
      EXPECT_NEAR(source[i], expected, 1e-12);
    }
  }
}
#endif

TEST_F(DivergenceDiamagneticDriftTest,
       DivergenceFormIsConservativeWithoutBoundaryFluxes) {
  Options options = diamagneticOptions();
  options["diamagnetic_drift"]["radial_boundary"] = "closed";
  options["diamagnetic_drift"]["form"] = "divergence";

  DiamagneticDrift component("diamagnetic_drift", options, nullptr);

  const Field3D quantity = makeField<Field3D>([](Ind3D& i) {
    const bool on_y_boundary = i.y() <= mesh->ystart || i.y() >= mesh->yend;
    return on_y_boundary ? 0.0 : 1.0 + 0.1 * i.z();
  });
  const Field3D temperature{1.0};

  Options state;
  state["species"]["d+"]["charge"] = 1.0;
  state["species"]["d+"]["temperature"] = temperature;
  state["species"]["d+"]["density"] = quantity;
  component.declareAllSpecies({"d+"});
  component.transform(state);

  const auto source = get<Field3D>(state["species"]["d+"]["density_source"]);
  ASSERT_NEAR(volumeIntegral(source), 0.0, 1e-12);
}

// With uniform boundary values and divergence-free curvature, open gradient transport
// has zero volume-integrated density source.
TEST_F(WideDiamagneticDriftTest, GradientFormOpenIsConservativeWhenBoundaryValueUniform) {
  setCurvature(0.0, 1.0, 0.0);
  Options options = diamagneticOptions();
  options["diamagnetic_drift"]["radial_boundary"] = "open";
  options["diamagnetic_drift"]["form"] = "gradient";

  DiamagneticDrift component("diamagnetic_drift", options, nullptr);

  // With open radial boundaries, the gradient form's implied core flux vanishes
  // only when g = quantity * temperature / charge is uniform around the ring in
  // both the x guard cell and the first interior cell.
  const Field3D quantity = makeField<Field3D>([](Ind3D& i) {
    const bool on_boundary = i.x() <= mesh->xstart || i.x() >= mesh->xend
                             || i.y() <= mesh->ystart || i.y() >= mesh->yend;
    return on_boundary ? 1.0 : 1.5 + (0.1 * i.y()) + (0.05 * i.z());
  });
  const Field3D temperature{1.0};

  Options state;
  state["species"]["d+"]["charge"] = 1.0;
  state["species"]["d+"]["temperature"] = temperature;
  state["species"]["d+"]["density"] = quantity;
  component.declareAllSpecies({"d+"});
  component.transform(state);

  const auto source = get<Field3D>(state["species"]["d+"]["density_source"]);
  ASSERT_NEAR(volumeIntegral(source), 0.0, 1e-12);
}

// For divergence-free curvature, both closed modes conserve every source in both forms;
// open radial faces give a nonzero integral, and interior transport remains nonzero.
TEST_F(RadialDiamagneticDriftTest, RadialBoundaryModesConserveClosedSources) {
  // The fixture's constant curvature and J are discretely divergence-free.
  const auto quantity = makeField<Field3D>([](Ind3D& i) {
    return 1.0 + (0.2 * i.x() * i.x()) + (0.05 * i.y()) + (0.01 * i.z());
  });

  for (const auto* form : {"gradient", "divergence"}) {
    for (const auto* boundary : {"open", "closed", "closed_average"}) {
      SCOPED_TRACE(std::string(form) + ", " + boundary);
      Options options = diamagneticOptions();
      options["diamagnetic_drift"]["form"] = form;
      options["diamagnetic_drift"]["radial_boundary"] = boundary;
      DiamagneticDrift component("diamagnetic_drift", options, nullptr);
      auto state = diamagneticState(quantity);
      component.declareAllSpecies({"d+"});
      component.transform(state);

      for (const auto* name : {"density_source", "energy_source", "momentum_source"}) {
        SCOPED_TRACE(name);
        const auto source = get<Field3D>(state["species"]["d+"][name]);
        // The radial variation produces a nonzero source away from either wall.
        EXPECT_NE(source(mesh->xstart + 1, mesh->ystart, mesh->zstart), 0.0);
        if (std::string(boundary) == "open") {
          EXPECT_GT(std::abs(volumeIntegral(source)), 1e-6);
        } else {
          EXPECT_NEAR(volumeIntegral(source), 0.0, 1e-12);
        }
      }
    }
  }
}

// Both forms conserve with unequal target guard temperatures in identity and shifted
// coordinates, across radial modes and charge signs, without changing input temperature.
TEST_F(RadialDiamagneticDriftTest, SheathFacesConserveWithNonuniformGuardTemperature) {
  // These x indices lie on open field lines with targets at both ends.
  dynamic_cast<FakeMesh*>(mesh)->ix_separatrix = 0;
  const BoutReal wave_number = TWOPI / mesh->LocalNz;

  for (const bool shifted : {false, true}) {
#if !BOUT_HAS_FFTW
    if (shifted) {
      continue;
    }
#endif
    SCOPED_TRACE(shifted ? "shifted" : "identity");
    if (shifted) {
      const auto zshift = makeField<Field2D>([](Ind2D& i) { return 0.4 * i.y(); });
      mesh->getCoordinates()->setParallelTransform(
          bout::utils::make_unique<ShiftedMetric>(*mesh, CELL_CENTRE, zshift,
                                                  static_cast<BoutReal>(mesh->LocalNz)));
    }
    setCurvature(0.0, 1.0, 0.0);

    // Define profiles in the aligned frame, including nonuniform target values
    // and deliberately different guard temperatures. The shifted case has
    // different z shifts in a target cell and its adjacent guard cell.
    auto quantity_al = makeField<Field3D>([wave_number](Ind3D& i) {
      return 1.0 + 0.15 * i.x() + 0.2 * i.y() * i.y()
             + 0.25 * std::sin(wave_number * i.z());
    });
    quantity_al.setDirectionY(YDirectionType::Aligned);
    const auto quantity = fromFieldAligned(quantity_al, "RGN_ALL");
    auto temperature_al = makeField<Field3D>([wave_number](Ind3D& i) {
      if (i.y() < mesh->ystart) {
        return 5.0 + (0.3 * std::cos(wave_number * i.z()));
      }
      if (i.y() > mesh->yend) {
        return 9.0 + (0.4 * std::sin(wave_number * i.z()));
      }
      return 2.0 + (0.4 * i.y()) + (0.2 * std::cos(wave_number * i.z()));
    });
    temperature_al.setDirectionY(YDirectionType::Aligned);
    const auto temperature = fromFieldAligned(temperature_al, "RGN_ALL");
    const auto original_temperature = copy(temperature);

    for (const auto* form : {"gradient", "divergence"}) {
      for (const auto* boundary : {"open", "closed", "closed_average"}) {
        for (const BoutReal charge : {1.0, -2.0}) {
          SCOPED_TRACE(std::string(form) + ", " + boundary);
          SCOPED_TRACE(charge);
          Options options = diamagneticOptions();
          options["diamagnetic_drift"]["form"] = form;
          options["diamagnetic_drift"]["radial_boundary"] = boundary;
          DiamagneticDrift component("diamagnetic_drift", options, nullptr);
          auto state = diamagneticState(quantity, 1.0, charge);
          state["species"]["d+"]["temperature"] = temperature;
          component.declareAllSpecies({"d+"});
          component.transform(state);

          EXPECT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["d+"]["temperature"]),
                                   original_temperature, "RGN_ALL"));
          for (const auto* name :
               {"density_source", "energy_source", "momentum_source"}) {
            SCOPED_TRACE(name);
            const auto source = get<Field3D>(state["species"]["d+"][name]);
            EXPECT_NEAR(volumeIntegral(source), 0.0, 1e-12);
            EXPECT_NE(source(mesh->xstart + 1, mesh->ystart, mesh->zstart), 0.0);
            EXPECT_NE(source(mesh->xstart + 1, mesh->yend, mesh->zstart), 0.0);
          }
        }
      }
    }
  }
}

// Curvature forms reject transforms without aligned conversion for every radial mode,
// with an error explaining the required capability.
TEST_F(UnsupportedDiamagneticDriftTest,
       CurvatureFormsRejectUnsupportedParallelTransform) {
  for (const auto* form : {"gradient", "divergence"}) {
    for (const auto* boundary : {"open", "closed", "closed_average"}) {
      SCOPED_TRACE(boundary);
      Options options = diamagneticOptions();
      options["diamagnetic_drift"]["form"] = form;
      options["diamagnetic_drift"]["radial_boundary"] = boundary;
      try {
        const DiamagneticDrift component("diamagnetic_drift", options, nullptr);
        FAIL() << "curvature forms must reject a transform without aligned conversion";
      } catch (const BoutException& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find(std::string(form) + " form is not supported"),
                  std::string::npos);
        EXPECT_NE(message.find("field-aligned coordinates"), std::string::npos);
        EXPECT_NE(message.find("Use form = bracket"), std::string::npos);
      }
    }
  }
}
