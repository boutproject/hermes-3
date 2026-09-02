#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/sound_speed.hxx"

#include <memory>
#include <vector>

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

namespace {

/// Construct normalized units and component options for SoundSpeed tests.
Options makeSoundSpeedOptions(BoutReal Bnorm = 1.0, BoutReal Nnorm = 1.0,
                              BoutReal Lnorm = 1.0, BoutReal tnorm = 1.0,
                              BoutReal eVnorm = 1.0, bool alfven_wave = false,
                              bool electron_dynamics = true,
                              BoutReal fastest_wave_factor = 1.0) {
  Options options;
  options["units"]["Tesla"] = Bnorm;
  options["units"]["inv_meters_cubed"] = Nnorm;
  options["units"]["meters"] = Lnorm;
  options["units"]["seconds"] = tnorm;
  options["units"]["eV"] = eVnorm;

  Options& component = options["test"];
  component["alfven_wave"] = alfven_wave;
  component["electron_dynamics"] = electron_dynamics;
  component["fastest_wave_factor"] = fastest_wave_factor;

  return options;
}

/// Create a Field3D that varies only in Y, with constant X and Z planes.
Field3D makeYField(const std::vector<BoutReal>& yvalues) {
  EXPECT_EQ(static_cast<int>(yvalues.size()), mesh->LocalNy);

  Field3D result(0.0);
  for (int x = 0; x < mesh->LocalNx; ++x) {
    for (int y = 0; y < mesh->LocalNy; ++y) {
      for (int z = 0; z < mesh->LocalNz; ++z) {
        result(x, y, z) = yvalues[static_cast<size_t>(y)];
      }
    }
  }

  return result;
}

/// Calculate the normalization factor used for the Alfven-speed comparison.
///
/// Takes `Tesla`, `inv_meters_cubed`, `meters`, and `seconds` from `options["units"]`
/// and returns the factor that multiplies `Bxy / sqrt(total_density)` in
/// `SoundSpeed::transform_impl()`.
BoutReal getBetaNorm(const Options& options) {
  const auto& units = options["units"];
  return get<BoutReal>(units["Tesla"])
         / sqrt(SI::mu0 * get<BoutReal>(units["inv_meters_cubed"]) * SI::Mp)
         / (get<BoutReal>(units["meters"]) / get<BoutReal>(units["seconds"]));
}

/// FakeMesh variant whose Y communication wraps edge values into guard cells.
///
/// This lets the periodic-boundary test verify the post-communication state of
/// `sound_speed` and `fastest_wave` without needing a full distributed mesh.
class PeriodicCommFakeMesh : public FakeMesh {
public:
  PeriodicCommFakeMesh(int nx, int ny, int nz, MpiWrapper& mpi_in)
      : FakeMesh(nx, ny, nz, mpi_in) {}

  comm_handle send(FieldGroup& g) override {
    pending_fields = g.field3d();
    return reinterpret_cast<comm_handle>(this);
  }

  comm_handle sendY(FieldGroup& g, comm_handle UNUSED(handle) = nullptr) override {
    pending_fields = g.field3d();
    return reinterpret_cast<comm_handle>(this);
  }

  int wait(comm_handle UNUSED(handle)) override {
    for (auto* field : pending_fields) {
      for (int x = 0; x < LocalNx; ++x) {
        if (!periodicY(x)) {
          continue;
        }
        for (int z = 0; z < LocalNz; ++z) {
          (*field)(x, ystart - 1, z) = (*field)(x, yend, z);
          (*field)(x, yend + 1, z) = (*field)(x, ystart, z);
        }
      }
    }

    pending_fields.clear();
    return 0;
  }

private:
  std::vector<Field3D*> pending_fields;
};

/// Fixture with a simple metric and periodic Y communication.
class PeriodicCommSoundSpeedTest : public ::testing::Test {
public:
  PeriodicCommSoundSpeedTest()
      : mesh_m(3, 5, 7, mpi), quiet_info(output_info), quiet_warn(output_warn) {
    bout::globals::mpi = &mpi;
    bout::globals::mesh = &mesh_m;
    bout::globals::mesh->createDefaultRegions();
    mesh_m.setCoordinates(nullptr);

    test_coords = std::make_shared<Coordinates>(
        bout::globals::mesh, Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0},
        Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0},
        Field2D{0.0}, Field2D{0.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0},
        Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0});

    test_coords->G1 = test_coords->G2 = test_coords->G3 = 0.1;
    test_coords->non_uniform = true;
    test_coords->d1_dx = test_coords->d1_dy = 0.2;
    test_coords->d1_dz = 0.0;
#if BOUT_USE_METRIC_3D
    test_coords->Bxy.splitParallelSlices();
    test_coords->Bxy.yup() = test_coords->Bxy.ydown() = test_coords->Bxy;
#endif

    mesh_m.setGridDataSource(new FakeGridDataSource());
    test_coords->setParallelTransform(
        bout::utils::make_unique<ParallelTransformIdentity>(*bout::globals::mesh));
    mesh_m.setCoordinates(test_coords);
    mesh_m.createBoundaryRegions();
  }

  ~PeriodicCommSoundSpeedTest() override {
    bout::globals::mesh = nullptr;
    bout::globals::mpi = nullptr;

    Options::cleanup();
  }

private:
  PeriodicCommFakeMesh mesh_m;
  std::shared_ptr<Coordinates> test_coords{nullptr};
  MpiWrapper mpi;
  WithQuietOutput quiet_info;
  WithQuietOutput quiet_warn;
};

} // namespace

// Reuse the "standard" fixture for FakeMesh
using SoundSpeedTest = FakeMeshFixture;

/// Baseline case: with one species, both outputs reduce to that species'
/// collective sound speed in the interior region (`RGN_NOBNDRY`).
TEST_F(SoundSpeedTest, OneSpeciesSetsSoundSpeedAndFastestWaveInInteriorCells) {
  Options options = makeSoundSpeedOptions();
  SoundSpeed component("test", options, nullptr);

  options["species"]["e"]["density"] = 2.0;
  options["species"]["e"]["pressure"] = 1.2;
  options["species"]["e"]["AA"] = 1.5;

  component.declareAllSpecies({"e"});
  component.transform(options);

  const BoutReal expected = sqrt(1.2 / (2.0 * 1.5));

  ASSERT_TRUE(options.isSet("sound_speed"));
  ASSERT_TRUE(
      IsFieldEqual(get<Field3D>(options["sound_speed"]), expected, "RGN_NOBNDRY"));
  ASSERT_TRUE(
      IsFieldEqual(get<Field3D>(options["fastest_wave"]), expected, "RGN_NOBNDRY"));
}

/// With multiple species and no temperature-driven species maximum, both outputs
/// use the collective sound speed sqrt(sum(pressure) / sum(density * AA)).
TEST_F(SoundSpeedTest, TwoSpeciesCollectiveSoundSpeed) {
  Options options = makeSoundSpeedOptions();
  SoundSpeed component("test", options, nullptr);

  options["species"]["e"]["density"] = 2.0;
  options["species"]["e"]["pressure"] = 1.2;
  options["species"]["e"]["AA"] = 1.5;

  options["species"]["h"]["density"] = 3.0;
  options["species"]["h"]["pressure"] = 2.5;
  options["species"]["h"]["AA"] = 0.9;

  component.declareAllSpecies({"e", "h"});
  component.transform(options);

  const BoutReal expected = sqrt((1.2 + 2.5) / (2.0 * 1.5 + 3.0 * 0.9));

  ASSERT_TRUE(options.isSet("sound_speed"));
  ASSERT_TRUE(
      IsFieldEqual(get<Field3D>(options["sound_speed"]), expected, "RGN_NOBNDRY"));
  ASSERT_TRUE(
      IsFieldEqual(get<Field3D>(options["fastest_wave"]), expected, "RGN_NOBNDRY"));
}

/// `sound_speed` stays equal to the collective value, while `fastest_wave`
/// tracks the pointwise maximum of the collective and per-species thermal speeds.
TEST_F(SoundSpeedTest, FastestWaveIsMaximumOfSpeciesAndCollectiveSoundSpeeds) {
  Options options = makeSoundSpeedOptions();
  SoundSpeed component("test", options, nullptr);

  options["species"]["h+"]["density"] = 1.0;
  options["species"]["h+"]["pressure"] = 8.0;
  options["species"]["h+"]["AA"] = 1.0;
  // Vary temperature while keeping pressure fixed so only the per-species thermal
  // speed changes; this makes `fastest_wave` exceed the collective sound speed
  // in selected interior cells.
  options["species"]["h+"]["temperature"] = makeYField({1.0, 1.0, 4.0, 9.0, 9.0});

  options["species"]["he+"]["density"] = 1.0;
  options["species"]["he+"]["pressure"] = 12.0;
  options["species"]["he+"]["AA"] = 4.0;
  options["species"]["he+"]["temperature"] = 1.0;

  component.declareAllSpecies({"h+", "he+"});
  component.transform(options);

  ASSERT_TRUE(IsFieldEqual(get<Field3D>(options["sound_speed"]), 2.0, "RGN_NOBNDRY"));
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(options["fastest_wave"]),
                           makeYField({0.0, 2.0, 2.0, 3.0, 0.0}), "RGN_NOBNDRY"));
}

/// When `alfven_wave` is enabled but the Alfven speed is below the collective
/// sound speed, `fastest_wave` should remain equal to `sound_speed`.
TEST_F(SoundSpeedTest, AlfvenWaveUsesSoundSpeedAtHighBeta) {
  Options options = makeSoundSpeedOptions(1.0, 1e38, 1e19, 1.0, 1.0, true);
  SoundSpeed component("test", options, nullptr);

  mesh->getCoordinates()->Bxy = 1.0;
  options["species"]["i"]["density"] = 1.0;
  options["species"]["i"]["pressure"] = 4.0;
  options["species"]["i"]["AA"] = 1.0;

  component.declareAllSpecies({"i"});
  component.transform(options);

  const BoutReal sound_speed = 2.0;
  const BoutReal alfven_speed = getBetaNorm(options);

  ASSERT_LT(alfven_speed, sound_speed);
  ASSERT_TRUE(
      IsFieldEqual(get<Field3D>(options["sound_speed"]), sound_speed, "RGN_NOBNDRY"));
  ASSERT_TRUE(
      IsFieldEqual(get<Field3D>(options["fastest_wave"]), sound_speed, "RGN_NOBNDRY"));
}

/// When the Alfven speed exceeds the collective sound speed, `fastest_wave`
/// should switch to the Alfven branch while `sound_speed` stays unchanged.
TEST_F(SoundSpeedTest, AlfvenWaveUsesAlfvenSpeedAtLowBeta) {
  Options options = makeSoundSpeedOptions(1.0, 1.0, 1.0, 1.0, 1.0, true);
  SoundSpeed component("test", options, nullptr);

  mesh->getCoordinates()->Bxy = 1.0;
  options["species"]["i"]["density"] = 1.0;
  options["species"]["i"]["pressure"] = 1.0;
  options["species"]["i"]["AA"] = 1.0;

  component.declareAllSpecies({"i"});
  component.transform(options);

  const BoutReal sound_speed = 1.0;
  const BoutReal alfven_speed = getBetaNorm(options);

  ASSERT_GT(alfven_speed, sound_speed);
  ASSERT_TRUE(
      IsFieldEqual(get<Field3D>(options["sound_speed"]), sound_speed, "RGN_NOBNDRY"));
  ASSERT_TRUE(
      IsFieldEqual(get<Field3D>(options["fastest_wave"]), alfven_speed, "RGN_NOBNDRY"));
}

/// Disabling `electron_dynamics` skips the electron species when building the
/// per-species maximum for `fastest_wave`, but still includes electron pressure
/// in the collective `sound_speed`.
TEST_F(SoundSpeedTest, ElectronDynamicsFalseExcludesElectronSpeciesSoundSpeed) {
  Options options = makeSoundSpeedOptions(1.0, 1.0, 1.0, 1.0, 1.0, false, false);
  SoundSpeed component("test", options, nullptr);

  options["species"]["e"]["density"] = 1.0;
  options["species"]["e"]["pressure"] = 8.0;
  options["species"]["e"]["AA"] = 5e-4;
  options["species"]["e"]["temperature"] = 100.0;

  options["species"]["i"]["density"] = 1.0;
  options["species"]["i"]["pressure"] = 1.0;
  options["species"]["i"]["AA"] = 1.0;
  options["species"]["i"]["temperature"] = 1.0;

  component.declareAllSpecies({"e", "i"});
  component.transform(options);

  ASSERT_TRUE(IsFieldEqual(get<Field3D>(options["sound_speed"]), 3.0, "RGN_NOBNDRY"));
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(options["fastest_wave"]), 3.0, "RGN_NOBNDRY"));
}

/// `fastest_wave_factor` is applied only to the final `fastest_wave` output and
/// must not change the separately reported `sound_speed`.
TEST_F(SoundSpeedTest, FastestWaveFactorScalesFastestWaveOnly) {
  Options options = makeSoundSpeedOptions(1.0, 1.0, 1.0, 1.0, 1.0, false, true, 1.7);
  SoundSpeed component("test", options, nullptr);

  options["species"]["i"]["density"] = 1.0;
  options["species"]["i"]["pressure"] = 4.0;
  options["species"]["i"]["AA"] = 1.0;

  component.declareAllSpecies({"i"});
  component.transform(options);

  ASSERT_TRUE(IsFieldEqual(get<Field3D>(options["sound_speed"]), 2.0, "RGN_NOBNDRY"));
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(options["fastest_wave"]), 3.4, "RGN_NOBNDRY"));
}

/// On a non-periodic mesh, Y boundary values should be overwritten by boundary
/// handling rather than preserving the supplied edge-cell inputs.
TEST_F(SoundSpeedTest, NonPeriodicBoundaryCellsAreSet) {
  auto* fake_mesh = dynamic_cast<FakeMesh*>(mesh);
  ASSERT_NE(fake_mesh, nullptr);
  fake_mesh->ix_separatrix = 0;

  Options options = makeSoundSpeedOptions();
  SoundSpeed component("test", options, nullptr);

  // Seed distinct edge values so the test can tell whether boundary handling
  // overwrites them with the expected interior-adjacent values.
  const Field3D input = makeYField({25.0, 4.0, 9.0, 16.0, 36.0});

  options["species"]["i"]["density"] = 1.0;
  options["species"]["i"]["pressure"] = input;
  options["species"]["i"]["AA"] = 1.0;
  options["species"]["i"]["temperature"] = input;

  component.declareAllSpecies({"i"});
  component.transform(options);

  const auto sound_speed = get<Field3D>(options["sound_speed"]);
  const auto fastest_wave = get<Field3D>(options["fastest_wave"]);

  for (int x = mesh->xstart; x <= mesh->xend; ++x) {
    for (int z = 0; z < mesh->LocalNz; ++z) {
      EXPECT_DOUBLE_EQ(sound_speed(x, mesh->ystart, z), 2.0);
      EXPECT_DOUBLE_EQ(sound_speed(x, mesh->yend, z), 4.0);
      EXPECT_DOUBLE_EQ(fastest_wave(x, mesh->ystart, z), 2.0);
      EXPECT_DOUBLE_EQ(fastest_wave(x, mesh->yend, z), 4.0);
    }
  }
}

/// On a periodic mesh, Y guard cells should be populated by communication from
/// the opposite interior edge after `sound_speed` and `fastest_wave` are computed.
TEST_F(PeriodicCommSoundSpeedTest, PeriodicGuardCellsAreCommunicated) {
  Options options = makeSoundSpeedOptions();
  SoundSpeed component("test", options, nullptr);

  const Field3D input = makeYField({100.0, 4.0, 9.0, 16.0, 121.0});

  options["species"]["i"]["density"] = 1.0;
  options["species"]["i"]["pressure"] = input;
  options["species"]["i"]["AA"] = 1.0;
  options["species"]["i"]["temperature"] = input;

  component.declareAllSpecies({"i"});
  component.transform(options);

  const auto sound_speed = get<Field3D>(options["sound_speed"]);
  const auto fastest_wave = get<Field3D>(options["fastest_wave"]);

  ASSERT_TRUE(
      IsFieldEqual(sound_speed, makeYField({0.0, 2.0, 3.0, 4.0, 0.0}), "RGN_NOBNDRY"));
  ASSERT_TRUE(
      IsFieldEqual(fastest_wave, makeYField({0.0, 2.0, 3.0, 4.0, 0.0}), "RGN_NOBNDRY"));

  for (int x = mesh->xstart; x <= mesh->xend; ++x) {
    for (int z = 0; z < mesh->LocalNz; ++z) {
      EXPECT_DOUBLE_EQ(sound_speed(x, mesh->ystart - 1, z), 4.0);
      EXPECT_DOUBLE_EQ(sound_speed(x, mesh->yend + 1, z), 2.0);
      EXPECT_DOUBLE_EQ(fastest_wave(x, mesh->ystart - 1, z), 4.0);
      EXPECT_DOUBLE_EQ(fastest_wave(x, mesh->yend + 1, z), 2.0);
    }
  }
}
