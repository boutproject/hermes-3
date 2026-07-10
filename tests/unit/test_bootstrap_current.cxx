#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/bootstrap_current.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using BootstrapCurrentTest = FakeMeshFixture;

TEST_F(BootstrapCurrentTest, CreateComponent) {
  Options options;

  options["units"]["eV"] = 1.0;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["inv_meters_cubed"] = 1e19;
  options["units"]["Tesla"] = 1.0;

  BootstrapCurrent component("test", options, nullptr);
}

TEST_F(BootstrapCurrentTest, L31) {
  Options options = {{"units",
                      {{"eV", 1.0},
                       {"meters", 1.0},
                       {"seconds", 1.0},
                       {"inv_meters_cubed", 1e19},
                       {"Tesla", 1.0}}}};

  BootstrapCurrent component("test", options, nullptr);

  // Values extracted from Fig 6 using Webplotdigitizer
  const std::vector<std::pair<BoutReal, BoutReal>> data_pairs = {
      {0.000024593943752273463, 0.7587339243351231},
      {0.0010156878372616458, 0.7476463978024723},
      {0.2568197705955592, 0.5704600181938178},
      {1.1867109054054632, 0.42775092306869067},
      {8.724594138300777, 0.2019174856857464},
      {47.253903077019295, 0.06978220930025159},
      {516.3968777088756, 0.013937891300858096},
  };

  const BoutReal Z = 1.0;
  const BoutReal f_t = 0.65;

  for (const auto& [nu_e, expected] : data_pairs) {
    EXPECT_NEAR(component.L31(Z, f_t, nu_e), expected, 0.02);
  }
}

TEST_F(BootstrapCurrentTest, L32) {
  Options options = {{"units",
                      {{"eV", 1.0},
                       {"meters", 1.0},
                       {"seconds", 1.0},
                       {"inv_meters_cubed", 1e19},
                       {"Tesla", 1.0}}}};

  BootstrapCurrent component("test", options, nullptr);

  // Figure 5
  const std::vector<std::pair<BoutReal, BoutReal>> data_pairs = {
      {0.000021768967843509877, -0.24141540384314963},
      {0.0012826498305280611, -0.22883143258865835},
      {0.012826498305280598, -0.19357131698172203},
      {0.07796360130405221, -0.12375410092173134},
      {0.5366976945540476, -0.003038587720669006},
      {3.0653952950565238, 0.07138728323699384},
      {13.65007806546011, 0.0440165599125133},
      {80.42765483057629, 0.016716138103420985},
      {710.1520602780329, 0.005702234025933062},
  };

  const BoutReal Z = 1.0;
  const BoutReal f_t = 0.453;

  for (const auto& [nu_e, expected] : data_pairs) {
    EXPECT_NEAR(component.L32(Z, f_t, nu_e), expected, 0.02);
  }
}

TEST_F(BootstrapCurrentTest, alpha) {
  Options options = {{"units",
                      {{"eV", 1.0},
                       {"meters", 1.0},
                       {"seconds", 1.0},
                       {"inv_meters_cubed", 1e19},
                       {"Tesla", 1.0}}}};

  BootstrapCurrent component("test", options, nullptr);

  // Figure 5
  const std::vector<std::pair<BoutReal, BoutReal>> data_pairs = {
      {0.0001654967074946295, -0.9455445544554455},
      {0.010038493711203604, -0.8811881188118811},
      {0.43532520353253096, -0.5990099009900987},
      {5.495847306930358, -0.19306930693069324},
      {34.25125177003207, 0.1386138613861383},
      {105.52605861234818, 0.6633663366336631},
      {300.6909621346027, 1.5742574257425739},
      {926.0187529908255, 2.0198019801980194},
  };

  const BoutReal f_t = 0.24;

  for (const auto& [nu_i, expected] : data_pairs) {
    EXPECT_NEAR(component.alpha(f_t, nu_i), expected, 0.04);
  }
}
