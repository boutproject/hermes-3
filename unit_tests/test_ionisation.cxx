
#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh

#include "../include/ionisation.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using IonisationTest = FakeMeshFixture;

TEST_F(IonisationTest, CreateComponent) {
  Options options;
  
  options["units"]["eV"] = 1.0;
  options["units"]["inv_meters_cubed"] = 1.0;
  options["units"]["seconds"] = 1.0;
  Ionisation component("test", options, nullptr);
}

TEST_F(IonisationTest, MissingData) {
  Options options;
  
  options["units"]["eV"] = 1.0;
  options["units"]["inv_meters_cubed"] = 1.0;
  options["units"]["seconds"] = 1.0;
  Ionisation component("test", options, nullptr);

  ASSERT_THROW(component.transform(options), BoutException);
}

TEST_F(IonisationTest, SourcesSanityCheck) {
  Options options;
  
  options["units"]["eV"] = 1.0;
  options["units"]["inv_meters_cubed"] = 1.0;
  options["units"]["seconds"] = 1.0;
  Ionisation component("test", options, nullptr);

  // Set atom properties
  options["species"]["h"]["density"] = 1e19;
  options["species"]["h"]["temperature"] = 3.5;
  options["species"]["h"]["velocity"] = 10.0;
  options["species"]["h"]["AA"] = 2.0;

  // Ion properties
  options["species"]["h+"]["AA"] = 2.0;
  
  // Electron properties
  options["species"]["e"]["density"] = 1e18;
  options["species"]["e"]["temperature"] = 5.;

  // Do the calculation
  component.transform(options);
  
  Field3D h_density_source = get<Field3D>(options["species"]["h"]["density_source"]);
  Field3D hp_density_source = get<Field3D>(options["species"]["h+"]["density_source"]);
  
  BOUT_FOR(i, h_density_source.getRegion("RGN_NOBNDRY")) {
    ASSERT_LT(h_density_source[i], 0.0) << "Atom (h) density source not negative at " << i;
    ASSERT_GT(hp_density_source[i], 0.0) << "Ion (h+) density source not positive at " << i;
  }

  Field3D h_momentum_source = get<Field3D>(options["species"]["h"]["momentum_source"]);
  Field3D hp_momentum_source = get<Field3D>(options["species"]["h+"]["momentum_source"]);

  BOUT_FOR(i, h_momentum_source.getRegion("RGN_NOBNDRY")) {
    ASSERT_LT(h_momentum_source[i], 0.0) << "Atom (h) momentum source not negative at " << i;
    ASSERT_GT(hp_momentum_source[i], 0.0) << "Ion (h+) momentum source not positive at " << i;
  }
  
  Field3D h_energy_source = get<Field3D>(options["species"]["h"]["energy_source"]);
  Field3D hp_energy_source = get<Field3D>(options["species"]["h+"]["energy_source"]);

  BOUT_FOR(i, h_energy_source.getRegion("RGN_NOBNDRY")) {
    ASSERT_LT(h_energy_source[i], 0.0) << "Atom (h) energy source not negative at " << i;
    ASSERT_GT(hp_energy_source[i], 0.0) << "Ion (h+) energy source not positive at " << i;
  }

  Field3D e_energy_source = get<Field3D>(options["species"]["e"]["energy_source"]);
  BOUT_FOR(i, h_energy_source.getRegion("RGN_NOBNDRY")) {
    ASSERT_LT(e_energy_source[i], 0.0) << "Electron (e) energy source not negative at " << i;
  }
}
