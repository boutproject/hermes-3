
#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh

#include "../../include/anomalous_diffusion.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

#include <field_factory.hxx>  // For generating functions

// Reuse the "standard" fixture for FakeMesh
using AnomalousDiffusionTest = FakeMeshFixture;

TEST_F(AnomalousDiffusionTest, CreateComponent) {
  Options options;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;

  AnomalousDiffusion component("test", options, nullptr);
}

TEST_F(AnomalousDiffusionTest, NoDiffusion) {
  Options options;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;

  AnomalousDiffusion component("h", options, nullptr);

  Field3D N = FieldFactory::get()->create3D("1 + y * (x - 0.5)", &options, mesh);
  mesh->communicate(N);
  
  Options state;
  state["species"]["h"]["density"] = N;
  
  // If D is not set, then the diffusion should not be calculated
  component.transform(state);

  ASSERT_FALSE(state["species"]["h"].isSet("density_source"));
  ASSERT_FALSE(state["species"]["h"].isSet("momentum_source"));
  ASSERT_FALSE(state["species"]["h"].isSet("energy_source"));
}

TEST_F(AnomalousDiffusionTest, ParticleDiffusion) {
  
  Coordinates *coords = mesh->getCoordinates();
  coords->Bxy = 1.0; // Note: This is non-finite or zero?
  
  Options options;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;

  options["h"]["anomalous_D"] = 1.0;  // Set particle diffusion for "h" species

  AnomalousDiffusion component("h", options, nullptr);

  Options state;
  state["species"]["h"]["density"] =
    FieldFactory::get()->create3D("1 + y * (x - 0.5)", &options, mesh);
  state["species"]["h"]["AA"] = 1.0; // Atomic mass number
  
  component.transform(state);

  // Expect all sources to be set
  ASSERT_TRUE(state["species"]["h"].isSet("density_source"));
  ASSERT_TRUE(state["species"]["h"].isSet("momentum_source"));
  ASSERT_TRUE(state["species"]["h"].isSet("energy_source"));

  // Expect momentum and energy sources to be zero
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["h"]["momentum_source"]), 0.0,
                           "RGN_NOBNDRY"));
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["h"]["energy_source"]), 0.0,
                           "RGN_NOBNDRY"));

  // Expect the sum over all cells of density source to be zero
  Field2D dV = coords->J * coords->dx * coords->dy * coords->dz; // Cell volume

  Field3D source = get<Field3D>(state["species"]["h"]["density_source"]);
  BoutReal integral = 0.0;
  BOUT_FOR_SERIAL(i, source.getRegion("RGN_NOBNDRY")) {
    ASSERT_TRUE(std::isfinite(dV[i])) << "Volume element not finite at " << i.x() << ", " << i.y() << ", " << i.z();
    ASSERT_TRUE(std::isfinite(source[i])) << "Density source not finite at " << i.x() << ", " << i.y() << ", " << i.z();
    integral += source[i] * dV[i];
  }
  ASSERT_LT(abs(integral), 1e-3) << "Integral of density source should be close to zero";
}
