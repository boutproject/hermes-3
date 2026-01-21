
#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh
#include "fake_mesh_fixture.hxx"

#include "../../include/classical_diffusion.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

#include <bout/field_factory.hxx>  // For generating functions

// Reuse the "standard" fixture for FakeMesh
using ClassicalDiffusionTest = FakeMeshFixture;

TEST_F(ClassicalDiffusionTest, CreateComponent) {
  Options options;
  ClassicalDiffusion component("test", options, nullptr);
}

TEST_F(ClassicalDiffusionTest, ParticleDiffusion) {
  
  Coordinates *coords = mesh->getCoordinates();
  coords->Bxy = 1.0; // Note: This is non-finite or zero?
  
  Options options;
  options["classical_diffusion"]["custom_D"] = -1.0;  // Set particle diffusion for "h" species

  ClassicalDiffusion component("classical_diffusion", options, nullptr);

  Options state;
  // an ion
  state["species"]["h"]["density"] =
    FieldFactory::get()->create3D("1 + y * (x - 0.5)", &options, mesh);
  state["species"]["h"]["pressure"] =
    FieldFactory::get()->create3D("1", &options, mesh);
  state["species"]["h"]["velocity"] =
    FieldFactory::get()->create3D("1", &options, mesh);
  state["species"]["h"]["temperature"] =
    FieldFactory::get()->create3D("1", &options, mesh);
  state["species"]["h"]["AA"] = 1.0; // Atomic mass number
  state["species"]["h"]["charge"] = 1.0; // Atomic mass number
  // the electrons
  state["species"]["e"]["density"] =
    FieldFactory::get()->create3D("1 + y * (x - 0.5)", &options, mesh);
  state["species"]["e"]["pressure"] =
    FieldFactory::get()->create3D("1", &options, mesh);
  state["species"]["e"]["velocity"] =
    FieldFactory::get()->create3D("1", &options, mesh);
  state["species"]["e"]["temperature"] =
    FieldFactory::get()->create3D("1", &options, mesh);
  state["species"]["e"]["AA"] = 1.0/1836.0; // Atomic mass number
  state["species"]["e"]["charge"] = -1.0; // Atomic mass number
  state["species"]["e"]["collision_frequency"] = 1.0; // nu_ei
  
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
