
#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh

#include "../../include/sound_speed.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using SoundSpeedTest = FakeMeshFixture;

TEST_F(SoundSpeedTest, CreateComponent) {
  Options options;
  
  SoundSpeed component("test", options, nullptr);
}

TEST_F(SoundSpeedTest, OneSpecies) {
  Options options;
  SoundSpeed component("test", options, nullptr);
  
  options["species"]["e"]["density"] = 2.0;
  options["species"]["e"]["pressure"] = 1.2;
  options["species"]["e"]["AA"] = 1.5;
  
  component.transform(options);

  ASSERT_TRUE(options.isSet("sound_speed"));
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(options["sound_speed"]), sqrt(1.2 / (2.0 * 1.5)),
                           "RGN_NOBNDRY"));
}

TEST_F(SoundSpeedTest, TwoSpecies) {
  Options options;
  SoundSpeed component("test", options, nullptr);
  
  options["species"]["e"]["density"] = 2.0;
  options["species"]["e"]["pressure"] = 1.2;
  options["species"]["e"]["AA"] = 1.5;

  options["species"]["h"]["density"] = 3.0;
  options["species"]["h"]["pressure"] = 2.5;
  options["species"]["h"]["AA"] = 0.9;
  
  component.transform(options);

  ASSERT_TRUE(options.isSet("sound_speed"));
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(options["sound_speed"]),
                           sqrt((1.2 + 2.5) / (2.0 * 1.5 + 3.0 * 0.9)), "RGN_NOBNDRY"));
}
