
#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh

#include "../../include/fixed_fraction_ions.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using FixedFractionIonsTest = FakeMeshFixture;

TEST_F(FixedFractionIonsTest, CreateComponent) {
  Options options;
  options["test"]["fractions"] = "h @ 1";
    
  FixedFractionIons component("test", options, nullptr);
}

TEST_F(FixedFractionIonsTest, MissingFractions) {
  Options options;

  ASSERT_THROW(FixedFractionIons component("test", options, nullptr),
               BoutException);
}

TEST_F(FixedFractionIonsTest, MalformedFraction) {
  Options options;
  options["test"]["fractions"] = "h 1";

  ASSERT_THROW(FixedFractionIons component("test", options, nullptr),
               BoutException);
}

TEST_F(FixedFractionIonsTest, TrailingCommaNoThrow) {
  Options options;
  options["test"]["fractions"] = "h @ 1,";
  
  FixedFractionIons component("test", options, nullptr);
}

TEST_F(FixedFractionIonsTest, MissingElecDensity) {
  Options options;
  options["test"]["fractions"] = "h @ 1";
    
  FixedFractionIons component("test", options, nullptr);

  Options state;
  ASSERT_THROW(component.transform(state), BoutException);
}

TEST_F(FixedFractionIonsTest, SingleIon) {
  Options options;
  options["test"]["fractions"] = "h @ 0.6";
    
  FixedFractionIons component("test", options, nullptr);

  Options state;
  state["species"]["e"]["density"] = 1.4;
  component.transform(state);

  ASSERT_TRUE(IsFieldEqual(state["species"]["h"]["density"].as<Field3D>(),
                           0.6 * 1.4, "RGN_NOBNDRY"));
}

TEST_F(FixedFractionIonsTest, TwoIons) {
  Options options;
  options["test"]["fractions"] = "some_thing @ 0.4, ne @ 0.05";
    
  FixedFractionIons component("test", options, nullptr);

  Options state;
  state["species"]["e"]["density"] = 1.4;
  component.transform(state);

  ASSERT_TRUE(IsFieldEqual(state["species"]["some_thing"]["density"].as<Field3D>(),
                           0.4 * 1.4, "RGN_NOBNDRY"));
  
  ASSERT_TRUE(IsFieldEqual(state["species"]["ne"]["density"].as<Field3D>(),
                           0.05 * 1.4, "RGN_NOBNDRY"));
}
