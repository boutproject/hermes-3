#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh                                                                                                                                                                     

#include "../../include/adhoc_potential.hxx"

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
using AdhocPotentialTest = FakeMeshFixture;

TEST_F(AdhocPotentialTest, CreateComponent) {
  Options options;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["Tesla"] = 1.0;


  
  AdhocPotential component("test", options, nullptr);
}


TEST_F(AdhocPotentialTest, SetPhiDefault) {
  Options options;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["Tesla"] = 1.0;

  AdhocPotential component("test", options, nullptr);

  Field3D T = FieldFactory::get()->create3D("1 + 0.2 * x + sin(z)", &options, mesh);
  
  mesh->communicate(T);

  Options state;
  state["species"]["e"]["temperature"] = T;
  //state["species"]["h+"]["temperature"] = T;
  state["species"]["e"]["charge"] = -1;
  
  component.transform(state);


  Field3D thisphi = get<Field3D>(state["fields"]["phi"]);

  
  BOUT_FOR_SERIAL(i, T.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(thisphi[i] ,3.0 * T[i]);    
  }
  
}
