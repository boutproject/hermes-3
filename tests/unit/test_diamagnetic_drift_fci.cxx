#include "gtest/gtest.h"

#include "test_extras.hxx" // FakeMesh                                                                                                                                                                     

#include "../../include/diamagnetic_drift_fci.hxx"

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
using DiamagneticDriftFCITest = FakeMeshFixture;

TEST_F(DiamagneticDriftFCITest, CreateComponent) {
  Options options;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["Tesla"] = 1.0;

  mesh->getCoordinates()->Bxy = 1.0;
  
  DiamagneticDriftFCI component("test", options, nullptr);
}

TEST_F(DiamagneticDriftFCITest, NoTemperature) {
  Options options;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["Tesla"] = 1.0;

  mesh->getCoordinates()->Bxy = 1.0;

  DiamagneticDriftFCI component("h+", options, nullptr);

  Field3D N = FieldFactory::get()->create3D("1 + y * (x - 0.5)", &options, mesh);

  Field3D T = FieldFactory::get()->create3D("1 + 0.2 * x + sin(z)", &options, mesh);
  
  mesh->communicate(N,T);

  Options state;
  state["species"]["h+"]["density"] = N;
  //state["species"]["h+"]["temperature"] = T;
  state["species"]["h+"]["charge"] = 1;
  
  component.transform(state);

  //ASSERT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["h+"]["density_source"]), 0.0, "RGN_NOBNDRY"));
  ASSERT_FALSE(state["species"]["h+"].isSet("density_source"));
  ASSERT_FALSE(state["species"]["h+"].isSet("momentum_source"));
  ASSERT_FALSE(state["species"]["h+"].isSet("energy_source"));
  //ASSERT_TRUE(IsFieldEqual(get<Field3D>(state["species"]["h+"]["energy_source"]), 0.0,"RGN_NOBNDRY"));
  
}


TEST_F(DiamagneticDriftFCITest, NoCharge) {
  Options options;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["Tesla"] = 1.0;

  mesh->getCoordinates()->Bxy = 1.0;

  DiamagneticDriftFCI component("h+", options, nullptr);

  Field3D N = FieldFactory::get()->create3D("1 + y * (x - 0.5)", &options, mesh);

  Field3D T = FieldFactory::get()->create3D("1 + 0.2 * x + sin(z)", &options, mesh);

  mesh->communicate(N,T);

  Options state;
  state["species"]["h+"]["density"] = N;
  state["species"]["h+"]["temperature"] = T;
  state["species"]["h+"]["pressure"] = N*T;
  state["species"]["h+"]["charge"] = 0;

  component.transform(state);

  ASSERT_FALSE(state["species"]["h+"].isSet("density_source"));
  ASSERT_FALSE(state["species"]["h+"].isSet("momentum_source"));
  ASSERT_FALSE(state["species"]["h+"].isSet("energy_source"));  
}



TEST_F(DiamagneticDriftFCITest, TwoSpeciesOppositeCharge) {
  Options options;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["Tesla"] = 1.0;


  Field3D B = FieldFactory::get()->create3D("1 + x", &options, mesh);
  mesh->getCoordinates()->Bxy = B;
  
  
  DiamagneticDriftFCI component("h+", options, nullptr);

  Field3D N = FieldFactory::get()->create3D("1 + y * (x - 0.5)", &options, mesh);  
  Field3D T = FieldFactory::get()->create3D("1 + 0.2 * x + sin(z)", &options, mesh);
  Field3D V = FieldFactory::get()->create3D("1 + 1.3 * x + cos(2.0*z)", &options, mesh);
  mesh->communicate(N,V,T);

  Options state;
  state["species"]["h+"]["density"] = N;
  state["species"]["h+"]["temperature"] = T;
  state["species"]["h+"]["pressure"] = N*T;
  state["species"]["h+"]["charge"] = 1;
  state["species"]["h+"]["velocity"] = V;
  state["species"]["h+"]["momentum"] = N*V;
  

  state["species"]["e"]["density"] = N;
  state["species"]["e"]["temperature"] = T;
  state["species"]["e"]["pressure"] = N*T;
  state["species"]["e"]["charge"] = -1;
  state["species"]["e"]["velocity"] = V;
  state["species"]["e"]["momentum"] = N*V;
  

  component.transform(state);

  Field3D ds1 = get<Field3D>(state["species"]["h+"]["density_source"]);
  Field3D ds2 =	get<Field3D>(state["species"]["e"]["density_source"]);

  Field3D ms1 = get<Field3D>(state["species"]["h+"]["momentum_source"]);
  Field3D ms2 = get<Field3D>(state["species"]["e"]["momentum_source"]);
  
  Field3D es1 = get<Field3D>(state["species"]["h+"]["energy_source"]);
  Field3D es2 = get<Field3D>(state["species"]["e"]["energy_source"]);


  
  BOUT_FOR_SERIAL(i, ds1.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(ds1[i],-ds2[i]);
    ASSERT_DOUBLE_EQ(ms1[i],-ms2[i]);
    ASSERT_DOUBLE_EQ(es1[i],-es2[i]);
  }  
}



TEST_F(DiamagneticDriftFCITest, NoGradientB) {
  Options options;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["Tesla"] = 1.0;


  Field3D B = FieldFactory::get()->create3D("1.0", &options, mesh);
  mesh->getCoordinates()->Bxy = B;


  DiamagneticDriftFCI component("h+", options, nullptr);

  Field3D N = FieldFactory::get()->create3D("1 + y * (x - 0.5)", &options, mesh);
  Field3D V = FieldFactory::get()->create3D("1 + 1.3 * x + cos(2.0*z)", &options, mesh);
  Field3D T = FieldFactory::get()->create3D("1 + 0.2 * x + sin(z)", &options, mesh);

  mesh->communicate(N,V,T);

  Options state;
  state["species"]["h+"]["density"] = N;
  state["species"]["h+"]["velocity"] = V;
  state["species"]["h+"]["temperature"] = T;
  state["species"]["h+"]["momentum"] = N*V;
  state["species"]["h+"]["pressure"] = N*T;
  state["species"]["h+"]["charge"] = 1;

  
  component.transform(state);

  Field3D ds1 = get<Field3D>(state["species"]["h+"]["density_source"]);
  Field3D ms1 = get<Field3D>(state["species"]["h+"]["momentum_source"]);
  Field3D es1 = get<Field3D>(state["species"]["h+"]["energy_source"]);


  BOUT_FOR_SERIAL(i, ds1.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(ds1[i],0.0);
    ASSERT_DOUBLE_EQ(ms1[i],0.0);
    ASSERT_DOUBLE_EQ(es1[i],0.0);
  }
  
}


TEST_F(DiamagneticDriftFCITest, DifferentCharge) {
  Options options;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["Tesla"] = 1.0;


  Field3D B = FieldFactory::get()->create3D("1 + x", &options, mesh);
  mesh->getCoordinates()->Bxy = B;


  DiamagneticDriftFCI component("h+", options, nullptr);

  Field3D N = FieldFactory::get()->create3D("1 + y * (x - 0.5)", &options, mesh);
  Field3D T = FieldFactory::get()->create3D("1 + 0.2 * x + sin(z)", &options, mesh);
  Field3D V = FieldFactory::get()->create3D("1 + 1.3 * x + cos(2.0*z)", &options, mesh);
  mesh->communicate(N,V,T);

  Options state;
  state["species"]["h+"]["density"] = N;
  state["species"]["h+"]["temperature"] = T;
  state["species"]["h+"]["pressure"] = N*T;
  state["species"]["h+"]["charge"] = 1;
  state["species"]["h+"]["velocity"] = V;
  state["species"]["h+"]["momentum"] = N*V;


  state["species"]["h++"]["density"] = N;
  state["species"]["h++"]["temperature"] = T;
  state["species"]["h++"]["pressure"] = N*T;
  state["species"]["h++"]["charge"] = 2;
  state["species"]["h++"]["velocity"] = V;
  state["species"]["h++"]["momentum"] = N*V;


  state["species"]["h++++"]["density"] = N;
  state["species"]["h++++"]["temperature"] = T;
  state["species"]["h++++"]["pressure"] = N*T;
  state["species"]["h++++"]["charge"] = 4;
  state["species"]["h++++"]["velocity"] = V;
  state["species"]["h++++"]["momentum"] = N*V;




  
  component.transform(state);

  Field3D ds1 = get<Field3D>(state["species"]["h+"]["density_source"]);
  Field3D ds2 = get<Field3D>(state["species"]["h++"]["density_source"]);
  Field3D ds3 = get<Field3D>(state["species"]["h++++"]["density_source"]);


  Field3D ms1 = get<Field3D>(state["species"]["h+"]["momentum_source"]);
  Field3D ms2 = get<Field3D>(state["species"]["h++"]["momentum_source"]);
  Field3D ms3 = get<Field3D>(state["species"]["h++++"]["momentum_source"]);
  
  Field3D es1 = get<Field3D>(state["species"]["h+"]["energy_source"]);
  Field3D es2 = get<Field3D>(state["species"]["h++"]["energy_source"]);
  Field3D es3 = get<Field3D>(state["species"]["h++++"]["energy_source"]);


  BOUT_FOR_SERIAL(i, ds1.getRegion("RGN_NOBNDRY")) {
    ASSERT_DOUBLE_EQ(ds1[i],2.0*ds2[i]);
    ASSERT_DOUBLE_EQ(ms1[i],2.0*ms2[i]);
    ASSERT_DOUBLE_EQ(es1[i],2.0*es2[i]);

    ASSERT_DOUBLE_EQ(ds1[i],4.0*ds3[i]);
    ASSERT_DOUBLE_EQ(ms1[i],4.0*ms3[i]);
    ASSERT_DOUBLE_EQ(es1[i],4.0*es3[i]);
    
  }
}




TEST_F(DiamagneticDriftFCITest, DoubleGradient) {
  Options options;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["Tesla"] = 1.0;


  Field3D B = FieldFactory::get()->create3D("1 + x", &options, mesh);
  mesh->getCoordinates()->Bxy = B;


  DiamagneticDriftFCI component("h+", options, nullptr);

  Field3D N = FieldFactory::get()->create3D("1 + y * (x - 0.5)", &options, mesh);
  Field3D T = FieldFactory::get()->create3D("1 + 0.2 * x + sin(z)", &options, mesh);
  Field3D V = FieldFactory::get()->create3D("1 + 1.3 * x + cos(2.0*z)", &options, mesh);
  mesh->communicate(N,V,T);

  Options state;
  state["species"]["h1"]["density"] = N;
  state["species"]["h1"]["temperature"] = T;
  state["species"]["h1"]["pressure"] = N*T;
  state["species"]["h1"]["charge"] = 1;
  state["species"]["h1"]["velocity"] = V;
  state["species"]["h1"]["momentum"] = N*V;

  // Multiply density by two to get the difference in gradient, then check if the diamagnetic drift scales accordingly
  state["species"]["h2"]["density"] = 2.0 * N;
  state["species"]["h2"]["temperature"] = T;
  state["species"]["h2"]["pressure"] = 2.0 * N * T;
  state["species"]["h2"]["charge"] = 1;
  state["species"]["h2"]["velocity"] = V;
  state["species"]["h2"]["momentum"] = 2.0 * N * V;

  component.transform(state);

  Field3D ds1 = get<Field3D>(state["species"]["h1"]["density_source"]);
  Field3D ds2 = get<Field3D>(state["species"]["h2"]["density_source"]);


  Field3D ms1 = get<Field3D>(state["species"]["h1"]["momentum_source"]);
  Field3D ms2 = get<Field3D>(state["species"]["h2"]["momentum_source"]);


  Field3D es1 = get<Field3D>(state["species"]["h1"]["energy_source"]);
  Field3D es2 = get<Field3D>(state["species"]["h2"]["energy_source"]);



  BOUT_FOR_SERIAL(i, ds1.getRegion("RGN_NOBNDRY")) {
    
    ASSERT_DOUBLE_EQ(ds1[i],0.5*ds2[i]);
    ASSERT_DOUBLE_EQ(ms1[i],0.5*ms2[i]);
    ASSERT_DOUBLE_EQ(es1[i],0.5*es2[i]);

  }
}

