
#ifndef __RADIATION_H__
#define __RADIATION_H__

#include <bout/field3d.hxx>
#include <bout/bout_types.hxx>

#include <vector>
#include <cmath>

class RadiatedPower {
public:
  const Field3D power(const Field3D &Te, const Field3D &Ne, const Field3D &Ni);
  
  virtual BoutReal power(BoutReal Te, BoutReal ne, BoutReal ni) = 0;
  
private:
};

class InterpRadiatedPower : public RadiatedPower {
public:
  InterpRadiatedPower(const std::string &file);
  
  BoutReal power(BoutReal, BoutReal, BoutReal) {return 0.0;}
  
private:
  std::vector<BoutReal> te_array;  // Te in eV
  std::vector<BoutReal> p_array;   // Radiative loss rate in Watts m^3
};


/// Rates supplied by Eva Havlicova
class HydrogenRadiatedPower : public RadiatedPower {
public:
  BoutReal power(BoutReal, BoutReal, BoutReal) {return 0.0;}
  
  // Collision rate coefficient <sigma*v> [m3/s]
  BoutReal ionisation(BoutReal Te);
  
  //<sigma*v> [m3/s]
  BoutReal recombination(BoutReal n, BoutReal Te);
  
  // <sigma*v> [m3/s]
  BoutReal chargeExchange(BoutReal Te);
  
  // <sigma*v> [m3/s]
  BoutReal excitation(BoutReal Te);
  
private:
  
};

/*!
 * Hydrogen rates, fitted by Hannah Willett May 2015
 * University of York
 */
class UpdatedRadiatedPower : public RadiatedPower {
public:
  BoutReal power(BoutReal, BoutReal, BoutReal) {return 0.0;}

  // Ionisation rate coefficient <sigma*v> [m3/s]
  BoutReal ionisation(BoutReal T); 
  
  // Recombination rate coefficient <sigma*v> [m3/s]
  BoutReal recombination(BoutReal n, BoutReal T);
  
  // Charge exchange rate coefficient <sigma*v> [m3/s]
  BoutReal chargeExchange(BoutReal Te);
  
  BoutReal excitation(BoutReal Te);
private:
  
};


/// Carbon in coronal equilibrium 
/// From I.H.Hutchinson Nucl. Fusion 34 (10) 1337 - 1348 (1994)
class HutchinsonCarbonRadiation : public RadiatedPower {
  BoutReal power(BoutReal Te, BoutReal ne, BoutReal ni) {
    if((Te < 0.0) || (ne < 0.0) || (ni < 0.0))
       return 0.0;
    return ne * ni * 2e-31*pow(Te/10., 3) / (1. + pow(Te/10., 4.5));
  }
};

/*
class PostJensen : public RadiatedPower {
public:
  BoutReal power(BoutReal Te, BoutReal ne, BoutReal ni) {
    if( (Te < Tmin) || (Te > Tmax) )
      return 0.0;
    
    return 0.0;
  }
protected:
  BoutReal Tmin, Tmax;
  BoutReal A[6];
  
  struct PJ_Data {
    const char* label; // Short name
    const char* name;  // Long name
    BoutReal Tmin, Tmax;
    BoutReal data[6];
  };
  
  static PJ_Data power_data[] = {
    {"C", "Carbon"},
    {0}
  };
};
*/

#endif // __RADIATION_H__
