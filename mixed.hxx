/*
 * Mixed diffusive-fluid model, similar to UEDGE model
 * 3D model, diffusive in X-Z and fluid in Y
 */

#ifndef __NEUTRAL_MIXED_H__
#define __NEUTRAL_MIXED_H__

#include "neutral-model.hxx"

class NeutralMixed : public NeutralModel {
public:
  NeutralMixed(Solver *solver, Mesh *mesh, Options &options);
  ~NeutralMixed() {}

  /// Update plasma quantities
  void update(const Field3D &Ne, const Field3D &Te, const Field3D &Ti, const Field3D &Vi);
  
  void addDensity(int x, int y, int z, BoutReal dndt) {ddt(Nn)(x,y,z) += dndt;}
  void addPressure(int x, int y, int z, BoutReal dpdt) {ddt(Pn)(x,y,z) += dpdt;}
  void addMomentum(int x, int y, int z, BoutReal dnvdt) {ddt(NVn)(x,y,z) += dnvdt;}
  
private:
  Field3D Nn, Pn, NVn; // Density, pressure and parallel momentum

  Field3D Pnlim;  // Limited pressure, used to calculate pressure-driven diffusive flows
  
  Field3D Dnn;
  
  bool sheath_ydown, sheath_yup;
  
  BoutReal neutral_gamma; // Heat transmission for neutrals
  
  BoutReal numdiff;  // Numerical dissipation

  BoutReal nn_floor; // Minimum Nn used when dividing NVn by Nn to get Vn.
};


#endif // __NEUTRAL_MIXED_H__
