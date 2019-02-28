/*
 * Full velocity Navier-Stokes model
 */

#ifndef __NEUTRAL_FULL_VELOCITY_H__
#define __NEUTRAL_FULL_VELOCITY_H__

#include "neutral-model.hxx"

class FullVelocity : public NeutralModel {
public:
  FullVelocity(Solver *solver, Mesh *mesh, Options &options);
  ~FullVelocity() {}

  /// Update plasma quantities
  void update(const Field3D &Ne, const Field3D &Te, const Field3D &Ti, const Field3D &Vi);
  
  void addDensity(int x, int y, int z, BoutReal dndt);
  void addPressure(int x, int y, int z, BoutReal dpdt);
  void addMomentum(int x, int y, int z, BoutReal dnvdt);
  
private:
  Field2D Nn2D;         // Neutral gas density (evolving)
  Field2D Pn2D;         // Neutral gas pressure (evolving)
  Vector2D Vn2D;        // Neutral gas velocity
  Field2D Tn2D;
  Field2D DivV2D;       // Divergence of gas velocity

  // Transformation to cylindrical coordinates
  
  // Grad x = Txr * Grad R + Txz * Grad Z
  // Grad y = Tyr * Grad R + Tyz * Grad Z
  Field2D Txr, Txz;
  Field2D Tyr, Tyz;
  
  // Grad R = Urx * Grad x + Ury * Grad y
  // Grad Z = Uzx * Grad x + Uzy * Grad y
  Field2D Urx, Ury;
  Field2D Uzx, Uzy;

  BoutReal gamma_ratio;    // Ratio of specific heats
  BoutReal neutral_viscosity; // Neutral gas viscosity
  BoutReal neutral_bulk;   // Neutral gas bulk viscosity
  BoutReal neutral_conduction; // Neutral gas thermal conduction
  BoutReal neutral_gamma; // Heat transmission for neutrals
  
  // Outflowing boundaries for neutrals
  bool outflow_ydown; // Allow neutral outflows?
};


#endif // __NEUTRAL_MODEL_H__
