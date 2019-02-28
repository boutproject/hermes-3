/*
 * Pseudo-recycling with Nn exponential approximation
 */

#ifndef __NEUTRAL_RECYCLING_H__
#define __NEUTRAL_RECYCLING_H__

#include "neutral-model.hxx"

class NeutralRecycling : public NeutralModel {
public:
  NeutralRecycling(Solver *solver, Mesh *mesh, Options &options);
  ~NeutralRecycling() {}

  /// Update plasma quantities
  void update(const Field3D &Ne, const Field3D &Te, const Field3D &Ti,
              const Field3D &Vi);

private:
  Field3D Nn, Tn;

  Field3D Nn0;                // scaling for neutral density approximation
  Field3D lambda_int, lambda; // for mean free path calculation
  BoutReal Lmax;              // Maximum mean free path [m]
  BoutReal frecycle;          // Recycling fraction

  Field2D hthe;

  // Utility functions
  const Field2D CumSumY2D(const Field2D &f, bool reverse);
  const Field3D CumSumY3D(const Field3D &f, bool reverse);
  const BoutReal bcast_lasty(const BoutReal f);
};

#endif // __NEUTRAL_MODEL_H__
