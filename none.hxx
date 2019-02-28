/*
 * Null neutral model which does nothing
 */

#ifndef __NEUTRAL_NONE_H__
#define __NEUTRAL_NONE_H__

#include "neutral-model.hxx"

class NeutralNone: public NeutralModel {
public:
  NeutralNone(Solver *solver, Mesh *mesh, Options &options) : NeutralModel(options) {}
  ~NeutralNone() {}

  /// Update plasma quantities
  void update(const Field3D &Ne, const Field3D &Te, const Field3D &Ti, const Field3D &Vi) {}
};


#endif // __NEUTRAL_NONE_H__
