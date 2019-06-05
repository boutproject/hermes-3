/*
 * Null neutral model which does nothing
 */

#ifndef __NEUTRAL_NONE_H__
#define __NEUTRAL_NONE_H__

#include "neutral-model.hxx"

class NeutralNone: public NeutralModel {
public:
  NeutralNone(Solver *, Mesh *, Options &options) : NeutralModel(options) {}
  ~NeutralNone() {}

  /// Update plasma quantities
  void update(const Field3D &, const Field3D &, const Field3D &,
              const Field3D &) {}
};


#endif // __NEUTRAL_NONE_H__
