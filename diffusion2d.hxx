
#ifndef __NEUTRAL_DIFFUSION2D_H__
#define __NEUTRAL_DIFFUSION2D_H__

#include "neutral-model.hxx"
#include <invert_laplace.hxx>

class Diffusion2D : public NeutralModel {
public:
  Diffusion2D(Solver *solver, Mesh *mesh, Options &options);
  ~Diffusion2D();

  void update(const Field3D &Ne, const Field3D &Te, const Field3D &Ti, const Field3D &Vi);
  
  void precon(BoutReal t, BoutReal gamma, BoutReal delta);
private:
  Field3D Nn;  // Neutral gas density (evolving)
  Field3D Pn;  // Neutral gas pressure (evolving)
  Field3D Dnn; // Neutral gas diffusion rate
  
  BoutReal Lmax; // Maximum mean free path [m]

  Laplacian *inv; // Laplacian inversion used for preconditioning
};

#endif // __NEUTRAL_DIFFUSION2D_H__
