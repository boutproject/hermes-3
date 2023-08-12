#pragma once
#ifndef CLASSICAL_DIFFUSION_H
#define CLASSICAL_DIFFUSION_H

#include "component.hxx"

struct ClassicalDiffusion : public Component {
  ClassicalDiffusion(std::string name, Options& alloptions, Solver*);

  void transform(Options &state) override;

  void outputVars(Options &state) override;
private:
  Field2D Bsq; // Magnetic field squared

  bool diagnose; ///< Output additional diagnostics?
  Field3D Dn; ///< Particle diffusion coefficient
};

namespace {
RegisterComponent<ClassicalDiffusion> registercomponentclassicaldiffusion("classical_diffusion");
}

#endif // CLASSICAL_DIFFUSION_H
