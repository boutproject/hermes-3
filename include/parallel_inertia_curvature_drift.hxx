#pragma once
#ifndef PARALLEL_INERTIA_CURVATURE_DRIFT_H
#define PARALLEL_INERTIA_CURVATURE_DRIFT_H

#include "component.hxx"

/// Calculate parallel inertia curvature flow

/// This term comes from the parallel inertia term in the momentum equation.
/// It is due to parallel inertia in a curved magnetic field (centrifugal-like force). 
/// See ~Vpar^2 term in the paper by in A.V. Chankin, 
/// Journal of Nuclear Materials 241-243 (1997) 199-213 

struct ParallelInertiaCurvatureDrift : public Component {
  ParallelInertiaCurvatureDrift(std::string name, Options &options, Solver *UNUSED(solver));

  /// For every species, if it has:
  ///  - parallel velocity 
  ///  - charge
  ///
  /// Modifies:
  ///  - density_source
  ///  - energy_source
  ///  - momentum_source
  ///  - currents
  void transform(Options &state) override;

  /// Save variables to the output
  void outputVars(Options& state) override;

private:
  Vector2D Curlb_B;
  bool bndry_flux;

  // Diagnostic outputs
  Field3D DivJicurv; // Divergence of parallel inertia curvature current

  bool diagnose; ///< Output additional diagnostics?

};

namespace {
RegisterComponent<ParallelInertiaCurvatureDrift> registercomponentparallelinertiacurvature("parallel_inertia_curvature_drift");
}

#endif // PARALLEL_INERTIA_CURVATURE_DRIFT_H


