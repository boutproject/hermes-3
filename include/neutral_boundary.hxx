#pragma once
#ifndef NEUTRAL_BOUNDARY_H
#define NEUTRAL_BOUNDARY_H

#include "component.hxx"

namespace {
Ind3D indexAt(const Field3D& f, int x, int y, int z) {
  int ny = f.getNy();
  int nz = f.getNz();
  return Ind3D{(x * ny + y) * nz + z, ny, nz};
}
}; // namespace

/// Per-species boundary condition for neutral particles at
/// sheath (Y) boundaries.
///
/// Sets boundary conditions:
/// - Free boundary conditions on logarithm
///   of density, temperature and pressure
/// - No-flow boundary conditions on velocity
///   and momentum.
///
/// Adds an energy sink corresponding to a flux of heat to the walls.
///
/// Heat flux into the wall is
///   q = gamma_heat * n * T * v_th
///
/// where v_th = sqrt(eT/m) is the thermal speed
///
struct NeutralBoundary : public Component {
  NeutralBoundary(std::string name, Options& options, Solver*);

  ///
  /// state
  ///  - species
  ///    - <name>
  ///      - density       Free boundary
  ///      - temperature   Free boundary
  ///      - pressure      Free boundary
  ///      - velocity [if set] Zero boundary
  ///      - momentum [if set] Zero boundary
  ///      - energy_source  Adds wall losses
  ///
  void transform(Options& state) override;

private:
  std::string name; ///< Short name of species e.g "d"

  BoutReal gamma_heat; ///< Heat flux coefficient

  bool lower_y; ///< Boundary condition at lower y?
  bool upper_y; ///< Boundary condition at upper y?
};

namespace {
RegisterComponent<NeutralBoundary>
    registercomponentneutralboundary("neutral_boundary");
}

#endif // NEUTRAL_BOUNDARY_H
