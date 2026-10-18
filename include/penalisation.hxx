#pragma once
#ifndef PENALISATION
#define PENALISATION

#include "component.hxx"

/// Use penalisation to set the perpendicular boundaries
/// Allows to create complex boundaries while only using slab grids
/// The idea is to change ddt(A) so the variable A
/// is evolved towards a certain boundary value A_b
/// ddt(A) = (1 - P) * ddt(A) + (P / tau) * (A_b - A)
/// P being a smooth function ranging from 0 to 1
/// 0 means no boundary, 1 means a point completely immersed
/// The choice of what A_b depends if neumann or dirichlet
/// boundaries are desired.
/// A_b can also enforce physical boundaries like sheath velocities

struct Penalisation : public NamedComponent<Penalisation> {
  /// Inputs
  /// ------
  ///
  /// @param name     Short name for species e.g. "e"
  /// @param alloptions  Component configuration options
  ///   - <name>
  ///     - charge   Required to have a particle charge
  ///     - AA       Atomic mass
  ///

  Penalisation(const std::string& name, Options& alloptions, Solver* UNUSED(solver));

  static constexpr auto type = "penalisation";

private:
  Field3D penalisation_wall;
  Field3D tau;

  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<Penalisation> registercomponentpenalisation;
}

namespace penalisation {

const Field3D neumann(const Field3D& f, const Field3D& ddt_f, const Field3D& pen,
                      const BoutReal pen_time, const BoutReal f_floor);

}

#endif // Penalisation
