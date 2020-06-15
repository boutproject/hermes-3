#pragma once
#ifndef SHEATH_BOUNDARY_H
#define SHEATH_BOUNDARY_H

#include "component.hxx"

/// Boundary condition at the wall in Y
///
/// These are based on
/// "Boundary conditions for the multi-ion magnetized plasma-wall transition"
///  by D.Tskhakaya, S.Kuhn. JNM 337-339 (2005), 405-409
///
/// Note: The approximation used here is for ions having similar
/// gyro-orbit sizes
struct SheathBoundary : public Component {
  SheathBoundary(std::string name, Options &options, Solver *);

  void transform(Options &state) override;
private:
  BoutReal Ge; // Secondary electron emission coefficient
  BoutReal sin_alpha; // sin of angle between magnetic field and wall.
  
  bool lower_y; // Boundary on lower y?
  bool upper_y; // Boundary on upper y?
};

namespace {
RegisterComponent<SheathBoundary>
    registercomponentsheathboundary("sheath_boundary");
}

#endif // SHEATH_BOUNDARY_H
