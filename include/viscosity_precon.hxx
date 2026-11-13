#pragma once

// Opt-in parallel viscosity preconditioner.
//
// BraginskiiIonViscosity and BraginskiiElectronViscosity publish their frozen
// coefficient eta here when their `precondition` option is on; EvolveMomentum then
// inverts the linearised viscous force of each published species with InvertParDiv,
// with the sheath-face closure at the targets. Nothing is published or used by default,
// so the default RHS and preconditioner are unchanged.
#include <bout/field3d.hxx>

#include <map>
#include <string>

namespace hermes::viscosity_precon {
struct Frozen {
  Field3D eta; ///< Parallel viscosity coefficient after the flux limiter
};
inline std::map<std::string, Frozen> frozen; ///< Refreshed by every RHS
/// Smooth Bohm condition (sheath_boundary:bohm_smoothing > 0): d(guard velocity) /
/// d(end-cell velocity) in the target guard cells (field-aligned), per ion species.
/// Refreshed by every RHS; nothing is published for the sharp condition
inline std::map<std::string, Field3D> sheath_closure;
} // namespace hermes::viscosity_precon
