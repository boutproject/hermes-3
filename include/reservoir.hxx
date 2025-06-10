#pragma once
#ifndef RESERVOIR_H
#define RESERVOIR_H

#include "component.hxx"

/// Simple 1D reservoir model
///
/// This is a simple reservoir model. It's meant for 1D but will work in 2D and 3D as
/// well. Define the reservoir location by using an analytical expression in the input
/// file e.g. reservoir_location = H(y - y_xpt) where y_xpt is the x-point location would
/// create one between the target and the X-point. The reservoir is set to a constant
/// density and your chosen species will equilibrate with it based on the density
/// difference and a provided timescale. The intended use is to provide some form of
/// cross-field transport for neutrals in 1D.
struct Reservoir : public Component {

  /// # Inputs
  /// - <name>
  ///   - reservoir_density     BoutReal, Constant density in m^-3
  ///   - reservoir_location    Field3D, >0 indicates reservoir location
  ///   - reservoir_timescale   BoutReal, timescale in seconds
  ///   - diagnose              bool. Save sources/sinks
  Reservoir(std::string name, Options& alloptions, Solver*);

  /// # Inputs
  /// - species
  ///   - <name>
  ///     - density
  ///     - pressure
  ///     - momentum
  ///
  /// # Outputs
  /// - species
  ///   - <name>
  ///     - density_source
  ///     - energy_source
  ///     - momentum_source
  void transform(Options& state) override;
  void outputVars(Options& state) override;

private:
  std::string name;           ///< Short name of the species e.g. h+
  Field3D reservoir_location; ///< Indicates reservoir if >0
  BoutReal density_div_sol, density_div_pfr, density_main_sol;
  BoutReal velocity_factor_div_sol, velocity_factor_div_pfr, velocity_factor_main_sol;
  Field3D location_div_sol, location_div_pfr, location_main_sol;

  bool diagnose, reservoir_sink_only;
  BoutReal baffle_position, xpoint_position;   // Parallel positions of reservoirs and xpoint

  /// Cell indices where reservoir_location > 0
  Region<Ind3D> region_div_sol;

  Field3D density_source, energy_source, momentum_source;
  Field3D density_source_main_sol, energy_source_main_sol, momentum_source_main_sol;
  Field3D density_source_div_sol, energy_source_div_sol, momentum_source_div_sol;
  Field3D density_source_div_pfr, energy_source_div_pfr, momentum_source_div_pfr;

  Field3D lpar;   // Parallel connection length, 0 at midplane edge
  Field3D area;   // Cross-sectional area in direction of cross-field transport
};

namespace {
RegisterComponent<Reservoir> registercomponentreservoir("reservoir");
}

#endif // RESERVOIR_H
