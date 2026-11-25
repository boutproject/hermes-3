#pragma once

#include "component.hxx"
#include "div_ops.hxx"
#include "domain_balance.hxx"

#include <map>
#include <string>

/// Whether a scheduled adaptive diffusion component has nonzero rate ceilings.
bool usesAdaptiveIndexDiffusion(const Options& options);

/// Bounded numerical mixing in perpendicular mesh-index space. The component
/// publishes common face conductances and solver-moment increments; equation
/// owners consume them in finally(), outside ordinary thermal source conversion.
struct AdaptiveIndexDiffusion : public NamedComponent<AdaptiveIndexDiffusion> {
  AdaptiveIndexDiffusion(std::string name, Options& options, Solver*);
  void outputVars(Options& state) override;
  static constexpr auto type = "adaptive_index_diffusion";

private:
  void transform_impl(GuardedOptions& state) override;
  bool enabled, density_trigger, temperature_trigger, diagnose;
  BoutReal density_threshold{0.0}, temperature_threshold{0.0};
  BoutReal frequency, norm_density, norm_temperature;
  std::string temperature_species;
  Field2D rate_x, rate_z;
  Field3D activation, loss_rate;
  std::map<std::string, Field3D> density_rates, momentum_rates, kinetic_rates;
  DomainBalance balance;
};

namespace {
RegisterComponent<AdaptiveIndexDiffusion> registercomponentadaptiveindexdiffusion;
}
