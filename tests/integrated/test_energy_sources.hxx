#pragma once

#include "component.hxx"
#include <bout/constants.hxx>

/// Constant SI sources used only by integrated-test builds. The options remain
/// scalar in BOUT.settings, allowing independent audits without volume outputs.
struct TestEnergySources : public NamedComponent<TestEnergySources> {
  TestEnergySources(std::string name, Options& alloptions, Solver*)
      : NamedComponent(name, {readWrite("species:{name}:density_source"),
                              readWrite("species:{name}:momentum_source"),
                              readWrite("species:{name}:energy_source")}),
        species_name(name) {
    auto& options = alloptions[name];
    const auto& units = alloptions["units"];
    const BoutReal density = units["inv_meters_cubed"];
    const BoutReal frequency = 1. / units["seconds"].as<BoutReal>();
    const BoutReal speed = units["meters"].as<BoutReal>() * frequency;
    const BoutReal pressure = SI::qe * density * units["eV"].as<BoutReal>();
    density_source =
        options["test_source_density"].as<BoutReal>() / (density * frequency);
    momentum_source = options["test_source_momentum"].as<BoutReal>()
                      / (SI::Mp * density * speed * frequency);
    thermal_source =
        options["test_source_energy"].as<BoutReal>() / (pressure * frequency);
    substitutePermissions("name", {name});
  }

  static constexpr auto type = "test_energy_sources";

private:
  std::string species_name;
  BoutReal density_source, momentum_source, thermal_source;

  void transform_impl(GuardedOptions& state) override {
    auto species = state["species"][species_name];
    add(species["density_source"], Field3D{density_source});
    add(species["momentum_source"], Field3D{momentum_source});
    add(species["energy_source"], Field3D{thermal_source});
  }
};

namespace {
RegisterComponent<TestEnergySources> register_test_energy_sources;
}
