#include "../include/adaptive_source.hxx"

#include <bout/constants.hxx>
#include <bout/fv_ops.hxx>
#include <bout/mesh.hxx>
#include <bout/options.hxx>
#include <bout/output_bout_types.hxx>

using bout::globals::mesh;

AdaptiveSource::AdaptiveSource(std::string name, Options& alloptions, Solver*)
    : NamedComponent(name, {readOnly("species:{name}:density", Regions::Interior),
                            readOnly("species:{name}:temperature", Regions::Interior),
                            readWrite("species:{name}:{output}")}) {

  Options& options = alloptions[name];
  // Normalisations
  const Options& units = alloptions["units"];
  const BoutReal Nnorm = units["inv_meters_cubed"];
  const BoutReal Tnorm = units["eV"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();

  diagnose = alloptions[name]["diagnose"]
                 .doc("Output additional diagnostics?")
                 .withDefault<bool>(false);

  BoutReal densitysource_normalisation = Nnorm * Omega_ci;

  adaptive_densityform = options["adaptive_densityform"]
                             .doc("Form of the adaptive density source.")
                             .withDefault<Field3D>(Field3D{0.0})
                         / densitysource_normalisation;

  BoutReal source_normalisation =
      SI::qe * Nnorm * Tnorm * Omega_ci; // [Pa/s] or [W/m^3] if converted to energy

  adaptive_pressureform = options["adaptive_pressureform"]
                              .doc("Form of the adaptive pressure source.")
                              .withDefault<Field3D>(Field3D{0.0});

  adaptive_density = options["adaptive_density"]
                         .doc("Density [m^-3] above which the source term will be "
                              "dampened. Only active when > 0.")
                         .withDefault<BoutReal>(-1.0)
                     / Nnorm;

  adaptive_temperature = options["adaptive_temperature"]
                             .doc("Temperature [eV] above which the source term will be "
                                  "dampened. Only active when > 0.")
                             .withDefault<BoutReal>(-1.0)
                         / Tnorm;

  adaptive_density_overshoot =
      options["adaptive_density_overshoot"]
          .doc("Density [m^-3] interval in which the source term is "
               "dampened when adaptive_density > 0.")
          .withDefault<BoutReal>(2e18)
      / Nnorm;

  adaptive_temperature_overshoot =
      options["adaptive_temperature_overshoot"]
          .doc("Temperature [eV] interval in which the source term is "
               "dampened when adaptive_temperature > 0.")
          .withDefault<BoutReal>(10.0)
      / Tnorm;

  substitutePermissions("name", {name});
  std::vector<std::string> output_vars;
  if (adaptive_density > 0.0) {
    output_vars.push_back("density_source");
  }
  if (adaptive_temperature > 0.0) {
    output_vars.push_back("energy_source");
  }

  substitutePermissions("output", output_vars);
}

void AdaptiveSource::transform_impl(GuardedOptions& state) {

  GuardedOptions species = state["species"][objectName()];

  if (species.isSet("density") && (adaptive_density > 0.0)) {
    const Field3D N = GET_NOBOUNDARY(Field3D, species["density"]);
    const Field3D density_source =
        adaptive_densityform
        * adaptive_sourceterm(N, adaptive_density, adaptive_density_overshoot);
    add(species["density_source"], density_source);
  }

  if (species.isSet("temperature") && (adaptive_temperature > 0.0)) {
    const Field3D T = GET_NOBOUNDARY(Field3D, species["temperature"]);
    const Field3D energy_source =
        (3.0 / 2.0) * adaptive_pressureform
        * adaptive_sourceterm(T, adaptive_temperature, adaptive_temperature_overshoot);
    add(species["energy_source"], energy_source);
  }
}

void AdaptiveSource::outputVars(Options& state) {

  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto rho_s0 = get<BoutReal>(state["rho_s0"]);

  const std::string& name = objectName();

  if (diagnose) {
  }
}
