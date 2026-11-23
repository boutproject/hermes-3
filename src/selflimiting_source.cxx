#include "../include/selflimiting_source.hxx"

#include <bout/constants.hxx>
#include <bout/fv_ops.hxx>
#include <bout/mesh.hxx>
#include <bout/options.hxx>
#include <bout/output_bout_types.hxx>

using bout::globals::mesh;

namespace {
// Selflimiting source function
// This function uses logistic growth to dampen this over the interval overshoot from 1 to 0
// k = 5.0/overshoot is chosen so the function goes from 1 to ~0.1 in the inverval of the overshoot
// This function just returns the correction, which has to be multiplied with the actual source
Field3D sourceterm(const Field3D& thisfield, const BoutReal maximum,
                   const BoutReal overshoot) {
  BoutReal k = 5.0 / overshoot;
  BoutReal x0 = maximum + overshoot / 2.0;
  return (1.0 - (1.0 / (1.0 + exp(-k * (thisfield - x0)))));
}
} // namespace

SelflimitingSource::SelflimitingSource(std::string name, Options& alloptions, Solver*)
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

  auto& soptions = options["selflimiting_source"];

  density_setpoint = soptions["density_setpoint"]
                         .doc("Density [m^-3] above which the source term will be "
                              "dampened. Only active when > 0.")
                         .withDefault<BoutReal>(-1.0)
                     / Nnorm;

  temperature_setpoint = soptions["temperature_setpoint"]
                             .doc("Temperature [eV] above which the source term will be "
                                  "dampened. Only active when > 0.")
                             .withDefault<BoutReal>(-1.0)
                         / Tnorm;

  if (density_setpoint > 0.0) {

    BoutReal densitysource_normalisation = Nnorm * Omega_ci;

    densityform = soptions["density_form"]
                      .doc("Form of the selflimiting density source.")
                      .withDefault<Field3D>(Field3D{0.0})
                  / densitysource_normalisation;

    density_overshoot = soptions["density_overshoot"]
                            .doc("Density [m^-3] interval in which the source term is "
                                 "dampened when density_setpoint > 0.")
                            .withDefault<BoutReal>(2e18)
                        / Nnorm;

    if (density_overshoot <= 0.0) {
      throw BoutException("Density overshoot must be positive but is {}",
                          density_overshoot);
    }
  }

  if (temperature_setpoint > 0.0) {
    BoutReal pressuresource_normalisation =
        SI::qe * Nnorm * Tnorm * Omega_ci; // [Pa/s] or [W/m^3] if converted to energy

    pressureform = soptions["pressure_form"]
                       .doc("Form of the selflimiting pressure source.")
                       .withDefault<Field3D>(Field3D{0.0})
                   / pressuresource_normalisation;

    temperature_overshoot =
        soptions["temperature_overshoot"]
            .doc("Temperature [eV] interval in which the source term is "
                 "dampened when temperature_setpoint > 0.")
            .withDefault<BoutReal>(10.0)
        / Tnorm;

    if (temperature_overshoot <= 0.0) {
      throw BoutException("Temperature overshoot must be positive but is {}",
                          temperature_overshoot);
    }
  }

  substitutePermissions("name", {name});
  std::vector<std::string> output_vars;
  if (density_setpoint > 0.0) {
    output_vars.push_back("density_source");
  }
  if (temperature_setpoint > 0.0) {
    output_vars.push_back("energy_source");
  }

  substitutePermissions("output", output_vars);
}

void SelflimitingSource::transform_impl(GuardedOptions& state) {

  GuardedOptions species = state["species"][objectName()];

  if (species.isSet("density") && (density_setpoint > 0.0)) {
    const Field3D N = GET_NOBOUNDARY(Field3D, species["density"]);
    density_source = densityform * sourceterm(N, density_setpoint, density_overshoot);
    add(species["density_source"], density_source);
  }

  if (species.isSet("temperature") && (temperature_setpoint > 0.0)) {
    const Field3D T = GET_NOBOUNDARY(Field3D, species["temperature"]);
    energy_source = (3.0 / 2.0) * pressureform
                    * sourceterm(T, temperature_setpoint, temperature_overshoot);
    add(species["energy_source"], energy_source);
  }
}

void SelflimitingSource::outputVars(Options& state) {

  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

  const std::string& name = objectName();

  if (diagnose) {
    if (density_source.isAllocated()) {
      set_with_attrs(state[name + std::string("_selflimiting_densitysource")],
                     density_source,
                     {{"time_dimension", "t"},
                      {"units", "m^-3 s^-1"},
                      {"conversion", Nnorm * Omega_ci},
                      {"standard_name", "selflimiting density source"},
                      {"species", name},
                      {"source", "selflimiting_source"}});
    }

    if (energy_source.isAllocated()) {

      set_with_attrs(state[name + std::string("_selflimiting_energysource")],
                     energy_source,
                     {{"time_dimension", "t"},
                      {"units", "W m^-3"},
                      {"conversion", Pnorm * Omega_ci},
                      {"standard_name", "selflimiting energy source"},
                      {"species", name},
                      {"source", "selflimiting_source"}});
    }
  }
}
