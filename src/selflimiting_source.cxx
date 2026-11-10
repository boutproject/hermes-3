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
Field3D selflimiting_sourceterm(const Field3D& thisfield, const BoutReal maximum,
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

  BoutReal densitysource_normalisation = Nnorm * Omega_ci;

  auto& soptions = options["selflimiting_source"];

  selflimiting_densityform = soptions["density_form"]
                                 .doc("Form of the selflimiting density source.")
                                 .withDefault<Field3D>(Field3D{0.0})
                             / densitysource_normalisation;

  BoutReal pressuresource_normalisation =
      SI::qe * Nnorm * Tnorm * Omega_ci; // [Pa/s] or [W/m^3] if converted to energy

  selflimiting_pressureform = soptions["pressure_form"]
                                  .doc("Form of the selflimiting pressure source.")
                                  .withDefault<Field3D>(Field3D{0.0})
                              / pressuresource_normalisation;

  selflimiting_density = soptions["density_setpoint"]
                             .doc("Density [m^-3] above which the source term will be "
                                  "dampened. Only active when > 0.")
                             .withDefault<BoutReal>(-1.0)
                         / Nnorm;

  selflimiting_temperature =
      soptions["temperature_setpoint"]
          .doc("Temperature [eV] above which the source term will be "
               "dampened. Only active when > 0.")
          .withDefault<BoutReal>(-1.0)
      / Tnorm;

  selflimiting_density_overshoot =
      soptions["density_overshoot"]
          .doc("Density [m^-3] interval in which the source term is "
               "dampened when selflimiting_density > 0.")
          .withDefault<BoutReal>(2e18)
      / Nnorm;

  selflimiting_temperature_overshoot =
      soptions["temperature_overshoot"]
          .doc("Temperature [eV] interval in which the source term is "
               "dampened when selflimiting_temperature > 0.")
          .withDefault<BoutReal>(10.0)
      / Tnorm;

  substitutePermissions("name", {name});
  std::vector<std::string> output_vars;
  if (selflimiting_density > 0.0) {
    output_vars.push_back("density_source");
  }
  if (selflimiting_temperature > 0.0) {
    output_vars.push_back("energy_source");
  }

  substitutePermissions("output", output_vars);
}

void SelflimitingSource::transform_impl(GuardedOptions& state) {

  GuardedOptions species = state["species"][objectName()];

  if (species.isSet("density") && (selflimiting_density > 0.0)) {
    const Field3D N = GET_NOBOUNDARY(Field3D, species["density"]);
    density_source = selflimiting_densityform
                     * selflimiting_sourceterm(N, selflimiting_density,
                                               selflimiting_density_overshoot);
    add(species["density_source"], density_source);
  }

  if (species.isSet("temperature") && (selflimiting_temperature > 0.0)) {
    const Field3D T = GET_NOBOUNDARY(Field3D, species["temperature"]);
    energy_source = (3.0 / 2.0) * selflimiting_pressureform
                    * selflimiting_sourceterm(T, selflimiting_temperature,
                                              selflimiting_temperature_overshoot);
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
                      {"units", "Pa s^-1"},
                      {"conversion", Pnorm * Omega_ci},
                      {"standard_name", "selflimiting energy source"},
                      {"species", name},
                      {"source", "selflimiting_source"}});
    }
  }
}
