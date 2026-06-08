
#include <bout/constants.hxx>

#include "../include/isothermal.hxx"

Isothermal::Isothermal(std::string name, Options& alloptions, Solver* UNUSED(solver))
    : Component({readIfSet(fmt::format("species:{}:density", name), Regions::Interior),
                 readWrite(fmt::format("species:{}:temperature", name)),
                 // FIXME: This is only written if density is set
                 readWrite(fmt::format("species:{}:pressure", name))}),
      name(name) {
  Options& options = alloptions[name];

  auto Tnorm = get<BoutReal>(alloptions["units"]["eV"]);
  T = options["temperature"].doc("Constant temperature [eV]").as<BoutReal>()
      / Tnorm; // Normalise
  is3D = options["is3D"].doc("Is the temperature 3D field?").withDefault<bool>(true);

  diagnose = options["diagnose"]
    .doc("Save additional output diagnostics")
    .withDefault<bool>(false);
}

void Isothermal::transform_impl(GuardedOptions& state) {

  GuardedOptions species = state["species"][name];
  if (is3D) {
    // If the temperature is 3D, we need to set it at every point in the domain
    Field3D Tfield = T;
    set(species["temperature"], Tfield);
  } else {
    // If the temperature is not 3D, we can just set it as a scalar and it will be broadcast to the whole domain
    set(species["temperature"], T);
  }

  // If density is set, also set pressure
  if (isSetFinalNoBoundary(species["density"])) {
    // Note: The boundary of N may not be set yet
    auto N = GET_NOBOUNDARY(Field3D, species["density"]);
    P = N * T;
    set(species["pressure"], P);
  }
}

void Isothermal::outputVars(Options& state) {
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Nnorm = get<BoutReal>(state["Nnorm"]);

  // Save the temperature to the output files
  set_with_attrs(state[std::string("T") + name], T,
                 {{"units", "eV"},
                  {"conversion", Tnorm},
                  {"long_name", name + " temperature"},
                  {"standard_name", "temperature"},
                  {"species", name},
                  {"source", "isothermal"}});

  if (diagnose) {
    // Save pressure as time-varying field
    set_with_attrs(state[std::string("P") + name], P,
                   {{"time_dimension", "t"},
                    {"units", "Pa"},
                    {"conversion", SI::qe * Tnorm * Nnorm},
                    {"long_name", name + " pressure"},
                    {"standard_name", "pressure"},
                    {"species", name},
                    {"source", "isothermal"}});
   }
}
