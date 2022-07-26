
#include "../include/isothermal.hxx"

Isothermal::Isothermal(std::string name, Options &alloptions,
                       Solver *UNUSED(solver))
    : name(name) {
  AUTO_TRACE();
  Options& options = alloptions[name];

  auto Tnorm = get<BoutReal>(alloptions["units"]["eV"]);
  T = options["temperature"].doc("Constant temperature [eV]").as<BoutReal>()
      / Tnorm; // Normalise

  // Save the temperature to the output files
  bout::globals::dump.addOnce(T, std::string("T") + name);

  if (options["diagnose"]
      .doc("Save additional output diagnostics")
      .withDefault<bool>(false)) {
    // Save pressure as time-varying field
    bout::globals::dump.addRepeat(P, std::string("P") + name);
    P = 0.0;
  }
}

void Isothermal::transform(Options &state) {
  AUTO_TRACE();

  Options& species = state["species"][name];

  set(species["temperature"], T);

  // If density is set, also set pressure
  if (isSetFinalNoBoundary(species["density"])) {
    // Note: The boundary of N may not be set yet
    auto N = GET_NOBOUNDARY(Field3D, species["density"]);
    P = N * T;
    set(species["pressure"], P);
  }
}
