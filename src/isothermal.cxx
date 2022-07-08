
#include "../include/isothermal.hxx"

Isothermal::Isothermal(std::string name, Options &alloptions,
                       Solver *UNUSED(solver))
    : name(name) {

  Options& options = alloptions[name];

  auto Tnorm = get<BoutReal>(alloptions["units"]["eV"]);
  T = options["temperature"]
          .doc("Constant temperature [eV]")
          .withDefault<Field3D>(Tnorm) /
      Tnorm; // Normalise

  // Save the temperature to the output files
  bout::globals::dump.addOnce(T, std::string("T") + name);

  if (options["diagnose"]
      .doc("Save additional output diagnostics")
      .withDefault<bool>(false)) {
    // Save pressure as time-varying field
    bout::globals::dump.addRepeat(P, std::string("P") + name);
  }
}

void Isothermal::transform(Options &state) {
  Options& species = state["species"][name];

  // Note: The boundary of N may not be set yet
  auto N = GET_NOBOUNDARY(Field3D, species["density"]);

  set(species["temperature"], T);
  P = N * T;
  set(species["pressure"], P);
}
