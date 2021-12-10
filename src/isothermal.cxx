
#include "../include/isothermal.hxx"

Isothermal::Isothermal(std::string name, Options &alloptions,
                       Solver *UNUSED(solver))
    : name(name) {

  Options& options = alloptions[name];

  auto Tnorm = get<BoutReal>(alloptions["units"]["eV"]);
  T = options["temperature"]
          .doc("Constant temperature [eV]")
          .withDefault(Tnorm) /
      Tnorm; // Normalise
}

void Isothermal::transform(Options &state) {

  Options& species = state["species"][name];

  // Note: The boundary of N may not be set yet
  auto N = GET_NOBOUNDARY(Field3D, species["density"]);

  set(species["temperature"], T);
  set(species["pressure"], N * T);
}
