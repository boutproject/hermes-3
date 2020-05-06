
#include "../include/isothermal_electrons.hxx"

IsothermalElectrons::IsothermalElectrons(std::string name, Options &alloptions, Solver *UNUSED(solver)) {

  Options& options = alloptions[name];

  auto Tnorm = get<BoutReal>(alloptions["units"]["eV"]);
  Te = options["temperature"]
           .doc("Electron temperature [eV]")
           .withDefault(Tnorm) /
       Tnorm; // Normalise
}

void IsothermalElectrons::transform(Options &state) {

  Options& electrons = state["species"]["e"];
  
  auto Ne = get<Field3D>(electrons["density"]);

  set(electrons["temperature"], Te);
  set(electrons["pressure"], Ne * Te);
}
