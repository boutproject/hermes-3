#include "reaction.hxx"

Reaction::Reaction(std::string name, Options& alloptions) : name(name) {

  // Extract common units to member vars
  const auto& units = alloptions["units"];
  Tnorm = get<BoutReal>(units["eV"]);
  Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
  FreqNorm = 1. / get<BoutReal>(units["seconds"]);

  // Define and extract common options 
  diagnose = alloptions[name]["diagnose"]
                 .doc("Output additional diagnostics?")
                 .withDefault<bool>(false);
}

