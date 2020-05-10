
#include "../include/fixed_fraction_ions.hxx"

FixedFractionIons::FixedFractionIons(std::string name, Options &alloptions,
                                     Solver *UNUSED(solver)) {

  std::string fractions_str =
      alloptions[name]["fractions"]
          .doc("Comma-separated list of pairs e.g. "
               "'species1@fraction1, species2@fraction2'")
          .as<std::string>();

  for (const auto &pair : strsplit(fractions_str, ',')) {
    auto species_frac = strsplit(pair, '@');
    if (species_frac.size() != 2) {
      throw BoutException("Expected 'species @ fraction', but got %s",
                          pair.c_str());
    }
    std::string species = trim(species_frac.front());
    BoutReal fraction = stringToReal(trim(species_frac.back()));

    fractions.emplace_back(std::pair<std::string, BoutReal>(species, fraction));
  }
  // Check that there are some species
  if (fractions.size() == 0) {
    throw BoutException("No ion species specified. Got fractions = '%s'",
                        fractions_str.c_str());
  }
}

void FixedFractionIons::transform(Options &state) {
  AUTO_TRACE();

  // Electron density
  auto Ne = get<Field3D>(state["species"]["e"]["density"]);

  for (const auto &spec_frac : fractions) {
    set(state["species"][spec_frac.first]["density"],
        Ne * spec_frac.second);
  }
}
