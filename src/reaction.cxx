#include <memory>
#include <regex>

#include "integrate.hxx"

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

  // Awful hack to extract the correct reaction expression from the params; depends on
  // instantiation order matching input file. There must be a better way...
  std::string reaction_grp_str = alloptions[name]["type"];
  std::regex match_parentheses("\\(|\\)");
  reaction_grp_str = std::regex_replace(reaction_grp_str, match_parentheses, "");
  std::string reaction_str;
  std::stringstream ss(reaction_grp_str);
  int inst_num = get_instance_num() + 1;
  for (auto ii = 0; ii < inst_num; ii++) {
    std::getline(ss, reaction_str, ',');
  }

  this->parser = std::make_unique<ReactionParser>(reaction_str);
}

void Reaction::transform(Options& state) {

  Field3D momentum_exchange, energy_exchange, energy_loss;

  std::vector<std::string> reactant_species = parser->get_reactant_species();

  // Restrict to 2 reactants for now;
  ASSERT1(reactant_species.size() == 2);
  Options& r1 = state["species"][reactant_species[0]];
  Options& r2 = state["species"][reactant_species[1]];

  Field3D n_r1 = get<Field3D>(r1["density"]);
  Field3D n_r2 = get<Field3D>(r2["density"]);
  Field3D Te = get<Field3D>(state["species"]["e"]["temperature"]);

  Field3D reaction_rate = cellAverage(
      [&](BoutReal ne, BoutReal n1, BoutReal te) {
        return ne * n1 * eval_rate_coeff(te * Tnorm, ne * Nnorm) * Nnorm / FreqNorm
               * rate_multiplier;
      },
      n_r1.getRegion("RGN_NOBNDRY"))(n_r1, n_r2, Te);

  // Use the Stoichiometry 'matrix' (vector) to set density sources for each species
  for (auto el : parser->get_stoich()) {
    std::string sp_name = el.first;
    int weight = el.second;
    if (weight != 0) {
      // Density sources
      add(state["species"][sp_name]["density_source"], weight * reaction_rate);
    }
  }

  if (this->diagnose) {
    set_diagnostic_fields(reaction_rate, momentum_exchange, energy_exchange, energy_loss);
  }
}