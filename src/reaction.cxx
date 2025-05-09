#include "reaction.hxx"

#include <memory>
#include <regex>

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