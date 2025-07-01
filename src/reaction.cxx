#include "integrate.hxx"
#include <iomanip>
#include <memory>
#include <regex>

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

/**
 * @brief Use the stoichiometry matrix to compute density, momentum and energy sources.
 * ASSUMES EXACTLY 2 REACTANTS AND 2 PRODUCTS for now.
 *
 * @param state
 */
void Reaction::transform(Options& state) {

  Field3D momentum_exchange, energy_exchange, energy_loss;

  std::vector<std::string> reactant_species = parser->get_reactant_species();

  // Restrict to 2 reactants for now;
  ASSERT1(reactant_species.size() == 2);
  Options& r1 = state["species"][reactant_species[0]];
  Options& r2 = state["species"][reactant_species[1]];

  // Electron (may be the same as r1 or r2)
  Options& electron = state["species"]["e"];

  Field3D n_r1 = get<Field3D>(r1["density"]);
  Field3D n_r2 = get<Field3D>(r2["density"]);

  Field3D n_e = get<Field3D>(electron["density"]);
  Field3D T_e = get<Field3D>(electron["temperature"]);

  // Calculate reaction rate using cell averaging. Optionally scale by multiplier
  Field3D reaction_rate = cellAverage(
      [&](BoutReal n1, BoutReal n2, BoutReal ne, BoutReal te) {
        return n1 * n2 * eval_reaction_rate(te * Tnorm, ne * Nnorm) * Nnorm / FreqNorm
               * rate_multiplier;
      },
      n_e.getRegion("RGN_NOBNDRY"))(n_r1, n_r2, n_e, T_e);

  // Use the Stoichiometry 'matrix' (vector) to set density sources for each species
  auto stoich = parser->get_stoich();
  for (auto el : stoich) {
    std::string sp_name = el.first;
    int pop_change = el.second;
    if (pop_change != 0) {
      // Density sources
      add(state["species"][sp_name]["density_source"], pop_change * reaction_rate);
    }
  }

  // Get heavy reactant species name(s)
  std::vector<std::string> heavy_reactant_species;
  std::copy_if(reactant_species.begin(), reactant_species.end(),
               std::back_inserter(heavy_reactant_species),
               [](std::string sp_name) { return sp_name != "e"; });

  // Get the mass and velocity of the heavy reactant
  Options& rh = state["species"][heavy_reactant_species[0]];
  Field3D v_rh = get<Field3D>(rh["velocity"]);
  BoutReal AA_rh = get<BoutReal>(rh["AA"]);
  Field3D T_rh = get<Field3D>(rh["temperature"]);

  // Get heavy product species name(s)
  std::vector<std::string> product_species = parser->get_product_species();
  std::vector<std::string> heavy_product_species;
  std::copy_if(product_species.begin(), product_species.end(),
               std::back_inserter(heavy_product_species),
               [](std::string sp_name) { return sp_name != "e"; });
  // Get the velocity of the heavy product
  Options& ph = state["species"][heavy_product_species[0]];
  Field3D v_ph = get<Field3D>(rh["velocity"]);

  // Momentum
  momentum_exchange = reaction_rate * AA_rh * v_rh;

  // Momentum sources
  for (auto el : stoich) {
    std::string sp_name = el.first;
    // Skip electron momentum
    if (sp_name.compare("e") == 0) {
      continue;
    }
    int pop_change = el.second;
    auto Gs = state["species"][sp_name]["AA"].as<BoutReal>()
              * state["species"][sp_name]["velocity"].as<Field3D>();
    if (pop_change < 0) {
      momentum_exchange = pop_change * reaction_rate * Gs;
    } else if (pop_change > 0) {
      momentum_exchange = 0;
      for (auto& rsp_name : heavy_reactant_species) {
        int rpop_change = stoich[rsp_name];
        if (rpop_change < 0) {
          auto Gr = state["species"][rsp_name]["AA"].as<BoutReal>()
                    * state["species"][rsp_name]["velocity"].as<Field3D>();
          momentum_exchange = reaction_rate * Gr;
        }
      }
    } else {
      continue;
    }
    add(state["species"][sp_name]["momentum_source"], momentum_exchange);
  }

  // Energy
  add(ph["energy_source"], 0.5 * AA_rh * reaction_rate * SQ(v_rh - v_ph));

  // Ion thermal energy transfer
  energy_exchange = reaction_rate * (3. / 2) * T_rh;
  subtract(rh["energy_source"], energy_exchange);
  add(ph["energy_source"], energy_exchange);

  // Subclasses perform any additional transform tasks
  transform_additional(state, reaction_rate, momentum_exchange, energy_exchange,
                       energy_loss);
}