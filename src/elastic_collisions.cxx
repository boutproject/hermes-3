#include "elastic_collisions.hxx"

namespace hermes {

///
ElasticCollision::ElasticCollision(std::string name, Options& options)
    : Reaction(name, options, false) {
  std::vector<std::string> reactants =
      this->parser->get_species(species_filter::reactants);
  std::vector<std::string> products = this->parser->get_species(species_filter::products);

  // Must have 2 reactants
  if (reactants.size() != 2) {
    throw BoutException(fmt::format("Expected exactly two reactant species in "
                                    "elastic collision reaction, but found {}: {}",
                                    reactants.size(), fmt::join(reactants, ", ")));
  }

  // Check that reactants and products match (order-independent)
  std::sort(reactants.begin(), reactants.end());
  std::sort(products.begin(), products.end());
  ASSERT2(reactants == products);
  if (reactants != products) {
    throw BoutException(fmt::format("Reactant and product species don't match in "
                                    "elastic collision reaction: reactants are {} but "
                                    "products are {}",
                                    fmt::join(reactants, ", "),
                                    fmt::join(products, ", ")));
  }

  // Store reactant names for convenience
  this->r1 = reactants[0];
  this->r2 = reactants[1];

  // Diagnostics
  add_diagnostic(
      this->r1, fmt::format("F{}{}", this->r1, this->r2),
      fmt::format("Momentum exchange due to elastic collisions between {} and {}",
                  this->r1, this->r2),
      ReactionDiagnosticType::momentum_src, this->rate_data->src_str());
}

///
void ElasticCollision::transform_additional(
    [[maybe_unused]] GuardedOptions& state,
    [[maybe_unused]] const RateData& rate_calc_results) {

  // Frictional force between the two species
  const BoutReal& m1 = state["species"][this->r1]["AA"].GetRef<BoutReal>();
  const BoutReal& m2 = state["species"][this->r2]["AA"].GetRef<BoutReal>();
  const Field3D& v1 = state["species"][this->r1]["velocity"].GetRef<Field3D>();
  const Field3D& v2 = state["species"][this->r2]["velocity"].GetRef<Field3D>();

  Field3D mu = (m1 * m2) / (m1 + m2);
  Field3D momentum_exchange = rate_calc_results.rate * mu * (v1 - v2);

  // Subtract from r1 and diagnostic
  update_source<subtract<Field3D>>(state, this->r1, ReactionDiagnosticType::momentum_src,
                                   momentum_exchange);
  // Add to r2
  add(state["species"][this->r2]["momentum_source"], momentum_exchange);
}

} // namespace hermes
