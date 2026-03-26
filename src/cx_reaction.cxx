#include "cx_reaction.hxx"
#include "hermes_utils.hxx"
#include "rate_helper.hxx"
#include "species_parser.hxx"

///
CXReaction::CXReaction(std::string name, Options& alloptions, Solver*)
    : CXReaction(name, alloptions) {}

///
CXReaction::CXReaction(std::string name, Options& alloptions)
    : Reaction(name, alloptions) {

  this->do_parallel_averaging = false;

  // Identify CX roles for each species and validate reactants and products
  set_species_and_validate();

  // Multiplier is associated with the lower ionisation state reactant
  rate_multiplier = alloptions[this->r1]["K_cx_multiplier"]
                        .doc("Scale the charge exchange rate by this factor")
                        .withDefault<BoutReal>(1.0);

  if (this->has_neutral_reactant) {
    /* This is useful for testing the impact of enabling the neutral momentum equation.
     * When set to true, neutral-ion CX behaves as if using diffusive neutrals but the
     * neutral transport still enjoys the full momentum equation treatment.
     */
    this->no_neutral_cx_mom_gain =
        alloptions[name]["no_neutral_cx_mom_gain"]
            .doc("If true, ion momentum in neutral-ion CX is still lost, but not "
                 "given to the neutrals")
            .withDefault<bool>(false);
  }

  // Energy exchange between r1 and p1
  set_energy_channel_weight(this->r1, this->p1, 1.0);
  set_energy_channel_weight(this->r1, this->p2, 0.0);

  // Energy exchange between r2 and p2
  set_energy_channel_weight(this->r2, this->p1, 0.0);
  set_energy_channel_weight(this->r2, this->p2, 1.0);

  //  Momentum exchange between r1 and p1
  set_momentum_channel_weight(this->r1, this->p1, 1.0);
  set_momentum_channel_weight(this->r1, this->p2, 0.0);

  // Momentum exchange between r2 and p2
  // (Can be turned off for neutral-ion CX using config option)
  if (this->has_neutral_reactant && this->no_neutral_cx_mom_gain) {
    set_momentum_channel_weight(this->r2, this->p2, 0.0);
  } else {
    set_momentum_channel_weight(this->r2, this->p2, 1.0);
  }
  set_momentum_channel_weight(this->r2, this->p1, 0.0);

  // Set permissions
  setPermissions(readOnly("species:{reactant}:{react_vals}"));
  setPermissions(readOnly("species:{sp}:{read_vals}"));
  setPermissions(readWrite("species:{sp}:{writevals}"));
  setPermissions(readWrite("species:{r1}:collision_frequencies:{r1}_{r2}_cx"));
  setPermissions(readWrite("species:{r2}:collision_frequencies:{r2}_{r1}_cx"));

  std::vector<std::string> writevals = {"momentum_source", "energy_source"};
  if (this->r1 != this->p2) {
    writevals.push_back("density_source");
  }
  substitutePermissions("reactant", {this->r1, this->r2});
  substitutePermissions("react_vals", {"density", "temperature"});
  substitutePermissions("read_vals", {"AA", "velocity"});
  substitutePermissions("sp", {this->r1, this->r2, this->p1, this->p2});
  substitutePermissions("writevals", writevals);
  substitutePermissions("r1", {this->r1});
  substitutePermissions("r2", {this->r2});

  if (diagnose) {
    // Set up diagnostics to be written to dump file
    DiagnosticTransformerType default_transformer = identity;
    std::string r1_long_name = SpeciesParser(this->r1).long_name();
    std::string standard_name = fmt::format("{}_charge_exchange", r1_long_name);

    // Asymmetric CX, density channel, two momentum channels, two energy channels
    if (!this->parser->is_symmetric()) {
      // Particle/density source diagnostic (keyed by reactant that loses an electron)
      add_diagnostic(
          this->r1, fmt::format("S{:s}{:s}_cx", this->r1, this->r2),
          fmt::format("Particle transfer to {:s} from {:s} due to CX with {:s}", this->r1,
                      this->p1, this->r2),
          ReactionDiagnosticType::density_src, standard_name, identity,
          "particle transfer");

      // p1 momentum diagnostic
      add_diagnostic(
          this->p1, fmt::format("F{:s}{:s}_cx", this->r1, this->r2),
          fmt::format("Momentum transfer to {:s} from {:s} due to CX with {:s}", this->r1,
                      this->p1, this->r2),
          ReactionDiagnosticType::momentum_src, standard_name, negate);

      // p1 energy diagnostic
      add_diagnostic(this->p1, fmt::format("E{:s}{:s}_cx", this->r1, this->r2),
                     fmt::format("Energy transfer to {:s} from {:s} due to CX with {:s}",
                                 this->r1, this->p1, this->r2),
                     ReactionDiagnosticType::energy_src, standard_name, negate,
                     "energy transfer");

      // p2 momentum diagnostic
      add_diagnostic(
          this->p2, fmt::format("F{:s}{:s}_cx", this->r2, this->r1),
          fmt::format("Momentum transfer to {:s} from {:s} due to CX with {:s}", this->r2,
                      this->p2, this->r1),
          ReactionDiagnosticType::momentum_src, standard_name, negate,
          "momentum transfer");
      // p2 energy diagnostic
      add_diagnostic(this->p2, fmt::format("E{:s}{:s}_cx", this->r2, this->r1),
                     fmt::format("Energy transfer to {:s} from {:s} due to CX with {:s}",
                                 this->r2, this->p2, this->r1),
                     ReactionDiagnosticType::energy_src, standard_name, negate,
                     "energy transfer");
    } else {

      // Momentum diagnostic (keyed by p1)
      add_diagnostic(
          this->p1, fmt::format("F{:s}{:s}_cx", this->r1, this->r2),
          fmt::format("Momentum transfer to {:s} from {:s} due to CX with {:s}", this->r1,
                      this->p1, this->r2),
          ReactionDiagnosticType::momentum_src, standard_name, negate);

      // Energy diagnostic (keyed by p1)
      add_diagnostic(this->p1, fmt::format("E{:s}{:s}_cx", this->r1, this->r2),
                     fmt::format("Energy transfer to {:s} from {:s} due to CX with {:s}",
                                 this->r1, this->p1, this->r2),
                     ReactionDiagnosticType::energy_src, standard_name, negate,
                     "energy transfer");
    }

    // Collision frequency added for both symmetric and asymmetric CX, keyed by r1
    add_diagnostic(this->r1, fmt::format("K{:s}{:s}_cx", this->r1, this->r2),
                   fmt::format("Collision frequency of CX of {:s} and {:s} producing "
                               "{:s} and {:s}. Note Kab != Kba",
                               this->r1, this->p1, this->r2, this->p2),
                   ReactionDiagnosticType::collision_freq, standard_name, identity);
  }
}

///
void CXReaction::set_species_and_validate() {
  std::vector<std::string> reactants =
      this->parser->get_species(species_filter::reactants);

  // Must have 2 reactants
  if (reactants.size() != 2) {
    throw BoutException(fmt::format(
        "Expected exactly two reactant species in CX reaction, but found {}: {}",
        reactants.size(), fmt::join(reactants, ", ")));
  }

  SpeciesParser firstr = SpeciesParser(reactants[0]);
  SpeciesParser secondr = SpeciesParser(reactants[1]);

  // r1 is the reactant in the lower charge state.
  if (firstr.get_charge() < secondr.get_charge()) {
    this->r1 = reactants[0];
    this->r2 = reactants[1];
  } else {
    this->r1 = reactants[1];
    this->r2 = reactants[0];
  }

  // Record whether we have a neutral reactant (special behaviour for neutral-ion CX)
  this->has_neutral_reactant = (firstr.get_charge() == 0 || secondr.get_charge() == 0);

  this->p1 = SpeciesParser(this->r1).ionised().get_str();
  this->p2 = SpeciesParser(this->r2).recombined().get_str();

  std::vector<std::string> products = this->parser->get_species(species_filter::products);
  // Must have 2 products
  if (products.size() != 2) {
    throw BoutException(fmt::format(
        "Expected exactly two product species in CX reaction, but found {}: {}",
        products.size(), fmt::join(products, ", ")));
  }

  // Verify that the products are as expected for a CX reaction between r1 and r2.
  if (std::find(products.begin(), products.end(), this->p1) == products.end()) {
    throw BoutException(
        fmt::format("Expected species '{}' not found in reaction products", this->p1));
  }

  if (std::find(products.begin(), products.end(), this->p2) == products.end()) {
    throw BoutException(
        fmt::format("Expected species '{}' not found in reaction products", this->p2));
  }
}

///
void CXReaction::transform_additional(GuardedOptions& state, const RateData& rate_data) {

  // Refs to species data for convenience
  GuardedOptions r1s = state["species"][this->r1];
  GuardedOptions r2s = state["species"][this->r2];
  GuardedOptions p1s = state["species"][this->p1];
  GuardedOptions p2s = state["species"][this->p2];

  auto r1_vel = get<Field3D>(r1s["velocity"]);

  // Frictional heating converts kinetic energy to thermal energy
  // Handles both symmetric and asymmetric CX
  const BoutReal r1_AA = get<BoutReal>(r1s["AA"]);
  auto p1_vel = get<Field3D>(p1s["velocity"]);
  add(p1s["energy_source"], 0.5 * r1_AA * rate_data.rate * SQ(p1_vel - r1_vel));

  const BoutReal r2_AA = get<BoutReal>(r2s["AA"]);
  auto r2_vel = get<Field3D>(r2s["velocity"]);
  auto p2_vel = get<Field3D>(p2s["velocity"]);
  add(p2s["energy_source"], 0.5 * r2_AA * rate_data.rate * SQ(p2_vel - r2_vel));

  // Set 'collision_frequencies' values
  std::string r1_coll_freq_key =
      fmt::format("collision_frequencies:{:s}_{:s}_cx", this->r1, this->r2);
  update_source<set>(state, this->r1, ReactionDiagnosticType::collision_freq,
                     r1_coll_freq_key, rate_data.coll_freq(this->r1));

  std::string r2_coll_freq_key =
      fmt::format("collision_frequencies:{:s}_{:s}_cx", this->r2, this->r1);
  update_source<set>(state, this->r2, ReactionDiagnosticType::collision_freq,
                     r2_coll_freq_key, rate_data.coll_freq(this->r2));
}