#include "molecular_reactions.hxx"

namespace hermes {

// Dissociation implementation
Dissociation::Dissociation(std::string name, Options& options) : Reaction(name, options) {
  // Copy collision frequency to state with suffix "_diss"
  const std::string dissociated_species =
      this->parser->get_single_species(species_filter::reactants, species_filter::heavy);
  std::string collfreq_lbl = fmt::format("{:s}_diss", dissociated_species);
  add_coll_freq(dissociated_species, collfreq_lbl, dissociated_species);
}

// DissociativeExc implementation
DissociativeExc::DissociativeExc(std::string name, Options& options)
    : Reaction(name, options) {
  // Copy collision frequency to state with suffix "_dissexc"
  const std::string dissociated_species =
      this->parser->get_single_species(species_filter::reactants, species_filter::heavy);
  std::string collfreq_lbl = fmt::format("{:s}_dissexc", dissociated_species);
  add_coll_freq(dissociated_species, collfreq_lbl, dissociated_species);
}

// DissociativeIzn implementation
DissociativeIzn::DissociativeIzn(std::string name, Options& options)
    : Reaction(name, options) {
  // Copy collision frequency to state with suffix "_dissizn"
  const std::string dissociated_species =
      this->parser->get_single_species(species_filter::reactants, species_filter::heavy);
  std::string collfreq_lbl = fmt::format("{:s}_dissizn", dissociated_species);
  add_coll_freq(dissociated_species, collfreq_lbl, dissociated_species);
}

// NonDissociativeIzn implementation
NonDissociativeIzn::NonDissociativeIzn(std::string name, Options& options)
    : Reaction(name, options) {
  // Copy collision frequency to state with suffix "_ndissizn"
  const std::string dissociated_species =
      this->parser->get_single_species(species_filter::reactants, species_filter::heavy);
  std::string collfreq_lbl = fmt::format("{:s}_ndissizn", dissociated_species);
  add_coll_freq(dissociated_species, collfreq_lbl, dissociated_species);
}

// DissociativeRec implementation
DissociativeRec::DissociativeRec(std::string name, Options& options)
    : Reaction(name, options) {
  // Copy collision frequency to state with suffix "_dissrec"
  const std::string dissociated_species =
      this->parser->get_single_species(species_filter::reactants, species_filter::heavy);
  std::string collfreq_lbl = fmt::format("{:s}_dissrec", dissociated_species);
  add_coll_freq(dissociated_species, collfreq_lbl, dissociated_species);
}

} // namespace hermes
