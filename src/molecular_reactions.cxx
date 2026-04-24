#include "molecular_reactions.hxx"

namespace hermes {

// Dissociation implementation
Dissociation::Dissociation(std::string name, Options& options)
    : Reaction(name, options) {}

Dissociation::Dissociation(std::string name, Options& options, Solver*)
    : Reaction(name, options) {}

// DissociativeExc implementation
DissociativeExc::DissociativeExc(std::string name, Options& options)
    : Reaction(name, options) {}

DissociativeExc::DissociativeExc(std::string name, Options& options, Solver*)
    : Reaction(name, options) {}

// DissociativeIzn implementation
DissociativeIzn::DissociativeIzn(std::string name, Options& options)
    : Reaction(name, options) {}

DissociativeIzn::DissociativeIzn(std::string name, Options& options, Solver*)
    : Reaction(name, options) {}

// NonDissociativeIzn implementation
NonDissociativeIzn::NonDissociativeIzn(std::string name, Options& options)
    : Reaction(name, options) {}

NonDissociativeIzn::NonDissociativeIzn(std::string name, Options& options, Solver*)
    : Reaction(name, options) {}

// DissociativeRec implementation
DissociativeRec::DissociativeRec(std::string name, Options& options)
    : Reaction(name, options) {}

DissociativeRec::DissociativeRec(std::string name, Options& options, Solver*)
    : Reaction(name, options) {}

} // namespace hermes
