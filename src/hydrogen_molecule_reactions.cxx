#include "hydrogen_molecule_reactions.hxx"

namespace hermes {

// MolHCX implementation
MolHCX::MolHCX(std::string name, Options& options) : CXReaction(name, options) {}

MolHCX::MolHCX(std::string name, Options& options, Solver* solver)
    : CXReaction(name, options, solver) {}

} // namespace hermes
