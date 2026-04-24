#include "hydrogen_molecule_reactions.hxx"

namespace hermes {

// MolHCX implementation
MolHCX::MolHCX(std::string name, Options& options) : CXReaction(name, options) {}

MolHCX::MolHCX(std::string name, Options& options, [[maybe_unused]] Solver* solver)
    : MolHCX(name, options) {}

} // namespace hermes
