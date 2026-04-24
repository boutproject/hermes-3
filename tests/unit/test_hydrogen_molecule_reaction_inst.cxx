#include "../../include/hydrogen_molecule_reactions.hxx"
#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx"
#include "gtest/gtest.h"

namespace bout {
namespace globals {
extern Mesh* mesh;
}
} // namespace bout
using namespace bout::globals;

using HydrogenMoleculeReactionInstTest = FakeMeshFixture;

static Options base_options{
    {"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}}};


TEST_F(HydrogenMoleculeReactionInstTest, CreateMolHDissociativeExc) {
  hermes::ReactionBase::reset_instance_counter();
  Options options = base_options.copy();
  options["test"]["type"] = "h2+ + e -> h + h+ + e";
  hermes::MolHDissociativeExc<hermes::HIsotope::h> component("test", options, nullptr);
}

TEST_F(HydrogenMoleculeReactionInstTest, CreateMolHDissociativeIzn) {
  hermes::ReactionBase::reset_instance_counter();
  Options options = base_options.copy();
  options["test"]["type"] = "h2 + e -> h + h+ + 2e";
  hermes::MolHDissociativeIzn<hermes::HIsotope::h> component("test", options, nullptr);
}
