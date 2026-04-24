#pragma once
#ifndef TEST_MOLECULE_REACTIONS_H
#define TEST_MOLECULE_REACTIONS_H

#include "hydrogen_molecule_reactions.hxx"

namespace hermes {

/**
 * @brief Class to test dissociative reactions.
 *
 */

template <typename RTYPE>
class DissociativeReactionTest : public ReactionTest<RTYPE> {

  static_assert(std::is_base_of<Reaction, RTYPE>(),
                "Template arg to DissociativeReactionTest must derive from Reaction");

protected:
  DissociativeReactionTest(std::string lbl, std::string reaction_str)
      : ReactionTest<RTYPE>(lbl, reaction_str) {

    // Generate sorted lists of heavy reactants and products
    this->heavy_reactants =
        this->parser.get_species(species_filter::reactants, species_filter::heavy);
    std::sort(this->heavy_reactants.begin(), this->heavy_reactants.end());
    this->heavy_products =
        this->parser.get_species(species_filter::products, species_filter::heavy);
    std::sort(this->heavy_products.begin(), this->heavy_products.end());
    this->heavy_species = this->parser.get_species(species_filter::heavy);
    std::sort(this->heavy_species.begin(), this->heavy_species.end());
  };

  /// Sorted lists of heavy reactants, products, all species, for use in generating the
  /// test state.
  std::vector<std::string> heavy_reactants;
  std::vector<std::string> heavy_products;
  std::vector<std::string> heavy_species;

  /**
   * @brief Generated state should work with any number of heavy reactants/products.
   *
   * @return Options
   */
  Options generate_state() final {
    // Must have at least one heavy reactant and one heavy product.
    if (heavy_reactants.size() == 0 || heavy_products.size() == 0) {
      throw BoutException(
          fmt::format("DissociativeReactionTest::generate_state expects reactions with "
                      "at least one heavy reactant and one heavy product."));
    }

    std::string comp_name("test" + this->lbl);
    Field3D e_vel(1.0);
    Options state{{comp_name,
                   {
                       {"type", this->parser.get_reaction_str()},
                       {"diagnose", true},
                   }},
                  {"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"species", {{"e", {{"AA", 1. / 1836}}}}}};

    // Linear functions for various fields that are inputs to the reaction transforms

    // Cycle through axes when setting fields, starting with x
    // (species names are sorted, so the state generated should be order independent)
    linfunc_axis axis = linfunc_axis::x;

    // Densities and temperatures for electrons
    state["species"]["e"]["density"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(this->logn_min, this->logn_max, axis++), &state, mesh);
    state["species"]["e"]["temperature"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(this->logT_min, this->logT_max, axis++), &state, mesh);

    // Densities and temperatures for all heavy reactants
    for (const auto& reactant : this->heavy_reactants) {
      state["species"][reactant]["density"] = FieldFactory::get()->create3D(
          this->gen_lin_field_str(this->logn_min, this->logn_max, axis++), &state, mesh);

      state["species"][reactant]["temperature"] = FieldFactory::get()->create3D(
          this->gen_lin_field_str(this->logT_min, this->logT_max, axis++), &state, mesh);
    }

    // Masses, velocities for all heavy species
    for (const auto& sp_name : this->heavy_species) {
      state["species"][sp_name]["AA"] = 1.0;
      state["species"][sp_name]["velocity"] = FieldFactory::get()->create3D(
          this->gen_lin_field_str(this->logv_min, this->logv_max, axis++), &state, mesh);
    }
    return state;
  }
};

// Dissociation of H2
class H2DissTest : public DissociativeReactionTest<MolHDissociation<HIsotope::h>> {
public:
  H2DissTest()
      : DissociativeReactionTest<MolHDissociation<HIsotope::h>>("H2Diss",
                                                                "h2 + e -> 2h + e") {}
};

// Dissociative excitation of T2
class T2pDissExcTest : public DissociativeReactionTest<MolHDissociativeExc<HIsotope::t>> {
public:
  T2pDissExcTest()
      : DissociativeReactionTest<MolHDissociativeExc<HIsotope::t>>(
            "T2pDissExc", "t2+ + e -> t + t+ + e") {}
};

// Dissociative ionisation of D2
class D2DissIznTest : public DissociativeReactionTest<MolHDissociativeIzn<HIsotope::d>> {
public:
  D2DissIznTest()
      : DissociativeReactionTest<MolHDissociativeIzn<HIsotope::d>>(
            "D2DissIzn", "d2 + e -> d + d+ + 2e") {}
};

// Non-dissociative ionisation of H2
class H2NonDissIznTest
    : public DissociativeReactionTest<MolHNonDissociativeIzn<HIsotope::h>> {
public:
  H2NonDissIznTest()
      : DissociativeReactionTest<MolHNonDissociativeIzn<HIsotope::h>>(
            "H2NonDissIzn", "h2 + e -> h2+ + 2e") {}
};

// Dissociative recombination of T2+
class T2pDissRecTest : public DissociativeReactionTest<MolHDissociativeRec<HIsotope::t>> {
public:
  T2pDissRecTest()
      : DissociativeReactionTest<MolHDissociativeRec<HIsotope::t>>("T2pDissRec",
                                                                   "t2+ + e -> 2t") {}
};

// CX involving molecules
class D2DpCXTest : public CXReactionTest<CXReaction> {
public:
  D2DpCXTest() : CXReactionTest<CXReaction>("D2DpCX", "d2 + d+ -> d2+ + d") {}
};

// Elastic collision between H2 and H+
// Use DissociativeReactionTest fixture for now
class H2HpElasticCollisionTest
    : public DissociativeReactionTest<MolHElasticCollision<HIsotope::h>> {
public:
  H2HpElasticCollisionTest()
      : DissociativeReactionTest<MolHElasticCollision<HIsotope::h>>(
            "H2HpElastic", "h2 + h+ -> h2 + h+") {}
};

} // namespace hermes

#endif // TEST_MOLECULE_REACTIONS_H
