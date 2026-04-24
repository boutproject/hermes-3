#pragma once
#ifndef TEST_IZNREC_REACTIONS_H
#define TEST_IZNREC_REACTIONS_H

#include "test_reactions.hxx"

namespace hermes {

/**
 * @brief Class to test reactions of the form heavy_species1 + <M>e -> heavy_species2 +
 * <N>e where M and N are positive integers.
 *
 */

template <typename RTYPE>
class IznRecReactionTest : public ReactionTest<RTYPE> {

  // The code below can be uncommented once impurity ionisation and recombination
  // reactions inherit from IznRecReaction.

  // static_assert(std::is_base_of<IznRecReaction, RTYPE>(), "Template arg to
  // ReactionTest must derive from IznRecReaction");

protected:
  IznRecReactionTest(std::string lbl, std::string reaction_str)
      : ReactionTest<RTYPE>(lbl, reaction_str) {

    this->heavy_reactant =
        this->parser.get_single_species(species_filter::reactants, species_filter::heavy);
    this->heavy_product =
        this->parser.get_single_species(species_filter::products, species_filter::heavy);
  };

  std::string heavy_reactant;
  std::string heavy_product;

  Options generate_state() final {
    // Assume 1 heavy reactant, 1 heavy product for izn/rec. Bail out otherwise.
    std::size_t num_heavy = this->parser.get_species(species_filter::heavy).size();
    if (num_heavy != 2) {
      throw BoutException(fmt::format(
          "IznRecReactionTest::generate_state assumes exactly one heavy reactant "
          "and one heavy product per reaction; found {:d} heavy species.",
          num_heavy));
    }

    // N.B. No attempt to set the correct masses for heavy species; always set to 1
    std::string comp_name("test" + this->lbl);
    Field3D e_vel(1.0);
    Field3D hp_dens(1.0);
    Field3D hp_temp(1.0);
    const int hr_q = SpeciesParser(this->heavy_reactant).get_charge();
    const int hp_q = SpeciesParser(this->heavy_product).get_charge();
    Options state{{comp_name,
                   {
                       {"type", this->parser.get_reaction_str()},
                       {"diagnose", true},
                   }},
                  {"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"species",
                   {{"e", {{"AA", 1. / 1836}, {"velocity", e_vel}}},
                    {this->heavy_reactant, {{"AA", 1.0}, {"charge", hr_q}}},
                    {this->heavy_product,
                     {{"AA", 1.0},
                      {"charge", hp_q},
                      {"density", hp_dens},
                      {"temperature", hp_temp}}}}}};

    // Linear functions for various fields that are inputs to the reaction transforms
    state["species"]["e"]["density"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(this->logn_min, this->logn_max, linfunc_axis::y), &state,
        mesh);
    state["species"]["e"]["temperature"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(this->logT_min, this->logT_max, linfunc_axis::z), &state,
        mesh);
    state["species"][this->heavy_reactant]["density"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(this->logn_min, this->logn_max, linfunc_axis::x), &state,
        mesh);
    state["species"][this->heavy_reactant]["temperature"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(this->logT_min, this->logT_max, linfunc_axis::y), &state,
        mesh);
    state["species"][this->heavy_reactant]["velocity"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(this->logv_min, this->logv_max, linfunc_axis::x), &state,
        mesh);
    state["species"][this->heavy_product]["velocity"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(this->logv_min, this->logv_max, linfunc_axis::z), &state,
        mesh);
    return state;
  }
};

class HIznTest : public IznRecReactionTest<IznReaction> {
public:
  HIznTest() : IznRecReactionTest<IznReaction>("HIzn", "h + e -> h+ + 2e") {}
};

class DIznTest : public IznRecReactionTest<IznReaction> {
public:
  DIznTest() : IznRecReactionTest<IznReaction>("DIzn", "d + e -> d+ + 2e") {}
};

class TIznTest : public IznRecReactionTest<IznReaction> {
public:
  TIznTest() : IznRecReactionTest<IznReaction>("TIzn", "t + e -> t+ + 2e") {}
};

class HRecTest : public IznRecReactionTest<RecReaction> {
public:
  HRecTest() : IznRecReactionTest<RecReaction>("HRec", "h+ + e -> h") {}
};

class DRecTest : public IznRecReactionTest<RecReaction> {
public:
  DRecTest() : IznRecReactionTest<RecReaction>("DRec", "d+ + e -> d") {}
};

class TRecTest : public IznRecReactionTest<RecReaction> {
public:
  TRecTest() : IznRecReactionTest<RecReaction>("TRec", "t+ + e -> t") {}
};

class HeIzn01Test : public IznRecReactionTest<IznReaction> {
public:
  HeIzn01Test() : IznRecReactionTest<IznReaction>("HeIzn01", "he + e -> he+ + 2e") {}
};

class HeRec10Test : public IznRecReactionTest<RecReaction> {
public:
  HeRec10Test() : IznRecReactionTest<RecReaction>("HeRec10", "he+ + e -> he") {}
};

} // namespace hermes

#endif // TEST_IZNREC_REACTIONS_H
