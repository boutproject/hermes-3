#pragma once
#ifndef TEST_CX_REACTIONS_H
#define TEST_CX_REACTIONS_H

#include "species_parser.hxx"
#include "test_reactions.hxx"

namespace hermes {

/**
 * @brief Class to test charge exchange reaction transforms.
 *
 */
template <typename RTYPE = CXReaction>
class CXReactionTest : public ReactionTest<RTYPE> {

protected:
  /**
   * @brief Set up a CX reaction test.
   *
   * @details Note that large parts of this constructor reproduces the logic in
   * CXReaction::set_species_and_validate(). Should be possible to remove this once all
   * charge exchange reactions inherit from CXReaction.
   *
   * @warning The order of terms in \p reaction_str matters because it affects the test
   * input state. So even though the CXReaction class recognises (h + d+ -> h+ + d) and
   * (d+ + h -> h+ + d) as the same reaction, the test result will change.
   *
   * @param lbl
   * @param reaction_str
   */
  CXReactionTest(const std::string& lbl, const std::string& reaction_str)
      : ReactionTest<RTYPE>(lbl, reaction_str) {

    // Confirm that the reaction has exactly 2 reactants
    ASSERT0(this->parser.get_species(species_filter::reactants).size() == 2);

    /*
    The convention used for CX tests is:
      - r1 = first reactant in reaction string
      - r2 = second reactant in reaction string
      - p1 = the product that lost an electron
      - p2 = the product that gained an electron
    */
    this->r1 = this->parser.get_reactant_by_position(1);
    this->r2 = this->parser.get_reactant_by_position(2);
    SpeciesParser r1_species(this->r1);
    SpeciesParser r2_species(this->r2);
    if (r1_species.get_charge() < r2_species.get_charge()) {
      this->p1 = r1_species.ionised().get_str();
      this->p2 = r2_species.recombined().get_str();
    } else {
      this->p1 = r2_species.ionised().get_str();
      this->p2 = r1_species.recombined().get_str();
    }

    // Confirm that the reaction has exactly 2 products, p1 and p2
    std::vector<std::string> products =
        this->parser.get_species(species_filter::products);
    ASSERT0(products.size() == 2);
    ASSERT0(std::find(products.begin(), products.end(), this->p1) != products.end());
    ASSERT0(std::find(products.begin(), products.end(), this->p2) != products.end());
  };

  std::string r1, r2, p1, p2;

  /**
   * @brief Input state for CX reaction tests. Linearly varying densities, temperatures
   * and velocities for reactants. Set product properties for non-symmetric reactions.
   *
   * @return Options
   */
  Options generate_state() override {
    std::string comp_name = "test" + this->lbl;
    Options state{{comp_name,
                   {
                       {"type", this->parser.get_reaction_str()},
                       {"diagnose", true},
                   }},
                  {"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"species",
                   {{this->r1, {{"AA", 1.0}, {"charge", 0.0}}},
                    {this->r2, {{"AA", 1.0}, {"charge", 1.0}}}}}};

    // Reactant properties
    state["species"][this->r1]["density"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(this->logn_min, this->logn_max, linfunc_axis::x), &state,
        mesh);
    state["species"][this->r2]["density"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(this->logn_min, this->logn_max, linfunc_axis::y), &state,
        mesh);
    state["species"][this->r1]["temperature"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(this->logT_min, this->logT_max, linfunc_axis::z), &state,
        mesh);
    state["species"][this->r2]["temperature"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(this->logT_min, this->logT_max, linfunc_axis::x), &state,
        mesh);
    state["species"][this->r1]["velocity"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(this->logv_min, this->logv_max, linfunc_axis::y), &state,
        mesh);
    state["species"][this->r2]["velocity"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(this->logv_min, this->logv_max, linfunc_axis::z), &state,
        mesh);

    // Product properties
    if (!this->parser.is_symmetric()) {
      // For non-symmetric CX, add velocities for product species
      // todo: Currently no attempt to set correct masses and charges for product species!
      state["species"][p1]["AA"] = 1.0;
      state["species"][p1]["charge"] = 1.0;
      state["species"][p1]["velocity"] = FieldFactory::get()->create3D(
          this->gen_lin_field_str(this->logv_min, this->logv_max, linfunc_axis::z),
          &state, mesh);
      state["species"][p2]["AA"] = 1.0;
      state["species"][p2]["charge"] = 0.0;
      state["species"][p2]["velocity"] = FieldFactory::get()->create3D(
          this->gen_lin_field_str(this->logv_min, this->logv_max, linfunc_axis::y),
          &state, mesh);
    }

    return state;
  }
};

class HHpCXTest_noNeutralMomGain : public CXReactionTest<CXReaction> {
public:
  HHpCXTest_noNeutralMomGain()
      : CXReactionTest<CXReaction>("HHpCX_noNeutralMomGain", "h + h+ -> h+ + h") {}
  virtual Options generate_state() override {
    Options state = CXReactionTest<CXReaction>::generate_state();
    state["testHHpCX_noNeutralMomGain"]["no_neutral_cx_mom_gain"] = true;
    return state;
  }
};

class HHpCXTest : public CXReactionTest<CXReaction> {
public:
  HHpCXTest() : CXReactionTest<CXReaction>("HHpCX", "h + h+ -> h+ + h") {}
};

class DDpCXTest : public CXReactionTest<CXReaction> {
public:
  DDpCXTest() : CXReactionTest<CXReaction>("DDpCX", "d + d+ -> d+ + d") {}
};

class TTpCXTest : public CXReactionTest<CXReaction> {
public:
  TTpCXTest() : CXReactionTest<CXReaction>("TTpCX", "t + t+ -> t+ + t") {}
};

class HDpCXTest : public CXReactionTest<CXReaction> {
public:
  HDpCXTest() : CXReactionTest<CXReaction>("HDpCX", "h + d+ -> h+ + d") {}
};

class THpCXTest : public CXReactionTest<CXReaction> {
public:
  THpCXTest() : CXReactionTest<CXReaction>("THpCX", "t + h+ -> t+ + h") {}
};
class DTpCXTest : public CXReactionTest<CXReaction> {
public:
  DTpCXTest() : CXReactionTest<CXReaction>("DTpCX", "d + t+ -> d+ + t") {}
};

} // namespace hermes

#endif // TEST_CX_REACTIONS_H