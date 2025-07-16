#ifndef TEST_AMJUEL_REACTIONS_H__
#define TEST_AMJUEL_REACTIONS_H__

#include "test_reactions.hxx"

#include "../../include/amjuel_helium.hxx"
#include "../../include/amjuel_hyd_ionisation.hxx"
#include "../../include/amjuel_hyd_recombination.hxx"
#include "../../include/hydrogen_charge_exchange.hxx"
/**
 * @brief Base fixture for tests of AmjuelReaction subclasses.
 *
 * @tparam RTYPE The reaction class; must derive from AmjuelReaction.
 */
template <typename RTYPE>
class AmjuelReactionTest : public IznRecReactionTest<RTYPE> {

  static_assert(std::is_base_of<AmjuelReaction, RTYPE>(),
                "Template arg to AmjuelReactionTest must derive from AmjuelReaction");

public:
  AmjuelReactionTest(std::string lbl, std::string reaction_str)
      : IznRecReactionTest<RTYPE>(lbl, reaction_str) {}
};

class AmjuelHRecTest : public AmjuelReactionTest<AmjuelHydRecombinationIsotope<'h'>> {
public:
  AmjuelHRecTest()
      : AmjuelReactionTest<AmjuelHydRecombinationIsotope<'h'>>("HRec", "h+ + e -> h") {}
};

class AmjuelDRecTest : public AmjuelReactionTest<AmjuelHydRecombinationIsotope<'d'>> {
public:
  AmjuelDRecTest()
      : AmjuelReactionTest<AmjuelHydRecombinationIsotope<'d'>>("DRec", "d+ + e -> d") {}
};

class AmjuelTRecTest : public AmjuelReactionTest<AmjuelHydRecombinationIsotope<'t'>> {
public:
  AmjuelTRecTest()
      : AmjuelReactionTest<AmjuelHydRecombinationIsotope<'t'>>("TRec", "t+ + e -> t") {}
};

class AmjuelHIznTest : public AmjuelReactionTest<AmjuelHydIonisationIsotope<'h'>> {
public:
  AmjuelHIznTest()
      : AmjuelReactionTest<AmjuelHydIonisationIsotope<'h'>>("HIzn", "h + e -> h+ + 2e") {}
};

class AmjuelDIznTest : public AmjuelReactionTest<AmjuelHydIonisationIsotope<'d'>> {
public:
  AmjuelDIznTest()
      : AmjuelReactionTest<AmjuelHydIonisationIsotope<'d'>>("DIzn", "d + e -> d+ + 2e") {}
};

class AmjuelTIznTest : public AmjuelReactionTest<AmjuelHydIonisationIsotope<'t'>> {
public:
  AmjuelTIznTest()
      : AmjuelReactionTest<AmjuelHydIonisationIsotope<'t'>>("TIzn", "t + e -> t+ + 2e") {}
};

class AmjuelHeIzn01Test : public AmjuelReactionTest<AmjuelHeIonisation01> {
public:
  AmjuelHeIzn01Test()
      : AmjuelReactionTest<AmjuelHeIonisation01>("HeIzn01", "he + e -> he+ + 2e") {}
};

class AmjuelHeRec10Test : public AmjuelReactionTest<AmjuelHeRecombination10> {
public:
  AmjuelHeRec10Test()
      : AmjuelReactionTest<AmjuelHeRecombination10>("HeRec10", "he+ + e -> he") {}
};

/**
 * @brief Test class for CX reactions of the form <neutral_1> + <ion_2> -> <ion_1> +
 * <neutral_2> where _1 and _2 may be the same atom.
 *
 *  *
 * @tparam Isotope1
 * @tparam Isotope2
 */
template <char Isotope1, char Isotope2>
class AmjuelCXTest
    : public ReactionTest<HydrogenChargeExchangeIsotope<Isotope1, Isotope2>> {
protected:
  virtual Options generate_state() override final {
    std::string neutral_sp_in{Isotope1};
    std::string ion_sp_in{Isotope2, '+'};
    std::string neutral_sp_out{Isotope2};
    std::string ion_sp_out{Isotope1, '+'};

    // N.B. No attempt to set the correct masses for heavy species; always set to 1
    // Assume neutral
    std::string comp_name = "test" + this->lbl;
    Options state{{comp_name, {{"type", this->reaction_str}}},
                  {"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"species",
                   {{neutral_sp_in, {{"AA", 1.0}, {"charge", 0.0}}},
                    {ion_sp_in, {{"AA", 1.0}, {"charge", 1.0}}}}}};

    // Density and Temperature ranges (log vals)
    const BoutReal logn_min = std::log(1e14), logn_max = std::log(1e22);
    const BoutReal logT_min = std::log(0.1), logT_max = std::log(2e4);
    const BoutReal logv_min = std::log(1), logv_max = std::log(100);

    // Linear functions for various fields that are inputs to the reaction transforms
    state["species"][neutral_sp_in]["density"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(logn_min, logn_max, linfunc_axis::x), &state, mesh);
    state["species"][ion_sp_in]["density"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(logn_min, logn_max, linfunc_axis::y), &state, mesh);
    state["species"][neutral_sp_in]["temperature"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(logT_min, logT_max, linfunc_axis::z), &state, mesh);
    state["species"][ion_sp_in]["temperature"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(logT_min, logT_max, linfunc_axis::x), &state, mesh);
    state["species"][neutral_sp_in]["velocity"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(logv_min, logv_max, linfunc_axis::y), &state, mesh);
    state["species"][ion_sp_in]["velocity"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(logv_min, logv_max, linfunc_axis::z), &state, mesh);

    // For non-symmetric CX, add charges, masses, velocities for product species
    if (neutral_sp_out.compare(neutral_sp_in) != 0) {
      state["species"][neutral_sp_out]["AA"] = 1.0;
      state["species"][neutral_sp_out]["charge"] = 0.0;
      state["species"][neutral_sp_out]["velocity"] = FieldFactory::get()->create3D(
          this->gen_lin_field_str(logv_min, logv_max, linfunc_axis::y), &state, mesh);
      ;
    }
    if (ion_sp_out.compare(ion_sp_in) != 0) {
      state["species"][ion_sp_out]["AA"] = 1.0;
      state["species"][ion_sp_out]["charge"] = 1.0;
      state["species"][ion_sp_out]["velocity"] = FieldFactory::get()->create3D(
          this->gen_lin_field_str(logv_min, logv_max, linfunc_axis::z), &state, mesh);
      ;
    }

    return state;
  }

public:
  AmjuelCXTest(std::string lbl, std::string reaction_str)
      : ReactionTest<HydrogenChargeExchangeIsotope<Isotope1, Isotope2>>(lbl,
                                                                        reaction_str) {}
};

class HHpCXTest : public AmjuelCXTest<'h', 'h'> {
public:
  HHpCXTest() : AmjuelCXTest<'h', 'h'>("HHpCX", "h + h+ -> h+ + h") {}
};

class DDpCXTest : public AmjuelCXTest<'d', 'd'> {
public:
  DDpCXTest() : AmjuelCXTest<'d', 'd'>("DDpCX", "d + d+ -> d+ + d") {}
};

class TTpCXTest : public AmjuelCXTest<'t', 't'> {
public:
  TTpCXTest() : AmjuelCXTest<'t', 't'>("TTpCX", "t + t+ -> t+ + t") {}
};

class HDpCXTest : public AmjuelCXTest<'h', 'd'> {
public:
  HDpCXTest() : AmjuelCXTest<'h', 'd'>("HDpCX", "h + d+ -> h+ + d") {}
};

class THpCXTest : public AmjuelCXTest<'t', 'h'> {
public:
  THpCXTest() : AmjuelCXTest<'t', 'h'>("THpCX", "t + h+ -> t+ + h") {}
};
class DTpCXTest : public AmjuelCXTest<'d', 't'> {
public:
  DTpCXTest() : AmjuelCXTest<'d', 't'>("DTpCX", "d + t+ -> d+ + t") {}
};

#endif