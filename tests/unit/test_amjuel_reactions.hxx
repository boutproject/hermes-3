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
    : public CXReactionTest<HydrogenChargeExchangeIsotope<Isotope1, Isotope2>> {
public:
  AmjuelCXTest(std::string lbl, std::string reaction_str)
      : CXReactionTest<HydrogenChargeExchangeIsotope<Isotope1, Isotope2>>(
          lbl, reaction_str, {Isotope1}, {Isotope2, '+'}, {Isotope2}, {Isotope1, '+'}) {}
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