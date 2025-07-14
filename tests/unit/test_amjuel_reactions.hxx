#ifndef TEST_AMJUEL_REACTIONS_H__
#define TEST_AMJUEL_REACTIONS_H__

#include "test_reactions.hxx"

#include "../../include/amjuel_helium.hxx"
#include "../../include/amjuel_hyd_ionisation.hxx"
#include "../../include/amjuel_hyd_recombination.hxx"

/**
 * @brief Base fixture for tests of AmjuelReaction subclasses.
 *
 * @tparam RTYPE The reaction class; must derive from AmjuelReaction.
 */
template <typename RTYPE>
class AmjuelReactionTest : public ReactionTest<RTYPE> {

  static_assert(std::is_base_of<AmjuelReaction, RTYPE>(),
                "Template arg to AmjuelReactionTest must derive from AmjuelReaction");

public:
  AmjuelReactionTest(std::string lbl, std::string reaction_str)
      : ReactionTest<RTYPE>(lbl, reaction_str) {}
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

#endif