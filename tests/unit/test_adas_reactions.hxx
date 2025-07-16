#ifndef TEST_ADAS_REACTIONS_H__
#define TEST_ADAS_REACTIONS_H__

#include "test_reactions.hxx"

#include "../../include/adas_carbon.hxx"
#include "../../include/adas_lithium.hxx"
#include "../../include/adas_neon.hxx"

/**
 * @brief Base fixture for tests of OpenADAS subclasses.
 *
 * @tparam RTYPE The reaction class; must derive from OpenADAS.
 */
template <typename RTYPE>
class ADASReactionTest : public IznRecReactionTest<RTYPE> {

  static_assert(std::is_base_of<OpenADAS, RTYPE>(),
                "Template arg to ADASReactionTest must derive from OpenADAS");

public:
  ADASReactionTest(std::string lbl, std::string reaction_str)
      : IznRecReactionTest<RTYPE>(lbl, reaction_str) {}
};

class ADASCIznTest : public ADASReactionTest<ADASCarbonIonisation<0>> {
public:
  ADASCIznTest()
      : ADASReactionTest<ADASCarbonIonisation<0>>("CIzn", "c + e -> c+ + 2e") {}
};

class ADASCpIznTest : public ADASReactionTest<ADASCarbonIonisation<1>> {
public:
  ADASCpIznTest() : ADASReactionTest("CpIzn", "c+ + e -> c+2 + 2e") {}
};

class ADASC5pIznTest : public ADASReactionTest<ADASCarbonIonisation<5>> {
public:
  ADASC5pIznTest() : ADASReactionTest("C5pIzn", "c+5 + e -> c+6 + 2e") {}
};

class ADASCRecTest : public ADASReactionTest<ADASCarbonRecombination<0>> {
public:
  ADASCRecTest() : ADASReactionTest<ADASCarbonRecombination<0>>("CRec", "c+ + e -> c") {}
};

class ADASCpRecTest : public ADASReactionTest<ADASCarbonRecombination<1>> {
public:
  ADASCpRecTest()
      : ADASReactionTest<ADASCarbonRecombination<1>>("CpRec", "c+2 + e -> c+") {}
};

class ADASC5pRecTest : public ADASReactionTest<ADASCarbonRecombination<5>> {
public:
  ADASC5pRecTest()
      : ADASReactionTest<ADASCarbonRecombination<5>>("C5pRec", "c+6 + e -> c+5") {}
};

class ADASLiIznTest : public ADASReactionTest<ADASLithiumIonisation<0>> {
public:
  ADASLiIznTest()
      : ADASReactionTest<ADASLithiumIonisation<0>>("LiIzn", "li + e -> li+ + 2e") {}
};

class ADASLipIznTest : public ADASReactionTest<ADASLithiumIonisation<1>> {
public:
  ADASLipIznTest()
      : ADASReactionTest<ADASLithiumIonisation<1>>("LipIzn", "li+ + e -> li+2 + 2e") {}
};

class ADASLi2pIznTest : public ADASReactionTest<ADASLithiumIonisation<2>> {
public:
  ADASLi2pIznTest()
      : ADASReactionTest<ADASLithiumIonisation<2>>("Li2pIzn", "li+2 + e -> li+3 + 2e") {}
};

class ADASLiRecTest : public ADASReactionTest<ADASLithiumRecombination<0>> {
public:
  ADASLiRecTest()
      : ADASReactionTest<ADASLithiumRecombination<0>>("LiRec", "li+ + e -> li") {}
};

class ADASLipRecTest : public ADASReactionTest<ADASLithiumRecombination<1>> {
public:
  ADASLipRecTest()
      : ADASReactionTest<ADASLithiumRecombination<1>>("LipRec", "li+2 + e -> li+") {}
};

class ADASLi2pRecTest : public ADASReactionTest<ADASLithiumRecombination<2>> {
public:
  ADASLi2pRecTest()
      : ADASReactionTest<ADASLithiumRecombination<2>>("Li2pRec", "li+3 + e -> li+2") {}
};

class ADASNeIznTest : public ADASReactionTest<ADASNeonIonisation<0>> {
public:
  ADASNeIznTest()
      : ADASReactionTest<ADASNeonIonisation<0>>("NeIzn", "ne + e -> ne+ + 2e") {}
};

class ADASNepIznTest : public ADASReactionTest<ADASNeonIonisation<1>> {
public:
  ADASNepIznTest()
      : ADASReactionTest<ADASNeonIonisation<1>>("NepIzn", "ne+ + e -> ne+2 + 2e") {}
};

class ADASNe9pIznTest : public ADASReactionTest<ADASNeonIonisation<9>> {
public:
  ADASNe9pIznTest()
      : ADASReactionTest<ADASNeonIonisation<9>>("Ne9pIzn", "ne+9 + e -> ne+10 + 2e") {}
};

class ADASNeRecTest : public ADASReactionTest<ADASNeonRecombination<0>> {
public:
  ADASNeRecTest()
      : ADASReactionTest<ADASNeonRecombination<0>>("NeRec", "ne+ + e -> ne") {}
};

class ADASNepRecTest : public ADASReactionTest<ADASNeonRecombination<1>> {
public:
  ADASNepRecTest()
      : ADASReactionTest<ADASNeonRecombination<1>>("NepRec", "ne+2 + e -> ne+") {}
};

class ADASNe9pRecTest : public ADASReactionTest<ADASNeonRecombination<9>> {
public:
  ADASNe9pRecTest()
      : ADASReactionTest<ADASNeonRecombination<9>>("Ne9pRec", "ne+10 + e -> ne+9") {}
};

#endif