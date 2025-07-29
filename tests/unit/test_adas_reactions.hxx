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
class ADASIznRecReactionTest : public IznRecReactionTest<RTYPE> {

  static_assert(std::is_base_of<OpenADAS, RTYPE>(),
                "Template arg to ADASReactionTest must derive from OpenADAS");

public:
  ADASIznRecReactionTest(std::string lbl, std::string reaction_str)
      : IznRecReactionTest<RTYPE>(lbl, reaction_str) {}
};

class ADASCIznTest : public ADASIznRecReactionTest<ADASCarbonIonisation<0>> {
public:
  ADASCIznTest()
      : ADASIznRecReactionTest<ADASCarbonIonisation<0>>("CIzn", "c + e -> c+ + 2e") {}
};

class ADASCpIznTest : public ADASIznRecReactionTest<ADASCarbonIonisation<1>> {
public:
  ADASCpIznTest() : ADASIznRecReactionTest("CpIzn", "c+ + e -> c+2 + 2e") {}
};

class ADASC5pIznTest : public ADASIznRecReactionTest<ADASCarbonIonisation<5>> {
public:
  ADASC5pIznTest() : ADASIznRecReactionTest("C5pIzn", "c+5 + e -> c+6 + 2e") {}
};

class ADASCRecTest : public ADASIznRecReactionTest<ADASCarbonRecombination<0>> {
public:
  ADASCRecTest()
      : ADASIznRecReactionTest<ADASCarbonRecombination<0>>("CRec", "c+ + e -> c") {}
};

class ADASCpRecTest : public ADASIznRecReactionTest<ADASCarbonRecombination<1>> {
public:
  ADASCpRecTest()
      : ADASIznRecReactionTest<ADASCarbonRecombination<1>>("CpRec", "c+2 + e -> c+") {}
};

class ADASC5pRecTest : public ADASIznRecReactionTest<ADASCarbonRecombination<5>> {
public:
  ADASC5pRecTest()
      : ADASIznRecReactionTest<ADASCarbonRecombination<5>>("C5pRec", "c+6 + e -> c+5") {}
};

class ADASLiIznTest : public ADASIznRecReactionTest<ADASLithiumIonisation<0>> {
public:
  ADASLiIznTest()
      : ADASIznRecReactionTest<ADASLithiumIonisation<0>>("LiIzn", "li + e -> li+ + 2e") {}
};

class ADASLipIznTest : public ADASIznRecReactionTest<ADASLithiumIonisation<1>> {
public:
  ADASLipIznTest()
      : ADASIznRecReactionTest<ADASLithiumIonisation<1>>("LipIzn",
                                                         "li+ + e -> li+2 + 2e") {}
};

class ADASLi2pIznTest : public ADASIznRecReactionTest<ADASLithiumIonisation<2>> {
public:
  ADASLi2pIznTest()
      : ADASIznRecReactionTest<ADASLithiumIonisation<2>>("Li2pIzn",
                                                         "li+2 + e -> li+3 + 2e") {}
};

class ADASLiRecTest : public ADASIznRecReactionTest<ADASLithiumRecombination<0>> {
public:
  ADASLiRecTest()
      : ADASIznRecReactionTest<ADASLithiumRecombination<0>>("LiRec", "li+ + e -> li") {}
};

class ADASLipRecTest : public ADASIznRecReactionTest<ADASLithiumRecombination<1>> {
public:
  ADASLipRecTest()
      : ADASIznRecReactionTest<ADASLithiumRecombination<1>>("LipRec", "li+2 + e -> li+") {
  }
};

class ADASLi2pRecTest : public ADASIznRecReactionTest<ADASLithiumRecombination<2>> {
public:
  ADASLi2pRecTest()
      : ADASIznRecReactionTest<ADASLithiumRecombination<2>>("Li2pRec",
                                                            "li+3 + e -> li+2") {}
};

class ADASNeIznTest : public ADASIznRecReactionTest<ADASNeonIonisation<0>> {
public:
  ADASNeIznTest()
      : ADASIznRecReactionTest<ADASNeonIonisation<0>>("NeIzn", "ne + e -> ne+ + 2e") {}
};

class ADASNepIznTest : public ADASIznRecReactionTest<ADASNeonIonisation<1>> {
public:
  ADASNepIznTest()
      : ADASIznRecReactionTest<ADASNeonIonisation<1>>("NepIzn", "ne+ + e -> ne+2 + 2e") {}
};

class ADASNe9pIznTest : public ADASIznRecReactionTest<ADASNeonIonisation<9>> {
public:
  ADASNe9pIznTest()
      : ADASIznRecReactionTest<ADASNeonIonisation<9>>("Ne9pIzn",
                                                      "ne+9 + e -> ne+10 + 2e") {}
};

class ADASNeRecTest : public ADASIznRecReactionTest<ADASNeonRecombination<0>> {
public:
  ADASNeRecTest()
      : ADASIznRecReactionTest<ADASNeonRecombination<0>>("NeRec", "ne+ + e -> ne") {}
};

class ADASNepRecTest : public ADASIznRecReactionTest<ADASNeonRecombination<1>> {
public:
  ADASNepRecTest()
      : ADASIznRecReactionTest<ADASNeonRecombination<1>>("NepRec", "ne+2 + e -> ne+") {}
};

class ADASNe9pRecTest : public ADASIznRecReactionTest<ADASNeonRecombination<9>> {
public:
  ADASNe9pRecTest()
      : ADASIznRecReactionTest<ADASNeonRecombination<9>>("Ne9pRec", "ne+10 + e -> ne+9") {
  }
};

#endif