#ifndef TEST_ADAS_REACTIONS_H__
#define TEST_ADAS_REACTIONS_H__

#include "test_reactions.hxx"

#include "../../include/adas_carbon.hxx"
#include "../../include/adas_lithium.hxx"
#include "../../include/adas_neon.hxx"
#include <cstddef>

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

// C izn
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

// C rec
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

// Li izn
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

// Li rec
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

// Ne izn
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

// Ne rec
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

/**
 * @brief Base fixture class to test ADAS CX reactions
 *
 * @tparam RTYPE the reaction class
 */
template <typename RTYPE>
class ADASCXReactionTest : public CXReactionTest<RTYPE> {
protected:
  ADASCXReactionTest(const std::string& lbl, const std::string& reaction_str,
                     const std::string& neutral_sp_in, const std::string& ion_sp_in,
                     const std::string& neutral_sp_out, const std::string& ion_sp_out)
      : CXReactionTest<RTYPE>(lbl, reaction_str, neutral_sp_in, ion_sp_in, neutral_sp_out,
                              ion_sp_out) {}
  // Add n_e, T_e to the input state
  virtual Options generate_state() override {
    Options state = CXReactionTest<RTYPE>::generate_state();
    state["species"]["e"]["density"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(this->logn_min, this->logn_max, linfunc_axis::z), &state,
        mesh);
    state["species"]["e"]["temperature"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(this->logT_min, this->logT_max, linfunc_axis::x), &state,
        mesh);
    return state;
  }
};

// C CX
template <int level, char Hisotope>
class ADASCCXReactionTest : public ADASCXReactionTest<ADASCarbonCX<level, Hisotope>> {
public:
  ADASCCXReactionTest(std::string lbl, std::string reaction_str)
      : ADASCXReactionTest<ADASCarbonCX<level, Hisotope>>(
          lbl, reaction_str, {carbon_species_name<level + 1>}, {Hisotope},
          {carbon_species_name<level>}, {Hisotope, '+'}) {}
};

class ADASCpHCXTest : public ADASCCXReactionTest<0, 'h'> {
public:
  ADASCpHCXTest() : ADASCCXReactionTest<0, 'h'>("CpHCXCX", "c+ + h -> c + h+") {}
};

class ADASCp3DCXTest : public ADASCCXReactionTest<2, 'd'> {
public:
  ADASCp3DCXTest() : ADASCCXReactionTest<2, 'd'>("Cp3DCX", "c+3 + d -> c+2 + d+") {}
};
class ADASCp5TCXTest : public ADASCCXReactionTest<4, 't'> {
public:
  ADASCp5TCXTest() : ADASCCXReactionTest<4, 't'>("Cp5TCX", "c+5 + t -> c+4 + t+") {}
};

// Li CX
template <std::size_t level, char Hisotope>
class ADASLiCXReactionTest : public ADASCXReactionTest<ADASLithiumCX<level, Hisotope>> {
public:
  ADASLiCXReactionTest(std::string lbl, std::string reaction_str)
      : ADASCXReactionTest<ADASLithiumCX<level, Hisotope>>(
          lbl, reaction_str, {lithium_species_name<level + 1>}, {Hisotope},
          {lithium_species_name<level>}, {Hisotope, '+'}) {}
};

class ADASLipHCXTest : public ADASLiCXReactionTest<0, 'h'> {
public:
  ADASLipHCXTest() : ADASLiCXReactionTest<0, 'h'>("LipHCX", "li+ + h -> li + h+") {}
};

class ADASLip2DCXTest : public ADASLiCXReactionTest<1, 'd'> {
public:
  ADASLip2DCXTest() : ADASLiCXReactionTest<1, 'd'>("Lip2DCX", "li+2 + d -> li+ + d+") {}
};

class ADASLip3TCXTest : public ADASLiCXReactionTest<2, 't'> {
public:
  ADASLip3TCXTest() : ADASLiCXReactionTest<2, 't'>("Lip3TCX", "li+3 + t -> li+2 + t+") {}
};

// Ne CX
template <std::size_t level, char Hisotope>
class ADASNeCXReactionTest : public ADASCXReactionTest<ADASNeonCX<level, Hisotope>> {
public:
  ADASNeCXReactionTest(std::string lbl, std::string reaction_str)
      : ADASCXReactionTest<ADASNeonCX<level, Hisotope>>(
          lbl, reaction_str, {neon_species_name<level + 1>}, {Hisotope},
          {neon_species_name<level>}, {Hisotope, '+'}) {}
};

class ADASNepHCXTest : public ADASNeCXReactionTest<0, 'h'> {
public:
  ADASNepHCXTest() : ADASNeCXReactionTest<0, 'h'>("NepHCX", "ne+ + h -> ne + h+") {}
};

class ADASNep5DCXTest : public ADASNeCXReactionTest<4, 'd'> {
public:
  ADASNep5DCXTest() : ADASNeCXReactionTest<4, 'd'>("Nep5DCX", "ne+5 + d -> ne+4 + d+") {}
};

class ADASNep9TCXTest : public ADASNeCXReactionTest<8, 't'> {
public:
  ADASNep9TCXTest() : ADASNeCXReactionTest<8, 't'>("Nep9TCX", "ne+9 + t -> ne+8 + t+") {}
};
#endif
