#ifndef TEST_ADAS_REACTIONS_H__
#define TEST_ADAS_REACTIONS_H__

#include "test_cx_reactions.hxx"
#include "test_izn_rec_reactions.hxx"

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
class ADASIznRecReactionTest : public hermes::IznRecReactionTest<RTYPE> {

  static_assert(std::is_base_of<OpenADAS, RTYPE>(),
                "Template arg to ADASReactionTest must derive from OpenADAS");

public:
  ADASIznRecReactionTest(std::string lbl, std::string reaction_str)
      : hermes::IznRecReactionTest<RTYPE>(lbl, reaction_str) {}

protected:
  /**
   * @brief Check radiated power in SI units against an ADAS radiation table node
   */
  void radiated_power_si_test() {
    // Non-trivial normalisations, so that the unit conversion is tested
    const BoutReal Tnorm = 5.0;  // eV
    const BoutReal Nnorm = 1e18; // m^-3
    const BoutReal tnorm = 1e-6; // s

    // Set up the state and overwrite the units with the test-specific values
    Options state = this->generate_state();
    state["units"]["eV"].force(Tnorm);
    state["units"]["inv_meters_cubed"].force(Nnorm);
    state["units"]["seconds"].force(tnorm);

    // Create the component
    hermes::ReactionBase::reset_instance_counter();
    RTYPE component("test" + this->lbl, state, nullptr);

    /* Use nodal values of ne,Te in the radiation coefficient table to avoid interpolation errors.
       Choose the node closest to Te = 10 eV, ne = 1e19 m^-3.
    */
    const auto [Te_SI, ne_SI, radiation_coefficient] =
        component.nearest_radiation_coef_node(10., 1e19);
    const BoutReal n1_SI = 1e17; // m^-3

    // Overwrite the species densities and temperatures in the state with the test-specific values
    state["species"]["e"]["density"].force(Field3D(ne_SI / Nnorm));
    state["species"]["e"]["temperature"].force(Field3D(Te_SI / Tnorm));
    state["species"][this->heavy_reactant]["density"].force(Field3D(n1_SI / Nnorm));

    // Apply the reaction transform
    component.transform(state);

    // Check that the radiated power matches the expected value
    const Field3D rate =
        get<Field3D>(state["species"][this->heavy_product]["density_source"]);
    const Field3D energy = get<Field3D>(state["species"]["e"]["energy_source"]);
    const BoutReal rate_to_SI = Nnorm / tnorm;                   // m^-3 s^-1
    const BoutReal power_to_SI = SI::qe * Tnorm * Nnorm / tnorm; // W m^-3

    // Electron energy loss = radiated power - electron heating * reaction rate
    const BoutReal expected = ne_SI * n1_SI * radiation_coefficient; // W m^-3
    const BoutReal electron_heating = component.get_electron_heating();
    BOUT_FOR_SERIAL(i, energy.getRegion("RGN_NOBNDRY")) {
      const BoutReal radiated =
          -energy[i] * power_to_SI + SI::qe * electron_heating * rate[i] * rate_to_SI;
      ASSERT_NEAR(radiated, expected, 1e-8 * expected);
    }
  }
};

// C izn
class ADASCIznTest : public ADASIznRecReactionTest<ADASCarbonIonisation<0>> {
public:
  ADASCIznTest()
      : ADASIznRecReactionTest<ADASCarbonIonisation<0>>("CIzn", "c + e -> c+ + 2e") {}
};

class ADASCpIznTest : public ADASIznRecReactionTest<ADASCarbonIonisation<1>> {
public:
  ADASCpIznTest()
      : ADASIznRecReactionTest<ADASCarbonIonisation<1>>("CpIzn", "c+ + e -> c+2 + 2e") {}
};

class ADASC5pIznTest : public ADASIznRecReactionTest<ADASCarbonIonisation<5>> {
public:
  ADASC5pIznTest()
      : ADASIznRecReactionTest<ADASCarbonIonisation<5>>("C5pIzn", "c+5 + e -> c+6 + 2e") {
  }
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
class ADASCXReactionTest : public hermes::CXReactionTest<RTYPE> {
protected:
  ADASCXReactionTest(const std::string& lbl, const std::string& reaction_str)
      : hermes::CXReactionTest<RTYPE>(lbl, reaction_str) {}
  // Add n_e, T_e to the input state
  virtual Options generate_state() override {
    Options state = hermes::CXReactionTest<RTYPE>::generate_state();
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
      : ADASCXReactionTest<ADASCarbonCX<level, Hisotope>>(lbl, reaction_str) {}
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
      : ADASCXReactionTest<ADASLithiumCX<level, Hisotope>>(lbl, reaction_str) {}
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
      : ADASCXReactionTest<ADASNeonCX<level, Hisotope>>(lbl, reaction_str) {}
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
