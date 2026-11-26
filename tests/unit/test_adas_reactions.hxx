#ifndef TEST_ADAS_REACTIONS_H__
#define TEST_ADAS_REACTIONS_H__

#include "test_cx_reactions.hxx"
#include "test_izn_rec_reactions.hxx"

#include "../../include/adas_argon.hxx"
#include "../../include/adas_carbon.hxx"
#include "../../include/adas_lithium.hxx"
#include "../../include/adas_neon.hxx"
#include "../../include/reaction_data.hxx"

#include "../../external/json.hxx"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <iterator>
#include <string>
#include <vector>

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
 * @brief Base fixture for argon tests. The argon OpenADAS files are not included
 * in Hermes-3, so the tests are skipped if they have not been created with
 * json_database/download_openadas.py. For the same reason there are no reference
 * data files for argon: the tests compare with the JSON files instead.
 */
template <typename RTYPE>
class ADASArgonReactionTest : public ADASIznRecReactionTest<RTYPE> {
public:
  using ADASIznRecReactionTest<RTYPE>::ADASIznRecReactionTest;

  void SetUp() override {
    Options options;
    const auto json_db_dir = get_json_db_dir(options);
    for (const auto* name :
         {"scd89_ar.json", "plt89_ar.json", "acd89_ar.json", "prb89_ar.json"}) {
      if (!std::filesystem::exists(json_db_dir / name)) {
        GTEST_SKIP() << (json_db_dir / name).string()
                     << " not found; create it with "
                        "'python3 json_database/download_openadas.py ar'";
      }
    }
  }

  /// Run the reaction with the input state of the regression tests (wide ranges
  /// of density and temperature): all sources must be finite, and the reaction
  /// rate (the density source of the product) must be positive.
  void sources_finite_test() {
    Options state = this->generate_state();
    hermes::ReactionBase::reset_instance_counter();
    RTYPE component("test" + this->lbl, state, nullptr);
    component.transform(state);

    for (const auto& species :
         {std::string("e"), this->heavy_reactant, this->heavy_product}) {
      for (const auto* source : {"density_source", "momentum_source", "energy_source"}) {
        if (!state["species"][species].isSet(source)) {
          continue;
        }
        const auto field = get<Field3D>(state["species"][species][source]);
        BOUT_FOR_SERIAL(i, field.getRegion("RGN_NOBNDRY")) {
          ASSERT_TRUE(std::isfinite(field[i])) << species << ":" << source;
        }
      }
    }
    const auto rate =
        get<Field3D>(state["species"][this->heavy_product]["density_source"]);
    BOUT_FOR_SERIAL(i, rate.getRegion("RGN_NOBNDRY")) { ASSERT_GT(rate[i], 0.0); }
  }

  /// At the node of the rate table closest to Te = 10 eV, ne = 1e19 m^-3, the
  /// reaction rate must be ne * n0 * <sigma v>, with <sigma v> read directly from
  /// the JSON file.
  ///
  /// @param rate_file  The JSON file with the rate coefficients (e.g. scd89_ar.json)
  /// @param level      The ionisation level used by the reaction
  void rate_at_table_node_test(const std::string& rate_file, std::size_t level) {
    // Non-trivial normalisations, so that the unit conversion is tested
    const BoutReal Tnorm = 5.0;  // eV
    const BoutReal Nnorm = 1e18; // m^-3
    const BoutReal tnorm = 1e-6; // s

    Options state = this->generate_state();
    state["units"]["eV"].force(Tnorm);
    state["units"]["inv_meters_cubed"].force(Nnorm);
    state["units"]["seconds"].force(tnorm);

    // Read the table from the JSON file
    const auto path = get_json_db_dir(state) / rate_file;
    std::ifstream json_file(path);
    ASSERT_TRUE(json_file.good()) << path.string();
    const auto data = nlohmann::json::parse(json_file);
    const auto log_temperature = data["log_temperature"].get<std::vector<BoutReal>>();
    const auto log_density = data["log_density"].get<std::vector<BoutReal>>();
    const auto closest = [](const std::vector<BoutReal>& v, BoutReal x) {
      return static_cast<std::size_t>(std::distance(
          v.begin(), std::min_element(v.begin(), v.end(), [x](BoutReal a, BoutReal b) {
            return std::abs(a - x) < std::abs(b - x);
          })));
    };
    const std::size_t iT = closest(log_temperature, 1.0);
    const std::size_t in = closest(log_density, 19.0);
    const auto log_coeff = data["log_coeff"][level][iT][in].get<BoutReal>();

    // The coefficient read by OpenADASRateCoefficient is the one in the file
    const OpenADASRateCoefficient coefficient(path.string(), level);
    EXPECT_EQ(coefficient.log_coeff[iT][in], log_coeff);

    const BoutReal Te_SI = std::pow(10., log_temperature[iT]); // eV
    const BoutReal ne_SI = std::pow(10., log_density[in]);     // m^-3
    const BoutReal n0_SI = 1e17;                               // m^-3
    state["species"]["e"]["density"].force(Field3D(ne_SI / Nnorm));
    state["species"]["e"]["temperature"].force(Field3D(Te_SI / Tnorm));
    state["species"][this->heavy_reactant]["density"].force(Field3D(n0_SI / Nnorm));

    hermes::ReactionBase::reset_instance_counter();
    RTYPE component("test" + this->lbl, state, nullptr);
    component.transform(state);

    // Density source of the product, in m^-3 s^-1
    const auto rate =
        get<Field3D>(state["species"][this->heavy_product]["density_source"]);
    const BoutReal expected = ne_SI * n0_SI * std::pow(10., log_coeff);
    BOUT_FOR_SERIAL(i, rate.getRegion("RGN_NOBNDRY")) {
      ASSERT_NEAR(rate[i] * Nnorm / tnorm, expected, 1e-8 * expected);
    }
  }
};

// Ar izn
class ADASArIznTest : public ADASArgonReactionTest<ADASArgonIonisation<0>> {
public:
  ADASArIznTest()
      : ADASArgonReactionTest<ADASArgonIonisation<0>>("ArIzn", "ar + e -> ar+ + 2e") {}
};

class ADASArpIznTest : public ADASArgonReactionTest<ADASArgonIonisation<1>> {
public:
  ADASArpIznTest()
      : ADASArgonReactionTest<ADASArgonIonisation<1>>("ArpIzn", "ar+ + e -> ar+2 + 2e") {}
};

class ADASAr17pIznTest : public ADASArgonReactionTest<ADASArgonIonisation<17>> {
public:
  ADASAr17pIznTest()
      : ADASArgonReactionTest<ADASArgonIonisation<17>>("Ar17pIzn",
                                                       "ar+17 + e -> ar+18 + 2e") {}
};

// Ar rec
class ADASArRecTest : public ADASArgonReactionTest<ADASArgonRecombination<0>> {
public:
  ADASArRecTest()
      : ADASArgonReactionTest<ADASArgonRecombination<0>>("ArRec", "ar+ + e -> ar") {}
};

class ADASArpRecTest : public ADASArgonReactionTest<ADASArgonRecombination<1>> {
public:
  ADASArpRecTest()
      : ADASArgonReactionTest<ADASArgonRecombination<1>>("ArpRec", "ar+2 + e -> ar+") {}
};

class ADASAr17pRecTest : public ADASArgonReactionTest<ADASArgonRecombination<17>> {
public:
  ADASAr17pRecTest()
      : ADASArgonReactionTest<ADASArgonRecombination<17>>("Ar17pRec",
                                                          "ar+18 + e -> ar+17") {}
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
