#ifndef TEST_REACTIONS_H__
#define TEST_REACTIONS_H__

#include <filesystem>
#include <regex>
#include <sstream>

#include <bout/constants.hxx>
#include <bout/field_factory.hxx> // For generating functions
#include <bout/options_io.hxx>
#include <gtest/gtest.h>

#include "component.hxx"
#include "test_extras.hxx" // FakeMesh

#include "fake_mesh_fixture.hxx" // IWYU pragma: export

enum class linfunc_axis { x, y, z };

/// Global mesh
namespace bout::globals {
extern Mesh* mesh;
} // namespace bout::globals

// The unit tests use the global mesh
using namespace bout::globals;

/**
 * @brief Base fixture for Reaction tests.
 *
 * @tparam RTYPE The reaction class; must derive from Component.
 *
 * @todo Change static_assert to require 'Reaction' base later.
 */
template <typename RTYPE>
class ReactionTest : public FakeMeshFixture_tmpl<8, 8, 8> {

  static_assert(std::is_base_of<Component, RTYPE>(),
                "Template arg to ReactionTest must derive from Component");

protected:
  ReactionTest(std::string lbl, std::string reaction_str)
      : lbl(lbl), reaction_str(reaction_str){};

  std::string lbl;
  std::string reaction_str;

  // Ranges to use for input functions
  const BoutReal logn_min = std::log(5e13), logn_max = std::log(2e22);
  const BoutReal logT_min = std::log(0.05), logT_max = std::log(4e4);
  const BoutReal logv_min = std::log(1), logv_max = std::log(100);

  // Subclasses must define a function to generate the test input state
  virtual Options generate_state() = 0;

  /**
   * @brief Util to generate an appropriate string to initialise a field with values that
   * increase linearly along axis \p axis_str between \p v_min and \p v_max.
   *
   * @param v_min
   * @param v_max
   * @param axis_str
   * @return std::string field definition suitable for passing to FieldFactory::create3D
   */
  std::string gen_lin_field_str(BoutReal v_min, BoutReal v_max, linfunc_axis axis) {
    // Set coordinate ranges (Use *start and *end indices to exclude guard cells)
    BoutReal axis_min, axis_max;
    std::string axis_str;
    switch (axis) {
    case linfunc_axis::x:
      axis_str = "x";
      axis_min = mesh->GlobalX(mesh->xstart);
      axis_max = mesh->GlobalX(mesh->xend);
      break;
    case linfunc_axis::y:
      axis_str = "y";
      axis_min = TWOPI * mesh->GlobalY(mesh->ystart);
      axis_max = TWOPI * mesh->GlobalY(mesh->yend);
      break;
    case linfunc_axis::z:
      axis_str = "z";
      axis_min = TWOPI * (mesh->zstart) / static_cast<BoutReal>(mesh->LocalNz);
      axis_max = TWOPI * (mesh->zend) / static_cast<BoutReal>(mesh->LocalNz);
      break;
    default:
      axis_min = axis_max = 0;
      throw BoutException("gen_lin_field_str: no such axis type");
      break;
    }

    BoutReal axis_range = axis_max - axis_min;
    BoutReal v_range = v_max - v_min;
    std::stringstream expression;
    expression << "exp(" << v_min << " + (" << axis_str << "-" << axis_min << ")/"
               << axis_range << "*" << v_range << ")";
    return expression.str();
  }

  std::filesystem::path ref_data_path() {
    return std::filesystem::path(__FILE__).parent_path() / "reactions" / (lbl + ".nc");
  }

  /**
   * @brief Compare the values of child nodes in a section of reference data with those in
   * a section of test data. If child nodes are sections themselves, recurse. Node names
   * are assumed to be the same in the two datasets.
   *
   * @param ref_section a section of reference data
   * @param test_section a section of test data
   *
   * @param compare_all_values See \ref sources_regression_test(bool, const int)
   * "sources_regression_test"
   * @param ignore_last_n_sigfigs See \ref sources_regression_test(bool, const int)
   * "sources_regression_test" "sources_regression_test"
   */
  void compare_child_values(const Options& ref_section, const Options& test_section,
                            const bool compare_all_values,
                            const int ignore_last_n_sigfigs) {
    for (auto ref_node : ref_section.getChildren()) {
      std::string ref_node_name = ref_node.first;
      // Using ref_node.second causes the 'as<Field3D>' cast to fail in some case!?
      const Options& ref_node_data = ref_section[ref_node_name];
      const Options& test_node_data = test_section[ref_node_name];

      if (ref_node_data.isValue()) {
        // If node is a value, compare ref, test
        if (compare_all_values || IsSubString(ref_node_name, "_source")) {
          Field3D test_field = test_node_data.as<Field3D>();
          Field3D ref_field = ref_node_data.as<Field3D>();
          ASSERT_TRUE(IsFieldEqualSigFigs(test_field, ref_field, "RGN_NOBNDRY",
                                          ignore_last_n_sigfigs))
              << "'" << this->lbl << "' [" << ref_node_data.str() << "]"
              << " differs from reference data (despite ignoring last "
              << ignore_last_n_sigfigs << " significant digits!)";
        }
      } else {
        // else recurse
        compare_child_values(ref_node_data, test_node_data, compare_all_values,
                             ignore_last_n_sigfigs);
      }
    }
  }

  /**
   * @brief Tests whether calling transform() on a reaction component reproduces the
   * source fields stored in committed reference data files. Subclasses must override
   * generate_state() in order to setup the input to the transform.
   *
   * @param compare_all_values If true, compare all fields in the reference/test states,
   * else just compare the *_source fields.
   *
   * @param ignore_last_n_sigfigs The reference and test fields must be equal at each
   * point, but the last \p ignore_last_n_sigfigs significant digits are ignored in the
   * comparison.
   */
  void sources_regression_test(bool compare_all_values = true,
                               const int ignore_last_n_sigfigs = 6) {

    // Read reference state
    Options ref_state = bout::OptionsIO::create(ref_data_path())->read();

    // Check that the reference state at least has a top-level species section
    const std::string sp_sec_name("species");
    if (!ref_state.isSection(sp_sec_name)) {
      FAIL() << this->lbl << " ref data doesn't have a [" << sp_sec_name << "] section!"
             << std::endl;
    }

    // Generate input state for test
    Options test_state = generate_state();

    // Run reaction
    RTYPE component = RTYPE("test" + lbl, test_state, nullptr);
    component.transform(test_state);

    compare_child_values(ref_state["species"], test_state["species"], compare_all_values,
                         ignore_last_n_sigfigs);
  }

  /**
   * @brief Generate test data by calling generate_state(), then running the reaction
   * transform.
   */
  void generate_data() {
    // Generate input state
    Options state = generate_state();

    // Run reaction
    RTYPE component = RTYPE("test" + lbl, state, nullptr);
    component.transform(state);

    // Write output state
    std::filesystem::path outpath = ref_data_path();
    std::filesystem::remove(outpath);
    bout::OptionsIO::create(std::string(outpath))->write(state);
  }
};

/**
 * @brief Class to test reactions of the form heavy_species1 + <M>e -> heavy_species2 +
 * <N>e where M and N are non-negative integers.
 *
 * @tparam RTYPE
 */
template <typename RTYPE>
class IznRecReactionTest : public ReactionTest<RTYPE> {
protected:
  IznRecReactionTest(std::string lbl, std::string reaction_str)
      : ReactionTest<RTYPE>(lbl, reaction_str) {
    parse_reaction_str(reaction_str);
  };

  std::string heavy_reactant;
  std::string heavy_product;
  std::map<std::string, int> heavy_sp_charges;

  virtual Options generate_state() override final {
    // Default state assumes one heavy reactant, one heavy product. Bail out if that's not
    // the case.
    std::size_t nsp = heavy_sp_charges.size();
    if (nsp != 2) {
      std::stringstream ss;
      ss << "ReactionTest::generate_state assumes exactly one heavy reactant and one "
            "heavy product per reaction;"
         << " found " << nsp << " heavy species."
         << " Override generate_state() in your test class.";
      throw BoutException(ss.str());
    }

    // N.B. No attempt to set the correct masses for heavy species; always set to 1
    std::string comp_name("test" + this->lbl);
    Options state{
        {comp_name, {{"type", this->reaction_str}}},
        {"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
        {"species",
         {{"e", {{"AA", 1. / 1836}, {"velocity", 1.0}}},
          {this->heavy_reactant,
           {{"AA", 1.0}, {"charge", this->heavy_sp_charges.at(this->heavy_reactant)}}},
          {this->heavy_product,
           {{"AA", 1.0},
            {"charge", this->heavy_sp_charges.at(this->heavy_product)},
            {"density", 1.0},
            {"temperature", 1.0}}}}}};

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

  /// Very clunky way to extract the heavy species names in the absence of a
  /// more robust parser
  virtual void parse_reaction_str(const std::string& reaction_str) {
    std::stringstream ss(reaction_str);
    std::string word;
    bool is_reactant = true;
    while (ss >> word) {
      if (word.compare("+") == 0) {
        continue;
      }
      if (word.compare("->") == 0) {
        is_reactant = false;
        continue;
      }

      std::regex pattern("([0-9]*)([a-zA-Z]*)(\\+?\\-?)([0-9]*)");
      std::smatch matches;
      bool has_matches = std::regex_search(word, matches, pattern);
      ASSERT_TRUE(has_matches) << "Unable to parse reaction term " << word << std::endl;
      std::string sp = matches[2].str() + matches[3].str() + matches[4].str();
      if (sp.compare("e") == 0) {
        continue;
      }
      heavy_sp_charges[sp] =
          (matches[4].length() == 0) ? matches[3].length() : stringToInt(matches[4]);

      if (is_reactant) {
        this->heavy_reactant = sp;
      } else {
        this->heavy_product = sp;
      }
    }
  }
};

/**
 * @brief Class to test charge exchange reactions.
 *
 * @tparam RTYPE
 */
template <typename RTYPE>
class CXReactionTest : public ReactionTest<RTYPE> {
protected:
  CXReactionTest(const std::string& lbl, const std::string& reaction_str,
                 const std::string& neutral_sp_in, const std::string& ion_sp_in,
                 const std::string& neutral_sp_out, const std::string& ion_sp_out)
      : ReactionTest<RTYPE>(lbl, reaction_str), neutral_sp_in(neutral_sp_in),
        ion_sp_in(ion_sp_in), neutral_sp_out(neutral_sp_out), ion_sp_out(ion_sp_out){};

  const std::string neutral_sp_in;
  const std::string ion_sp_in;
  const std::string neutral_sp_out;
  const std::string ion_sp_out;

  virtual Options generate_state() override {
    // N.B. No attempt to set the correct masses for heavy species; always set to 1
    // Assume neutral
    std::string comp_name = "test" + this->lbl;
    Options state{{comp_name, {{"type", this->reaction_str}}},
                  {"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"species",
                   {{neutral_sp_in, {{"AA", 1.0}, {"charge", 0.0}}},
                    {ion_sp_in, {{"AA", 1.0}, {"charge", 1.0}}}}}};

    // Linear functions for various fields that are inputs to the reaction transforms
    state["species"][neutral_sp_in]["density"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(this->logn_min, this->logn_max, linfunc_axis::x), &state,
        mesh);
    state["species"][ion_sp_in]["density"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(this->logn_min, this->logn_max, linfunc_axis::y), &state,
        mesh);
    state["species"][neutral_sp_in]["temperature"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(this->logT_min, this->logT_max, linfunc_axis::z), &state,
        mesh);
    state["species"][ion_sp_in]["temperature"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(this->logT_min, this->logT_max, linfunc_axis::x), &state,
        mesh);
    state["species"][neutral_sp_in]["velocity"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(this->logv_min, this->logv_max, linfunc_axis::y), &state,
        mesh);
    state["species"][ion_sp_in]["velocity"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(this->logv_min, this->logv_max, linfunc_axis::z), &state,
        mesh);

    // For non-symmetric CX, add charges, masses, velocities for product species
    if (neutral_sp_out.compare(neutral_sp_in) != 0) {
      state["species"][neutral_sp_out]["AA"] = 1.0;
      state["species"][neutral_sp_out]["charge"] = 0.0;
      state["species"][neutral_sp_out]["velocity"] = FieldFactory::get()->create3D(
          this->gen_lin_field_str(this->logv_min, this->logv_max, linfunc_axis::y),
          &state, mesh);
    }
    if (ion_sp_out.compare(ion_sp_in) != 0) {
      state["species"][ion_sp_out]["AA"] = 1.0;
      state["species"][ion_sp_out]["charge"] = 1.0;
      state["species"][ion_sp_out]["velocity"] = FieldFactory::get()->create3D(
          this->gen_lin_field_str(this->logv_min, this->logv_max, linfunc_axis::z),
          &state, mesh);
    }

    return state;
  }
};

#endif