#pragma once
#ifndef TEST_REACTIONS_H
#define TEST_REACTIONS_H

#include <filesystem>
#include <regex>
#include <sstream>

#include <bout/constants.hxx>
#include <bout/field_factory.hxx>
#include <bout/options_io.hxx>
#include <gtest/gtest.h>

#include "component.hxx"
#include "cx_reaction.hxx"
#include "izn_rec_reaction.hxx"
#include "reaction.hxx"
#include "reaction_parser.hxx"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx"

enum class linfunc_axis { x, y, z };

/**
 * @brief Prefix increment operator for linfunc_axis.
 *
 * @param current_val The current value of the enum.
 * @return linfunc_axis& The incremented value of the enum.
 */
inline linfunc_axis& operator++(linfunc_axis& current) {
  current = (current == linfunc_axis::z)
                ? linfunc_axis::x
                : static_cast<linfunc_axis>(static_cast<int>(current) + 1);
  return current;
}

/**
 * @brief Postfix increment operator for linfunc_axis.
 *
 * @param current The current value of the enum.
 * @return linfunc_axis The old value of the enum.
 */
inline linfunc_axis operator++(linfunc_axis& current, int) {
  linfunc_axis old = current;
  ++current;
  return old;
}

/// Global mesh
namespace bout::globals {
extern Mesh* mesh;
} // namespace bout::globals

// The unit tests use the global mesh
using namespace bout::globals;

namespace hermes {

/**
 * @brief Base fixture for Reaction tests.
 *
 * @tparam RTYPE The reaction class; must derive from ReactionBase.
 *
 * @todo Change static_assert to require 'Reaction' base later.
 */
template <typename RTYPE>
class ReactionTest : public FakeMeshFixture_tmpl<8, 8, 8> {

  static_assert(std::is_base_of<ReactionBase, RTYPE>(),
                "Template arg to ReactionTest must derive from ReactionBase");

protected:
  ReactionTest(std::string lbl, std::string reaction_str)
      : lbl(lbl), parser(reaction_str) {};

  std::string lbl;
  ReactionParser parser;

  // Ranges to use for input functions
  const BoutReal logn_min = std::log(5e13), logn_max = std::log(2e22);
  const BoutReal logT_min = std::log(0.05), logT_max = std::log(4e4);
  const BoutReal logv_min = std::log(1), logv_max = std::log(100);

  /**
   * @brief Subclasses must override this to generate the input state for the reaction
   * transform.
   *
   * @return Options
   */
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

    // For now we need to manually reset reaction instance counter in tests, otherwise
    // Reaction component construction will fail
    ReactionBase::reset_instance_counter();

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

    // For now we need to manually reset reaction instance counter in tests, otherwise
    // Reaction component construction will fail
    ReactionBase::reset_instance_counter();

    // Run reaction
    RTYPE component = RTYPE("test" + lbl, state, nullptr);
    component.transform(state);

    // Write output state
    std::filesystem::path outpath = ref_data_path();
    std::filesystem::remove(outpath);
    bout::OptionsIO::create(std::string(outpath))->write(state);
  }
};

} // namespace hermes

#endif
