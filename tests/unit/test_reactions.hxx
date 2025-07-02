#ifndef TEST_REACTIONS_H__
#define TEST_REACTIONS_H__

#include "gtest/gtest.h"
#include <filesystem>
#include <sstream>

#include "component.hxx"
#include "test_extras.hxx" // FakeMesh
#include "bout/options_io.hxx"
#include <bout/constants.hxx>
#include <bout/field_factory.hxx> // For generating functions

#include "fake_mesh_fixture.hxx" // IWYU pragma: export

/// Global mesh
namespace bout::globals {
extern Mesh* mesh;
} // namespace bout::globals

// The unit tests use the global mesh
using namespace bout::globals;

static Options default_opts =
    Options({{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}}});

/**
 * @brief Base fixture for Reaction tests.
 *
 * @tparam RTYPE The reaction class; must derive from Component.
 */
template <typename RTYPE>
class ReactionTest : public FakeMeshFixture {

  static_assert(std::is_base_of<Component, RTYPE>(),
                "Template arg to ReactionTest must derive from Component");

protected:
  ReactionTest(std::string lbl, std::string reaction_str)
      : lbl(lbl), reaction_str(reaction_str){};

  std::string lbl;
  std::string reaction_str;

  virtual Options generate_state() = 0;

  /**
   * @brief Util to generate an appropriate string to initialise a field with values that
   * increase linearly along axis \p axis_str (which has \p axis_ngrid elements) between
   * \p v_min and \p v_max.
   *
   * @param v_min
   * @param v_max
   * @param axis_str
   * @param axis_ngrid
   * @return std::string
   */
  std::string gen_lin_field_str(BoutReal v_min, BoutReal v_max, std::string axis_str,
                                BoutReal axis_min, BoutReal axis_max) {

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
   * @brief Tests whether calling transform() on a reaction component reproduces the
   * source fields stored in committed reference data files. Subclasses must override
   * generate_state() in order to setup the input to the transform.
   *
   * @param check_input_fields Whether to check that the input fields to the transform (in
   * addition to the output source fields) match the reference data.
   *
   * @param ignore_last_n_sigfigs The reference and test fields must be equal at each
   * point, but the last \p ignore_last_n_sigfigs significant digits are ignored in the
   * comparison.
   */
  void sources_regression_test(bool check_input_fields = true,
                               const int ignore_last_n_sigfigs = 3) {

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

    // Loop over all ref_state fields checking that the corresponding test_state field
    // matches
    for (auto sp : ref_state["species"].getChildren()) {
      for (auto fld : sp.second.getChildren()) {

        std::string sp_name = sp.first;
        std::string fld_name = fld.first;

        // Collision frequencies are more complicated; skip them for now.
        if (fld_name.compare("collision_frequencies") == 0) {
          continue;
        }

        if (check_input_fields || IsSubString(fld_name, "_source")) {
          Field3D test_field = test_state["species"][sp.first][fld.first].as<Field3D>();
          Field3D ref_field = ref_state["species"][sp.first][fld.first].as<Field3D>();

          ASSERT_TRUE(IsFieldEqualSigFigs(test_field, ref_field, "RGN_NOBNDRY",
                                          ignore_last_n_sigfigs))
              << "'" << this->lbl << "' [" << sp_name << "/" << fld_name << "]"
              << " differs from reference data (despite ignoring last "
              << ignore_last_n_sigfigs << "significant digits!)";
        }
      }
    }
  }

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

#endif