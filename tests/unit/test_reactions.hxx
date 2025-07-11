#ifndef TEST_REACTIONS_H__
#define TEST_REACTIONS_H__

#include "gtest/gtest.h"
#include <filesystem>
#include <regex>
#include <sstream>

#include <bout/constants.hxx>
#include <bout/field_factory.hxx> // For generating functions
#include <bout/options_io.hxx>

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
      : lbl(lbl), reaction_str(reaction_str) {
    parse_reaction_str(reaction_str);
  };

  std::string lbl;
  std::string heavy_reactant;
  std::string heavy_product;
  std::string reaction_str;
  std::map<std::string, int> heavy_sp_charges;

  virtual Options generate_state() {
    // N.B. No attempt to set the correct masses for heavy species; always set to 1
    std::string comp_name("test" + this->lbl);
    Options state{{comp_name, {{"type", this->reaction_str}}},
                  {"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"species",
                   {{"e", {{"AA", 1. / 1836}, {"velocity", 1.0}}},
                    {this->heavy_reactant,
                     {{"AA", 1.0},
                      {"charge", this->heavy_sp_charges.at(this->heavy_reactant)},
                      {"density", 1.0},
                      {"temperature", 1.0},
                      {"velocity", 1.0}}},
                    {this->heavy_product,
                     {{"AA", 1.0},
                      {"charge", this->heavy_sp_charges.at(this->heavy_product)},
                      {"density", 1.0},
                      {"temperature", 1.0},
                      {"velocity", 1.0}}}}}};

    // Density and Temperature ranges (log vals)
    const BoutReal logn_min = std::log(1e14), logn_max = std::log(1e22);
    const BoutReal logT_min = std::log(0.1), logT_max = std::log(2e4);
    const BoutReal logv_min = std::log(1), logv_max = std::log(100);

    // Linear functions for various fields that are inputs to the reaction transforms
    state["species"]["e"]["density"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(logn_min, logn_max, linfunc_axis::y), &state, mesh);
    state["species"]["e"]["temperature"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(logT_min, logT_max, linfunc_axis::z), &state, mesh);
    state["species"][this->heavy_reactant]["density"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(logn_min, logn_max, linfunc_axis::x), &state, mesh);
    state["species"][this->heavy_reactant]["temperature"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(logT_min, logT_max, linfunc_axis::y), &state, mesh);
    state["species"][this->heavy_reactant]["velocity"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(logv_min, logv_max, linfunc_axis::x), &state, mesh);
    state["species"][this->heavy_product]["velocity"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(logv_min, logv_max, linfunc_axis::z), &state, mesh);
    return state;
  }

  /**
   * @brief Util to generate an appropriate string to initialise a field with values that
   * increase linearly along axis \p axis_str (which has \p axis_ngrid elements) between
   * \p v_min and \p v_max.
   *
   * @param v_min
   * @param v_max
   * @param axis_str
   * @return std::string
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
    }

    BoutReal axis_range = axis_max - axis_min;
    BoutReal v_range = v_max - v_min;
    std::stringstream expression;
    expression << "exp(" << v_min << " + (" << axis_str << "-" << axis_min << ")/"
               << axis_range << "*" << v_range << ")";
    return expression.str();
  }

  /// Very clunky way to extract the heavy species names in the absence of a
  /// more robust parser
  void parse_reaction_str(const std::string& reaction_str) {
    std::stringstream ss(reaction_str);
    std::string word;
    bool is_reactant = true;
    int max_charge_state = std::numeric_limits<int>::min();
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
    auto nsp = heavy_sp_charges.size();
    if (nsp != 2) {
      FAIL() << "ReactionTest currently assumes exactly two heavy species per reaction; "
                "found "
             << nsp << " in the reaction string!";
    }
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