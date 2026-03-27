#include "amjuel_data.hxx"
#include <bout/options.hxx>
#include <filesystem>
#include <gtest/gtest.h>

// Location containing a valid Amjuel json file
static std::filesystem::path test_json_db_path =
    std::filesystem::path(__FILE__).parent_path() / "reactions";

static Options valid_options{
    {"test", {{"type", "x + y+ -> y+ + x"}}},
    {"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
    {"json_database_dir", test_json_db_path}};

/// @brief Test that setting a non-existent json db dir throws.
TEST(AmjuelDataTest, BadCustomDataDir) {
  Options options{{"json_database_dir", "/nonexistent/file/path"}};
  ASSERT_THROW(AmjuelData("dummy_name", options), BoutException);
}

/// @brief Test that setting an invalid label/ID (and therefore json filename) throws.
TEST(AmjuelDataTest, BadFilename) {
  std::string lbl = "non_existent_fname";
  if (std::filesystem::is_directory(test_json_db_path)) {
    ASSERT_THROW(AmjuelData(lbl, valid_options), BoutException);
  } else {
    // If tests are run on a filesystem where repo path isn't accessible, just skip
    GTEST_SKIP() << "Couldn't access test json db dir at " << test_json_db_path
                 << ", skipping!";
  }
}

/// @brief Test that trying to read invalid data throws.
TEST(AmjuelDataTest, InValidData) {
  if (std::filesystem::is_directory(test_json_db_path)) {
    ASSERT_THROW(AmjuelData("invalid", valid_options), BoutException);
  } else {
    // If tests are run on a filesystem where repo path isn't accessible, just skip
    GTEST_SKIP() << "Couldn't access test json db dir at " << test_json_db_path
                 << ", skipping!";
  }
}

/**
 * Test that
 * - data can be read by specifying a non-default json_database_dir (in
 * valid_options)
 * - the content of the extracted data (fit type, coeffs) is as expected
 * - the evaluation functions return the expected results
 * - metadata can be read correctly
 */
TEST(AmjuelDataTest, ValidDataDir) {
  if (std::filesystem::is_directory(test_json_db_path)) {
    const std::string metadata_key = "some_metadata";
    std::vector<std::string> metadata_keys{metadata_key};
    AmjuelData instance = AmjuelData("valid", valid_options, metadata_keys);
    ASSERT_EQ(instance.get_fit_type(), RateParamsTypes::nT);
    ASSERT_EQ(instance.get_coeffs().size(), 2);
    ASSERT_EQ(instance.get_coeffs()[0].size(), 2);
    ASSERT_DOUBLE_EQ(instance.get_coeffs()[0][0], 1.0);
    ASSERT_DOUBLE_EQ(instance.get_coeffs()[0][1], 2.0);
    ASSERT_DOUBLE_EQ(instance.get_coeffs()[1][0], 2.0);
    ASSERT_DOUBLE_EQ(instance.get_coeffs()[1][1], 4.0);
    // Test eval funcs.
    // Data type is nT, but test the other funcs anyway as if it contained the appropriate
    // coeffs.
    const BoutReal nt_eval_expected = 2.7182818284590449e-06;
    BoutReal nT_eval = instance.eval_sigma_vE_nT_impl(1.0, 1.0);
    ASSERT_DOUBLE_EQ(nT_eval, nt_eval_expected);

    const BoutReal T_eval_expected = 2.7182818284590451;
    BoutReal T_eval = instance.eval_sigma_v_T_impl(1.0);
    ASSERT_DOUBLE_EQ(T_eval, T_eval_expected);

    // eval_sigma_v_ET_impl not implemented yet
    ASSERT_THROW(instance.eval_sigma_v_ET_impl(1.0, 1.0), BoutException);

    // Test metadata extraction
    const BoutReal expected_metadata_value = 999.9;
    ASSERT_DOUBLE_EQ(instance.get_metadata(metadata_key), expected_metadata_value);
  } else {
    // If tests are run on a filesystem where repo path isn't accessible, just skip
    GTEST_SKIP() << "Couldn't access test json db dir at " << test_json_db_path
                 << ", skipping!";
  }
}
