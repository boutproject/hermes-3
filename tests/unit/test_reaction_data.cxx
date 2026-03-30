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
 * - Invalid metadata key throws at construction
 * - Valid metadata key can be read correctly
 * - Trying to get metadata that wasn't requested at construction throws
 */
TEST(AmjuelDataTest, MetaData) {
  if (std::filesystem::is_directory(test_json_db_path)) {

    std::vector<std::string> invalid_metadata_keys{"some_other_metadata"};
    ASSERT_THROW(AmjuelData("valid_nTfit", valid_options, invalid_metadata_keys),
                 BoutException);

    const std::string valid_metadata_key = "some_metadata";
    std::vector<std::string> valid_metadata_keys{valid_metadata_key};
    AmjuelData instance = AmjuelData("valid_nTfit", valid_options, valid_metadata_keys);

    // Test metadata extraction
    const BoutReal expected_metadata_value = 999.9;
    ASSERT_DOUBLE_EQ(instance.get_metadata(valid_metadata_key), expected_metadata_value);
    // Should throw if we try to get metadata that wasn't requested at construction
    ASSERT_THROW(instance.get_metadata("non_existent_key"), BoutException);
  } else {
    // If tests are run on a filesystem where repo path isn't accessible, just skip
    GTEST_SKIP() << "Couldn't access test json db dir at " << test_json_db_path
                 << ", skipping!";
  }
}

/// @brief Test that src_str returns the expected value
TEST(AmjuelDataTest, Src_str) {
  if (std::filesystem::is_directory(test_json_db_path)) {
    std::string data_ID = "valid_Tfit";
    AmjuelData instance = AmjuelData(data_ID, valid_options, {});
    ASSERT_EQ(instance.src_str(), fmt::format("amjuel_{}", data_ID));
  } else {
    // If tests are run on a filesystem where repo path isn't accessible, just skip
    GTEST_SKIP() << "Couldn't access test json db dir at " << test_json_db_path
                 << ", skipping!";
  }
}

/**
 * Test that
 * - ET fit data can be read
 * - the content of the extracted data (fit type, coeffs) is as expected
 * - the evaluation function throws(not implemented yet)
 */
TEST(AmjuelDataTest, ValidDataDir_ETfit) {
  if (std::filesystem::is_directory(test_json_db_path)) {
    std::vector<std::string> metadata_keys{};
    AmjuelData instance = AmjuelData("valid_ETfit", valid_options, metadata_keys);
    ASSERT_EQ(instance.get_fit_type(), RateParamsTypes::ET);
    ASSERT_EQ(instance.get_coeffs().size(), 2);
    ASSERT_EQ(instance.get_coeffs()[0].size(), 2);
    ASSERT_DOUBLE_EQ(instance.get_coeffs()[0][0], 1.0);
    ASSERT_DOUBLE_EQ(instance.get_coeffs()[0][1], 2.0);
    ASSERT_DOUBLE_EQ(instance.get_coeffs()[1][0], 2.0);
    ASSERT_DOUBLE_EQ(instance.get_coeffs()[1][1], 4.0);

    // Evalulating (E,T) fit not implemented yet, should throw
    ASSERT_THROW(instance.eval_sigma_v_ET(1.0, 1.0), BoutException);
  }
}

/**
 * Test that
 * - 2D coefficient data can be read by specifying a non-default json_database_dir (via
 * the Options object)
 * - the content of the extracted data (fit type, coeffs) is as expected
 * - the evaluation functions return the expected results
 */
TEST(AmjuelDataTest, ValidDataDir_nTfit) {
  if (std::filesystem::is_directory(test_json_db_path)) {
    std::vector<std::string> metadata_keys{};
    AmjuelData instance = AmjuelData("valid_nTfit", valid_options, metadata_keys);
    ASSERT_EQ(instance.get_fit_type(), RateParamsTypes::nT);
    ASSERT_EQ(instance.get_coeffs().size(), 2);
    ASSERT_EQ(instance.get_coeffs()[0].size(), 2);
    ASSERT_DOUBLE_EQ(instance.get_coeffs()[0][0], 1.0);
    ASSERT_DOUBLE_EQ(instance.get_coeffs()[0][1], 2.0);
    ASSERT_DOUBLE_EQ(instance.get_coeffs()[1][0], 2.0);
    ASSERT_DOUBLE_EQ(instance.get_coeffs()[1][1], 4.0);

    // Test evalulation of functions requiring params (n,T)
    const BoutReal sigma_v_expected = 2.7182818284590449e-06;
    BoutReal sigma_v = instance.eval_sigma_v_nT(1.0, 1.0);
    ASSERT_DOUBLE_EQ(sigma_v, sigma_v_expected);

    const BoutReal sigma_vE_expected = 2.7182818284590449e-06;
    BoutReal sigma_vE = instance.eval_sigma_vE_nT(1.0, 1.0);
    ASSERT_DOUBLE_EQ(sigma_vE, sigma_vE_expected);

    // Test evaluation of functions requiring other params; should throw
    ASSERT_THROW(instance.eval_sigma_v_ET(1.0, 1.0), BoutException);
    ASSERT_THROW(instance.eval_sigma_v_T(1.0), BoutException);
  } else {
    // If tests are run on a filesystem where repo path isn't accessible, just skip
    GTEST_SKIP() << "Couldn't access test json db dir at " << test_json_db_path
                 << ", skipping!";
  }
}

/**
 * Test that
 * - the content of the extracted data 1D data (fit type, coeffs) is as expected
 * - the evaluation functions return the expected results
 * - trying to evaluate functions with required parameters that don't match the fit type
 * throws an exception
 */
TEST(AmjuelDataTest, ValidDataDir_Tfit) {
  if (std::filesystem::is_directory(test_json_db_path)) {
    std::vector<std::string> metadata_keys{};
    AmjuelData instance = AmjuelData("valid_Tfit", valid_options, metadata_keys);
    ASSERT_EQ(instance.get_fit_type(), RateParamsTypes::T);
    ASSERT_EQ(instance.get_coeffs().size(), 1);
    ASSERT_EQ(instance.get_coeffs()[0].size(), 2);
    ASSERT_DOUBLE_EQ(instance.get_coeffs()[0][0], 1.0);
    ASSERT_DOUBLE_EQ(instance.get_coeffs()[0][1], 2.0);

    // Test evalulation of functions requiring param (T)
    const BoutReal sigma_v_expected = 2.7182818284590451;
    BoutReal sigma_v = instance.eval_sigma_v_T(1.0);
    ASSERT_DOUBLE_EQ(sigma_v, sigma_v_expected);

    // Test evaluation of functions requiring other params; should throw
    ASSERT_THROW(instance.eval_sigma_v_ET(1.0, 1.0), BoutException);
    ASSERT_THROW(instance.eval_sigma_v_nT(1.0, 1.0), BoutException);
    ASSERT_THROW(instance.eval_sigma_vE_nT(1.0, 1.0), BoutException);
  } else {
    // If tests are run on a filesystem where repo path isn't accessible, just skip
    GTEST_SKIP() << "Couldn't access test json db dir at " << test_json_db_path
                 << ", skipping!";
  }
}
