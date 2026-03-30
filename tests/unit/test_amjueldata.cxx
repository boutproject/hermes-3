#include <filesystem>

#include <bout/boutexception.hxx>
#include <bout/options.hxx>
#include <gtest/gtest.h>

#include "test_amjuel_reactions.hxx"

/*
 * N.B. AmjuelData is tested via parent AmjuelReaction instances. Testing it
 * directly would require intrusive 'friend' declarations)
 */

// Location containing a valid Amjuel json file
static std::filesystem::path test_json_db_path =
    std::filesystem::path(__FILE__).parent_path() / "reactions";

static Options valid_options{
    {"test", {{"type", "x + y+ -> y+ + x"}}},
    {"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
    {"json_database_dir", test_json_db_path}};

/// @brief Test that setting an invalid json db dir throws.
TEST(AmjuelDataTest, InvalidCustomDataDir) {
  Options options{{"json_database_dir", "/invalid/file/path"}};
  ASSERT_THROW(AmjuelReaction("dummy_name", "iz", "test", options), BoutException);
}

/// @brief Test that setting an invalid Amjuel label (and therefore json filename) throws.
TEST(AmjuelDataTest, InValidFilename) {
  std::string invalid_amjuel_lbl = "invalid_lbl";
  if (std::filesystem::is_directory(test_json_db_path)) {
    ASSERT_THROW(AmjuelReaction("test", "valid", invalid_amjuel_lbl, valid_options),
                 BoutException);
  } else {
    // If tests are run on a filesystem where repo path isn't accessible, just skip
    GTEST_SKIP() << "Couldn't access test json db dir at " << test_json_db_path
                 << ", skipping!";
  }
}

/// @brief Test that trying to read invalid data throws.
TEST(AmjuelDataTest, InValidData) {
  if (std::filesystem::is_directory(test_json_db_path)) {
    ASSERT_THROW(AmjuelReaction("test", "invalid", "test", valid_options), BoutException);
  } else {
    // If tests are run on a filesystem where repo path isn't accessible, just skip
    GTEST_SKIP() << "Couldn't access test json db dir at " << test_json_db_path
                 << ", skipping!";
  }
}

/// Test that json_database_dir can be overridden with a valid path
TEST(AmjuelDataTest, ValidCustomDataDir) {
  if (std::filesystem::is_directory(test_json_db_path)) {
    ASSERT_NO_THROW(AmjuelReaction("test", "valid", "test", valid_options));
  } else {
    // If tests are run on a filesystem where repo path isn't accessible, just skip
    GTEST_SKIP() << "Couldn't access test json db dir at " << test_json_db_path
                 << ", skipping!";
  }
}

/// Test that reading data without <sigma v E> coefficients works
TEST(AmjuelDataTest, ValidNoSigmavEData) {
  if (std::filesystem::is_directory(test_json_db_path)) {
    ASSERT_NO_THROW(AmjuelReaction("test", "valid_no-sigma-v-E", "test", valid_options));
  } else {
    // If tests are run on a filesystem where repo path isn't accessible, just skip
    GTEST_SKIP() << "Couldn't access test json db dir at " << test_json_db_path
                 << ", skipping!";
  }
}