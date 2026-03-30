#include "test_reaction_diagnostic.hxx"

#include <bout/options.hxx>
#include <optional>

// ============================================================================
// ReactionDiagnosticTest Implementation
// ============================================================================

ReactionDiagnostic ReactionDiagnosticTest::default_diag() { return make_diag(); }

ReactionDiagnostic
ReactionDiagnosticTest::make_diag(const std::string& name, const std::string& long_name,
                                  std::optional<ReactionDiagnosticType> type,
                                  const std::string& source,
                                  const std::string& standard_name,
                                  std::optional<DiagnosticTransformerType> transformer) {
  /*
  Set up the root options with normalisation units (required by ReactionDiagnostic
  constructor
   */
  set_root_state();

  // Use provided value if not empty/nullopt, otherwise use default
  std::string actual_name = name.empty() ? this->default_name : name;
  std::string actual_long_name = long_name.empty() ? this->default_long_name : long_name;
  ReactionDiagnosticType actual_type =
      type.has_value() ? type.value() : this->default_type;
  std::string actual_source = source.empty() ? this->default_source : source;

  // For standard_name, if empty use the computed default for the type
  std::string actual_standard_name =
      standard_name.empty() ? ReactionDiagnostic::default_std_name(actual_type)
                            : standard_name;

  // For transformer, if nullopt use the default transformer
  DiagnosticTransformerType actual_transformer =
      transformer.has_value() ? transformer.value() : this->default_transformer;

  return ReactionDiagnostic(actual_name, actual_long_name, actual_type, actual_source,
                            actual_standard_name, actual_transformer);
}

ReactionDiagnostic ReactionDiagnosticTest::make_diag_with_transformer(
    const DiagnosticTransformerType& transformer) {
  return make_diag("", "", std::nullopt, "", "", transformer);
}

void ReactionDiagnosticTest::set_root_state() {
  Options state;
  // Set up the normalisation units required by ReactionDiagnostic constructor, values
  // aren't important for testing.
  state["units"]["inv_meters_cubed"] = 1e19;
  state["units"]["eV"] = 1.0;
  state["units"]["seconds"] = 1e-6;
  Options::root() = state.copy();
}

// ============================================================================
// Tests for ReactionDiagnostic
// ============================================================================

/**
 * @brief Test that add_to_state() correctly adds the diagnostic field to an Options
 * object with the expected attributes and data. Also tests set_data().
 *
 */
TEST_F(ReactionDiagnosticTest, AddToState) {

  // Some of the test expects a collision freq type diagnostic; check the default type
  // hasn't changed.
  EXPECT_TRUE(this->default_type == ReactionDiagnosticType::collision_freq)
      << "Default diagnostic type has changed! Should be "
         "ReactionDiagnosticType::collision_freq.";

  // Create default diagnostic
  ReactionDiagnostic diag = default_diag();

  // Values expexted in the state after adding the diagnostic
  const BoutReal expected_conversion =
      1 / Options::root()["units"]["seconds"].as<BoutReal>();
  const std::string expected_long_name = this->default_long_name;
  const std::string expected_source = this->default_source;
  const std::string expected_standard_name =
      ReactionDiagnostic::default_std_name(this->default_type);
  const std::string expected_units = "s^-1";

  // Set the diagnostic data
  Field3D fld_to_add(1.0);
  diag.set_data(fld_to_add);

  // Create an output state and add the diagnostic to it
  Options output_state;
  diag.add_to_state(output_state);

  // Check that something got added to the output
  EXPECT_TRUE(output_state.isSet(diag.get_name())) << "Diagnostic not added to state";

  // Check that the stored field matches the input field
  Options& data_in_state = output_state[diag.get_name()];
  const Field3D extracted_fld = get<Field3D>(data_in_state);

  EXPECT_TRUE(extracted_fld == fld_to_add)
      << "Field in state does not match the diagnostic field.";

  // Check that the attributes were stored correctly in the state.
  EXPECT_EQ(data_in_state.attributes["long_name"].as<std::string>(), expected_long_name)
      << "Diagnostic long_name attribute doesn't match.";
  EXPECT_EQ(data_in_state.attributes["source"].as<std::string>(), expected_source)
      << "Diagnostic source attribute doesn't match.";
  EXPECT_EQ(data_in_state.attributes["standard_name"].as<std::string>(),
            expected_standard_name)
      << "Diagnostic standard_name attribute doesn't match.";
  EXPECT_EQ(data_in_state.attributes["units"].as<std::string>(), expected_units)
      << "Diagnostic units attribute doesn't match expected value for collision_freq.";
  EXPECT_EQ(data_in_state.attributes["conversion"].as<BoutReal>(), expected_conversion)
      << "Diagnostic conversion attribute doesn't match expected value for "
         "collision_freq.";
}

/**
 * @brief Test that transform() correctly applies a custom transform function
 * to the diagnostic data.
 *
 */
TEST_F(ReactionDiagnosticTest, Transform) {
  Field3D input(1.0);
  Field3D expected_output(2.0);

  // Create a custom transformer that multiplies by 2
  auto multiply_by_2 = [](const Field3D& input) -> Field3D { return 2.0 * input; };

  // Create the diagnostic
  ReactionDiagnostic diag = make_diag_with_transformer(multiply_by_2);

  // Call transform and check the output
  EXPECT_TRUE(diag.transform(input) == expected_output)
      << "Transform output does not match expected value.";
}