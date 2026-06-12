#ifndef TEST_REACTION_DIAGNOSTIC_H
#define TEST_REACTION_DIAGNOSTIC_H

#include <gtest/gtest.h>
#include <optional>

#include "fake_mesh_fixture.hxx"
#include "reaction_diagnostic.hxx"

/// @brief Base fixture for ReactionDiagnostic tests
class ReactionDiagnosticTest : public FakeMeshFixture {
public:
  ReactionDiagnosticTest() {
    // Set up the normalisation units required by ReactionDiagnostic constructor, values
    // aren't important for testing.
    this->units["inv_meters_cubed"] = 1e19;
    this->units["eV"] = 1.0;
    this->units["seconds"] = 1e-6;
  }

protected:
  // Default values for test diagnostics
  std::string default_name{"test_diag"};
  std::string default_long_name{"Test Diagnostic"};
  ReactionDiagnosticType default_type{ReactionDiagnosticType::collision_freq};
  std::string default_source{"test_source"};
  std::string default_standard_name{"test standard name"};
  DiagnosticTransformerType default_transformer{identity};

  /// Units required for ReactionDiagnostic construction
  Options units;

  /// Create a ReactionDiagnostic with default values for all properties
  ReactionDiagnostic default_diag();

  /// Create a ReactionDiagnostic, optionally overriding default values
  ReactionDiagnostic
  make_diag(const std::string& name = "", const std::string& long_name = "",
            std::optional<ReactionDiagnosticType> type = std::nullopt,
            const std::string& source = "", const std::string& standard_name = "",
            std::optional<DiagnosticTransformerType> transformer = std::nullopt);

  /// Create a ReactionDiagnostic with a custom transformer
  ReactionDiagnostic
  make_diag_with_transformer(const DiagnosticTransformerType& transformer);
};

#endif // TEST_REACTION_DIAGNOSTIC_H
