#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

#include "../../include/fieldline_geometry.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using FieldlineGeometryTest = FakeMeshFixture;

namespace {
/// Build an Options tree with the units and fieldline_geometry sections
/// required to construct a FieldlineGeometry component.
Options makeOptions(const std::string& lambda_int, const std::string& fieldline_radius,
                    const std::string& poloidal_magnetic_field, bool compute_Btor_from_R,
                    const std::string& toroidal_field_value, bool diagnose = false) {
  Options options;
  options["units"]["meters"] = 1.0;
  options["units"]["Tesla"] = 1.0;

  Options& geo = options["fieldline_geometry"];
  geo["lambda_int"] = lambda_int;
  geo["fieldline_radius"] = fieldline_radius;
  geo["poloidal_magnetic_field"] = poloidal_magnetic_field;
  geo["compute_Btor_from_R"] = compute_Btor_from_R;
  if (compute_Btor_from_R) {
    geo["upstream_toroidal_magnetic_field"] = toroidal_field_value;
  } else {
    geo["toroidal_magnetic_field"] = toroidal_field_value;
  }
  geo["diagnose"] = diagnose;
  return options;
}
} // namespace

TEST_F(FieldlineGeometryTest, CreateComponent) {
  Options options = makeOptions("0.01", "2.0", "0.5", false, "1.0");
  FieldlineGeometry component("test", options, nullptr);
}

TEST_F(FieldlineGeometryTest, ConstantInputsGivenBtor) {
  // compute_Btor_from_R = false: Btor is taken directly from the given function
  Options options = makeOptions("0.01", "2.0", "0.5", false, "1.0", true);
  FieldlineGeometry component("test", options, nullptr);

  Options outputs = {{"rho_s0", 1.0}, {"Bnorm", 1.0}};
  component.outputVars(outputs);

  const BoutReal lambda_int = 0.01;
  const BoutReal R = 2.0;
  const BoutReal Bpol = 0.5;
  const BoutReal Btor = 1.0;
  const BoutReal Btotal = sqrt(Bpol * Bpol + Btor * Btor);
  const BoutReal pitch = Bpol / Btotal;

  ASSERT_TRUE(IsFieldEqual(get<Field3D>(outputs["fieldline_geometry_lambda_int"]),
                           lambda_int, "RGN_NOBNDRY"));
  ASSERT_TRUE(
      IsFieldEqual(get<Field3D>(outputs["fieldline_geometry_Rxy"]), R, "RGN_NOBNDRY"));
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(outputs["fieldline_geometry_Bpxy"]), Bpol,
                           "RGN_NOBNDRY"));
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(outputs["fieldline_geometry_Btxy"]), Btor,
                           "RGN_NOBNDRY"));
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(outputs["fieldline_geometry_Bxy"]), Btotal,
                           "RGN_NOBNDRY"));
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(outputs["fieldline_geometry_magnetic_pitch"]),
                           pitch, "RGN_NOBNDRY"));

  // Inputs are constant in y, so there is no broadening or flux expansion
  ASSERT_TRUE(
      IsFieldEqual(get<Field3D>(outputs["fieldline_geometry_transport_broadening"]), 1.0,
                   "RGN_NOBNDRY"));
  ASSERT_TRUE(
      IsFieldEqual(get<Field3D>(outputs["fieldline_geometry_f_R"]), 1.0, "RGN_NOBNDRY"));
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(outputs["fieldline_geometry_flux_tube_width"]),
                           lambda_int, "RGN_NOBNDRY"));

  // dy = 1 in the fake mesh, so dlpol = dy * Bpol/B = pitch
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(outputs["fieldline_geometry_dlpol"]), pitch,
                           "RGN_NOBNDRY"));

  const BoutReal dlpol = pitch;
  const BoutReal side_area = dlpol * 2.0 * PI * R;
  const BoutReal volume = side_area * lambda_int;
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(outputs["fieldline_geometry_cell_side_area"]),
                           side_area, "RGN_NOBNDRY"));
  ASSERT_TRUE(IsFieldEqual(get<Field3D>(outputs["fieldline_geometry_cell_volume"]),
                           volume, "RGN_NOBNDRY"));
}

TEST_F(FieldlineGeometryTest, ComputeBtorFromConstantR) {
  // compute_Btor_from_R = true, with constant R: Btor = Btor,up * R_up / R = Btor,up
  Options options = makeOptions("0.01", "2.0", "0.5", true, "3.0", true);
  FieldlineGeometry component("test", options, nullptr);

  Options outputs = {{"rho_s0", 1.0}, {"Bnorm", 1.0}};
  component.outputVars(outputs);

  ASSERT_TRUE(
      IsFieldEqual(get<Field3D>(outputs["fieldline_geometry_Btxy"]), 3.0, "RGN_NOBNDRY"));
}

TEST_F(FieldlineGeometryTest, DiagnoseDefaultsToFalse) {
  // diagnose defaults to false: outputVars should not add any variables
  Options options = makeOptions("0.01", "2.0", "0.5", false, "1.0");
  FieldlineGeometry component("test", options, nullptr);

  Options outputs = {{"rho_s0", 1.0}, {"Bnorm", 1.0}};
  component.outputVars(outputs);

  ASSERT_FALSE(outputs.isSet("fieldline_geometry_lambda_int"));
  ASSERT_FALSE(outputs.isSet("fieldline_geometry_cell_volume"));
}

TEST_F(FieldlineGeometryTest, SetsCoordinatesJacobianAndBxy) {
  Options options = makeOptions("0.01", "2.0", "0.5", false, "1.0");
  FieldlineGeometry component("test", options, nullptr);

  Coordinates* coord = mesh->getCoordinates();

  // Bxy is no longer consistent with the new Jacobian, so should be NaN
  // everywhere to stop anyone accidentally using it.
  for (auto i : coord->Bxy.getRegion("RGN_ALL")) {
    ASSERT_TRUE(std::isnan(coord->Bxy[i]));
  }

  // J = 1 / (Beff) / Lnorm. With constant inputs, transport_broadening = 1
  // and Lnorm = 1, so J = 1/Btotal.
  const BoutReal Btotal = sqrt(0.5 * 0.5 + 1.0 * 1.0);
  for (int j = mesh->ystart; j <= mesh->yend; ++j) {
    ASSERT_NEAR(coord->J(0, j), 1.0 / Btotal, 1e-12);
  }

  // J must also be set in the guard cells (not left at the default mesh value),
  // otherwise the metric is discontinuous at domain / inter-processor
  // boundaries. Check the full local range including guard cells.
  for (int j = 0; j < mesh->LocalNy; ++j) {
    ASSERT_NEAR(coord->J(0, j), 1.0 / Btotal, 1e-12);
  }
}

TEST_F(FieldlineGeometryTest, ParallelLengthIncreasesAlongFieldline) {
  Options options = makeOptions("0.01", "2.0", "0.5", false, "1.0", true);
  FieldlineGeometry component("test", options, nullptr);

  Options outputs = {{"rho_s0", 1.0}, {"Bnorm", 1.0}};
  component.outputVars(outputs);

  auto lpar = get<Field3D>(outputs["fieldline_geometry_lpar"]);

  // lpar should increase monotonically from ystart to yend, with dy = 1
  for (int j = mesh->ystart; j < mesh->yend; ++j) {
    ASSERT_NEAR(lpar(0, j + 1, 0) - lpar(0, j, 0), 1.0, 1e-12);
  }
}

TEST_F(FieldlineGeometryTest, LparExpressionSubstitution) {
  // lambda_int can be written in terms of {lpar}, the parallel distance
  // from the upstream boundary. Check that the substitution is applied
  // per-cell, rather than just evaluated once.
  Options options =
      makeOptions("where({lpar} > 1.0, 5.0, 1.0)", "2.0", "0.5", false, "1.0", true);
  FieldlineGeometry component("test", options, nullptr);

  Options outputs = {{"rho_s0", 1.0}, {"Bnorm", 1.0}};
  component.outputVars(outputs);

  auto lpar = get<Field3D>(outputs["fieldline_geometry_lpar"]);
  auto lambda_int = get<Field3D>(outputs["fieldline_geometry_lambda_int"]);

  // Near the upstream boundary, lpar is small, so lambda_int should take the
  // "near" branch; further along the fieldline it should switch to the "far"
  // branch once lpar exceeds 1.0.
  ASSERT_LT(lpar(0, mesh->ystart, 0), 1.0);
  ASSERT_DOUBLE_EQ(lambda_int(0, mesh->ystart, 0), 1.0);

  ASSERT_GT(lpar(0, mesh->yend, 0), 1.0);
  ASSERT_DOUBLE_EQ(lambda_int(0, mesh->yend, 0), 5.0);
}

// The following tests reproduce findings from the manual validation done in
// https://github.com/boutproject/hermes-3/pull/232 (test_fieldline_geometry.pdf):
// the Jacobian is only normalised to 1 for a constant geometry if Bnorm equals
// the upstream field strength, and the geometry output variables satisfy
// self-consistent relationships even for y-varying profiles.

TEST_F(FieldlineGeometryTest, JacobianNotOneWhenBnormDoesNotMatchUpstreamField) {
  // Bpol = 1.0, Btor,up = 3.0 -> upstream |B| = sqrt(10), but Bnorm = 1.0.
  // J = 1 / (Beff / Bnorm) / Lnorm, so J should be 1/sqrt(10), not 1.
  Options options = makeOptions("0.01", "2.0", "1.0", true, "3.0");
  FieldlineGeometry component("test", options, nullptr);

  Coordinates* coord = mesh->getCoordinates();
  const BoutReal expected_J = 1.0 / sqrt(10.0);
  for (int j = mesh->ystart; j <= mesh->yend; ++j) {
    ASSERT_NEAR(coord->J(0, j), expected_J, 1e-12);
  }
}

TEST_F(FieldlineGeometryTest, JacobianIsOneWhenBnormMatchesUpstreamField) {
  // If Bnorm (units:Tesla) is set to the upstream field strength, the
  // Jacobian for a constant geometry should be exactly 1.
  Options options = makeOptions("0.01", "2.0", "1.0", true, "3.0");
  options["units"]["Tesla"].force(sqrt(10.0));

  FieldlineGeometry component("test", options, nullptr);

  Coordinates* coord = mesh->getCoordinates();
  for (int j = mesh->ystart; j <= mesh->yend; ++j) {
    ASSERT_NEAR(coord->J(0, j), 1.0, 1e-12);
  }
}

TEST_F(FieldlineGeometryTest, GeometryFactorsAreSelfConsistentForNonTrivialProfiles) {
  // lambda_int, fieldline_radius and poloidal_magnetic_field all vary with
  // {lpar}, so none of pitch_angle, f_R or transport_broadening are trivially
  // equal to 1. Check that the derived geometry factors are nonetheless
  // self-consistent with each other, as verified manually in the PR.
  Options options = makeOptions("1.0 + 0.1 * {lpar}", "2.0 + 0.1 * {lpar}",
                                "0.5 + 0.01 * {lpar}", true, "1.0", true);
  FieldlineGeometry component("test", options, nullptr);

  Options outputs = {{"rho_s0", 1.0}, {"Bnorm", 1.0}};
  component.outputVars(outputs);

  Coordinates* coord = mesh->getCoordinates();
  auto dy = coord->dy;
  auto pitch_angle = get<Field3D>(outputs["fieldline_geometry_magnetic_pitch"]);
  auto Rxy = get<Field3D>(outputs["fieldline_geometry_Rxy"]);
  auto lambda_int = get<Field3D>(outputs["fieldline_geometry_lambda_int"]);
  auto f_R = get<Field3D>(outputs["fieldline_geometry_f_R"]);
  auto flux_tube_width = get<Field3D>(outputs["fieldline_geometry_flux_tube_width"]);
  auto dlpol = get<Field3D>(outputs["fieldline_geometry_dlpol"]);
  auto cell_side_area = get<Field3D>(outputs["fieldline_geometry_cell_side_area"]);
  auto cell_volume = get<Field3D>(outputs["fieldline_geometry_cell_volume"]);

  // dlpol = dy * Bpol/B
  ASSERT_TRUE(IsFieldEqual(dlpol, dy * pitch_angle, "RGN_NOBNDRY"));

  // flux_tube_width = lambda_int * f_R
  ASSERT_TRUE(IsFieldEqual(flux_tube_width, lambda_int * f_R, "RGN_NOBNDRY"));

  // cell_side_area = dlpol * 2 * pi * R
  ASSERT_TRUE(IsFieldEqual(cell_side_area, dlpol * 2.0 * PI * Rxy, "RGN_NOBNDRY"));

  // cell_volume = cell_side_area * flux_tube_width
  ASSERT_TRUE(IsFieldEqual(cell_volume, cell_side_area * flux_tube_width, "RGN_NOBNDRY"));

  // The profiles are non-trivial, so these factors should not all just be 1
  // (otherwise this test wouldn't be exercising anything the constant-input
  // tests don't already cover).
  ASSERT_FALSE(IsFieldEqual(f_R, 1.0, "RGN_NOBNDRY"));
  ASSERT_FALSE(
      IsFieldEqual(get<Field3D>(outputs["fieldline_geometry_transport_broadening"]), 1.0,
                   "RGN_NOBNDRY"));
}

TEST_F(FieldlineGeometryTest, TransformPublishesGeometryToState) {
  // transform() should publish the geometry quantities into the shared state,
  // so that other components (running afterwards) can read them via get<>().
  Options options = makeOptions("0.01", "2.0", "0.5", false, "1.0");
  FieldlineGeometry component("test", options, nullptr);

  component.transform(options);

  for (const auto& name :
       {"fieldline_geometry_lpar", "fieldline_geometry_lambda_int",
        "fieldline_geometry_magnetic_pitch", "fieldline_geometry_Rxy",
        "fieldline_geometry_Bpxy", "fieldline_geometry_Btxy", "fieldline_geometry_Bxy",
        "fieldline_geometry_transport_broadening", "fieldline_geometry_f_R",
        "fieldline_geometry_flux_tube_width", "fieldline_geometry_dlpol",
        "fieldline_geometry_cell_side_area", "fieldline_geometry_cell_volume"}) {
    ASSERT_TRUE(options.isSet(name)) << name << " was not published to the state";
  }

  // The published major radius should match the requested constant value.
  ASSERT_TRUE(
      IsFieldEqual(get<Field3D>(options["fieldline_geometry_Rxy"]), 2.0, "RGN_NOBNDRY"));
}
