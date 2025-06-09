#include <bout/fv_ops.hxx>
#include <bout/vecops.hxx>

#include "../include/centrifugal_drift.hxx"

using bout::globals::mesh;

namespace {
// Convert a Vector2D to Vector3D
Vector3D toVector3D(const Vector2D& v2d) {
  Vector3D v3d;
  v3d.covariant = v2d.covariant;
  v3d.x = v2d.x;
  v3d.y = v2d.y;
  v3d.z = v2d.z;
  return v3d;
}
} // namespace

CentrifugalDrift::CentrifugalDrift(std::string name, Options& alloptions,
                                   Solver* UNUSED(solver)) {
  AUTO_TRACE();

  // Get the units
  const auto& units = alloptions["units"];
  BoutReal seconds = get<BoutReal>(units["seconds"]);
  BoutReal meters = get<BoutReal>(units["meters"]);

  // Get options for this component
  auto& options = alloptions[name];

  bndry_flux =
      options["bndry_flux"].doc("Allow fluxes through boundary?").withDefault<bool>(true);

  // Read acceleration in the X direction and normalise
  auto acceleration_x =
      options["acceleration"].doc("Acceleration in the X direction [m/s^2]").as<Field2D>()
    * (SQ(seconds) / meters);

  const auto metric = mesh->getCoordinates();

  // Convert to dimensional basis vector
  // e.g. X component in terms of the input:
  // a_x = a dot e_x
  //     = e_x^hat sqrt(e_x dot e_x)
  //     = (a dot e_x^hat) sqrt(g_xx)
  acceleration.covariant = true; // Covariant components
  acceleration.x = acceleration_x * sqrt(metric->g_11);
  acceleration.y = 0.0;
  acceleration.z = 0.0;

  Vector2D b_B; // b unit vector divided by magnitude of B
  b_B.covariant = false;
  b_B.x = 0.;
  b_B.y = 1. / (metric->J * metric->Bxy);
  b_B.z = 0.;

  // Calculate cross product to get (acceleration x b) / B drift vector
  axb_B = cross(acceleration, b_B);
  axb_B.toContravariant(); // Make sure it's contravariant
}

void CentrifugalDrift::transform(Options& state) {
  // Iterate through all species
  Options& allspecies = state["species"];

  for (auto& kv : allspecies.getChildren()) {
    Options& species = allspecies[kv.first]; // Note: Need non-const

    if (!(IS_SET(species["charge"]) and IS_SET(species["AA"]))) {
      // No charge or no mass
      continue; // Skip, go to next species
    }

    auto q = get<BoutReal>(species["charge"]);
    if (fabs(q) < 1e-5) {
      continue;
    }
    auto M = get<BoutReal>(species["AA"]);

    // Centrifugal drift velocity
    Vector3D vC = toVector3D((M / q) * axb_B);

    if (IS_SET(species["density"])) {
      // Drift of particles
      auto N = GET_VALUE(Field3D, species["density"]);
      Field3D div_flow = FV::Div_f_v(N, vC, bndry_flux);
      subtract(species["density_source"], div_flow);

      // Divergence of current that will drive vorticity
      add(state["fields"]["DivJextra"], q * div_flow);
    }

    if (IS_SET(species["pressure"])) {
      // Drift of thermal energy
      auto P = get<Field3D>(species["pressure"]);
      subtract(species["energy_source"], (5. / 2) * FV::Div_f_v(P, vC, bndry_flux));
    }

    if (IS_SET(species["momentum"])) {
      // Drift of parallel momentum
      auto NV = get<Field3D>(species["momentum"]);
      subtract(species["momentum_source"], FV::Div_f_v(NV, vC, bndry_flux));
    }
  }
}

void CentrifugalDrift::outputVars(Options& state) {
  AUTO_TRACE();

  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto rho_s0 = get<BoutReal>(state["rho_s0"]);
  auto Bnorm = get<BoutReal>(state["Bnorm"]);

  // Note: Vectors can't be added to Options, but components could be

  // set_with_attrs(state["acceleration"], acceleration,
  //                  {{"units", "m/s^2"},
  //                   {"conversion", rho_s0 * SQ(Omega_ci)},
  //                   {"long_name", "Acceleration"},
  //                   {"source", "centrifugal_drift"}});

  // set_with_attrs(state["axb_B"], axb_B,
  //                {{"units", "m/s^2/T"},
  //                 {"conversion", rho_s0 * SQ(Omega_ci) / Bnorm},
  //                 {"long_name", "acceleration x b / B"},
  //                 {"source", "centrifugal_drift"}});
}
