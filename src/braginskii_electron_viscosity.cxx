// Braginskii electron viscosity

#include <cmath>
#include <string>

#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/constants.hxx>
#include <bout/coordinates.hxx>
#include <bout/difops.hxx>
#include <bout/field2d.hxx>
#include <bout/field3d.hxx>
#include <bout/fv_ops.hxx>
#include <bout/mesh.hxx>
#include <bout/options.hxx>
#include <bout/solver.hxx>

#include "../include/braginskii_electron_viscosity.hxx"
#include "../include/component.hxx"

BraginskiiElectronViscosity::BraginskiiElectronViscosity(const std::string& name,
                                                         Options& alloptions, Solver*) {
  auto& options = alloptions[name];

  eta_limit_alpha = options["eta_limit_alpha"]
                        .doc("Viscosity flux limiter coefficient. <0 = turned off")
                        .withDefault(-1.0);

  density_floor = options["density_floor"].doc("Minimum density floor").withDefault(1e-8);

  BoutReal temperature_floor = options["temperature_floor"].doc("Low temperature scale for low_T_diffuse_perp")
    .withDefault<BoutReal>(0.1) / get<BoutReal>(alloptions["units"]["eV"]);

  pressure_floor = density_floor * temperature_floor;                        

  diagnose = options["diagnose"].doc("Output diagnostics?").withDefault<bool>(false);
}

void BraginskiiElectronViscosity::transform(Options& state) {
  AUTO_TRACE();

  Options& species = state["species"]["e"];

  if (!isSetFinal(species["pressure"], "electron_viscosity")) {
    throw BoutException("No electron pressure => Can't calculate electron viscosity");
  }

  if (!isSetFinal(species["velocity"], "electron_viscosity")) {
    throw BoutException("No electron velocity => Can't calculate electron viscosity");
  }

  const Field3D tau = 1. / floor(get<Field3D>(species["collision_frequency"]),1e-10);
  const Field3D P = floor(get<Field3D>(species["pressure"]), pressure_floor);
  const Field3D V = get<Field3D>(species["velocity"]);

  Coordinates* coord = P.getCoordinates();
  const Field3D Bxy = coord->Bxy;
  const Field3D sqrtB = sqrt(Bxy);

  // Parallel electron viscosity
  Field3D eta = (4. / 3) * 0.73 * P * tau;

  if (eta_limit_alpha > 0.) {
    // SOLPS-style flux limiter
    // Values of alpha ~ 0.5 typically

    const Field3D q_cl = eta * abs(Grad_par(V));   // Collisional value
    const Field3D q_fl = eta_limit_alpha * P; // Flux limit

    eta = eta / (1. + floor(q_cl,1e-10) / floor(q_fl,1e-10));

    eta.getMesh()->communicate(eta);
    eta.applyBoundary("neumann");
  }

  // Save term for output diagnostic
  viscosity = sqrtB * FV::Div_par_K_Grad_par(eta / Bxy, sqrtB * V);
  add(species["momentum_source"], viscosity);
}

void BraginskiiElectronViscosity::outputVars(Options& state) {
  AUTO_TRACE();
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto Cs0 = get<BoutReal>(state["Cs0"]);

  if (diagnose) {
    set_with_attrs(state["SNVe_viscosity"], viscosity,
                   {{"time_dimension", "t"},
                    {"units", "kg m^-2 s^-2"},
                    {"conversion", SI::Mp * Nnorm * Cs0 * Omega_ci},
                    {"standard_name", "momentum source"},
                    {"long_name", "electron parallel viscosity"},
                    {"species", "e"},
                    {"source", "electron_viscosity"}});
  }
}
