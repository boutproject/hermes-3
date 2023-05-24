/// Braginskii electron viscosity

#include <bout/fv_ops.hxx>
#include <bout/mesh.hxx>
#include <bout/difops.hxx>
#include <bout/constants.hxx>

#include "../include/electron_viscosity.hxx"

ElectronViscosity::ElectronViscosity(std::string name, Options& alloptions, Solver*) {
  auto& options = alloptions[name];

  eta_limit_alpha = options["eta_limit_alpha"]
                        .doc("Viscosity flux limiter coefficient. <0 = turned off")
                        .withDefault(-1.0);

  diagnose = options["diagnose"].doc("Output diagnostics?").withDefault<bool>(false);
}

void ElectronViscosity::transform(Options& state) {
  AUTO_TRACE();

  Options& species = state["species"]["e"];

  if (!isSetFinal(species["pressure"], "electron_viscosity")) {
    throw BoutException("No electron pressure => Can't calculate electron viscosity");
  }

  if (!isSetFinal(species["velocity"], "electron_viscosity")) {
    throw BoutException("No electron velocity => Can't calculate electron viscosity");
  }

  const Field3D tau = 1. / get<Field3D>(species["collision_frequency"]);
  const Field3D P = get<Field3D>(species["pressure"]);
  const Field3D V = get<Field3D>(species["velocity"]);

  Coordinates* coord = P.getCoordinates();
  const Field3D Bxy = coord->Bxy;
  const Field3D sqrtB = sqrt(Bxy);

  // Parallel electron viscosity
  Field3D eta = (4. / 3) * 0.73 * P * tau;

  if (eta_limit_alpha > 0.) {
    // SOLPS-style flux limiter
    // Values of alpha ~ 0.5 typically

    const Field3D q_cl = eta * Grad_par(V);   // Collisional value
    const Field3D q_fl = eta_limit_alpha * P; // Flux limit

    eta = eta / (1. + abs(q_cl / q_fl));

    eta.getMesh()->communicate(eta);
    eta.applyBoundary("neumann");
  }

  // Save term for output diagnostic
  viscosity = sqrtB * FV::Div_par_K_Grad_par(eta / Bxy, sqrtB * V);
  add(species["momentum_source"], viscosity);
}

void ElectronViscosity::outputVars(Options& state) {
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


