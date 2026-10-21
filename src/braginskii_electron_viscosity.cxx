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
#include "../include/free_streaming.hxx"
#include "../include/hermes_utils.hxx"

BraginskiiElectronViscosity::BraginskiiElectronViscosity(const std::string& name,
                                                         Options& alloptions, Solver*)
    : NamedComponent(name,
                     {readIfSet("species:e:pressure"), readIfSet("species:e:velocity"),
                      readOnly("species:e:collision_frequency"),
                      readWrite("species:e:momentum_source")}) {
  auto& options = alloptions[name];

  eta_limit_alpha = options["eta_limit_alpha"]
                        .doc("Viscosity flux limiter coefficient. <0 = turned off")
                        .withDefault(-1.0);

  eta_free_streaming =
      options["eta_free_streaming"]
          .doc("Reduce eta by 1 / (1 + lambda_e / (c L_par)) for collisionless electrons? 1/L_par = 1/L_line (+ |Grad_par Te|/Te, see eta_free_streaming_gradient)")
          .withDefault<bool>(false);
  if (eta_free_streaming) {
    eta_free_streaming_c = options["eta_free_streaming_c"]
                               .doc("Coefficient c in 1 / (1 + lambda_e / (c L_par))")
                               .withDefault(1.0);
    eta_free_streaming_gradient =
        options["eta_free_streaming_gradient"]
            .doc("Add the local gradient length |Grad_par Te|/Te to 1/L_par?")
            .withDefault<bool>(false);
    const BoutReal rho_s0 = alloptions["units"]["meters"];
    if (options.isSet("eta_free_streaming_length")) {
      L_line = options["eta_free_streaming_length"]
                  .doc("Field line length L_line [m]. Default: computed from the grid")
                  .as<BoutReal>()
              / rho_s0;
    } else {
      L_line = fieldLineLength(); // grids with branch cuts: set eta_free_streaming_length
    }
    output.write("\tElectron viscosity field line length L_line: {:e} - {:e} m\n",
                 min(L_line, true) * rho_s0, max(L_line, true) * rho_s0);
    setPermissions(readOnly("species:e:temperature"));
    setPermissions(readOnly("species:e:AA"));
  }

  diagnose = options["diagnose"].doc("Output diagnostics?").withDefault<bool>(false);
}

void BraginskiiElectronViscosity::transform_impl(GuardedOptions& state) {

  GuardedOptions species = state["species"]["e"];

  if (!isSetFinal(species["pressure"], "electron_viscosity")) {
    throw BoutException("No electron pressure => Can't calculate electron viscosity");
  }

  if (!isSetFinal(species["velocity"], "electron_viscosity")) {
    throw BoutException("No electron velocity => Can't calculate electron viscosity");
  }

  // Floor: nu = 0 (e.g. density <= 0 in a trial state) would give tau = inf and NaN eta
  const Field3D tau = 1. / softFloor(get<Field3D>(species["collision_frequency"]), 1e-10);
  const Field3D P = get<Field3D>(species["pressure"]);
  const Field3D V = get<Field3D>(species["velocity"]);

  Coordinates* coord = P.getCoordinates();

  const bout::FieldMetricParallel Bxy = coord->Bxy();

  // If not allocated calculate, otherwise skip as already done
  if (!sqrtB.isAllocated()) {
    sqrtB = sqrt(Bxy);
  }

  // Parallel electron viscosity
  Field3D eta = (4. / 3) * 0.73 * P * tau;

  if (eta_free_streaming) {
    // Braginskii eta assumes lambda_e << L_par; reduce it smoothly when collisionless
    const Field3D Te = get<Field3D>(species["temperature"]);
    const BoutReal AA = get<BoutReal>(species["AA"]);
    const Field3D lambda_e = sqrt(floor(Te, 0.0) / AA) * tau; // v_te * tau_e
    eta *= freeStreamingFactor(Te, lambda_e, L_line, eta_free_streaming_c,
                               eta_free_streaming_gradient);
  }

  if (eta_limit_alpha > 0.) {
    // SOLPS-style flux limiter
    // Values of alpha ~ 0.5 typically

    const Field3D q_cl = eta * Grad_par(V);   // Collisional value
    // Flux limit; floor avoids q_cl / q_fl = 0/0 where P = 0
    const Field3D q_fl = eta_limit_alpha * floor(P, 1e-6);

    eta = eta / (1. + abs(q_cl / q_fl));
  }

  // Fci needs communication for the parallel slices; Grad_par in the limiter and
  // free-streaming factor leaves the guard cells unset
  if (P.isFci() || eta_limit_alpha > 0.0 || eta_free_streaming) {
    eta.applyBoundary("neumann");
    eta.getMesh()->communicate(eta);
    eta.applyParallelBoundary("parallel_neumann_o2");
  }

  // Save term for output diagnostic
  viscosity =
      sqrtB
      * FV::Div_par_K_Grad_par(Field3DParallel{eta / Bxy}, Field3DParallel{sqrtB * V});
  add(species["momentum_source"], viscosity);
}

void BraginskiiElectronViscosity::outputVars(Options& state) {
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
