/// Ion viscosity model

#include <bout/fv_ops.hxx>
#include <difops.hxx>
#include <bout/mesh.hxx>

#include "../include/ion_viscosity.hxx"

IonViscosity::IonViscosity(std::string name, Options& alloptions, Solver*) {
  auto& options = alloptions[name];

  eta_limit_alpha = options["eta_limit_alpha"]
    .doc("Viscosity flux limiter coefficient. <0 = turned off")
    .withDefault(-1.0);
}

void IonViscosity::transform(Options &state) {
  AUTO_TRACE();

  Options& allspecies = state["species"];

  // Loop through all species
  for (auto& kv : allspecies.getChildren()) {
    if (kv.first == "e") {
      continue; // Skip electrons -> only ions
    }
    Options& species = allspecies[kv.first];

    if (!(isSetFinal(species["pressure"], "ion_viscosity") and
          isSetFinal(species["velocity"], "ion_viscosity"))) {
      // Species doesn't have a pressure and velocity => Skip
      continue;
    }

    const Field3D tau = 1. / get<Field3D>(species["collision_frequency"]);
    const Field3D P = get<Field3D>(species["pressure"]);
    const Field3D V = get<Field3D>(species["velocity"]);

    Coordinates *coord = P.getCoordinates();
    const Field3D Bxy = coord->Bxy;
    const Field3D sqrtB = sqrt(Bxy);

    // Parallel ion viscosity (4/3 * 0.96 coefficient)
    Field3D eta = 1.28 * P * tau;

    if (eta_limit_alpha > 0.) {
      // SOLPS-style flux limiter
      // Values of alpha ~ 0.5 typically

      const Field3D q_cl = eta * Grad_par(V);   // Collisional value
      const Field3D q_fl = eta_limit_alpha * P; // Flux limit

      eta = eta / (1. + abs(q_cl / q_fl));

      eta.getMesh()->communicate(eta);
      eta.applyBoundary("neumann");
    }

    add(species["momentum_source"],
        sqrtB * FV::Div_par_K_Grad_par(eta / Bxy, sqrtB * V));
  }
}
