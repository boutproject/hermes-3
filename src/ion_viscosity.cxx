/// Ion viscosity model

#include <bout/fv_ops.hxx>

#include "../include/ion_viscosity.hxx"

IonViscosity::IonViscosity(std::string name, Options& alloptions, Solver*) {

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

    add(species["momentum_source"], 1.28 * sqrtB *
        FV::Div_par_K_Grad_par(P * tau / Bxy, sqrtB * V));
  }
}
