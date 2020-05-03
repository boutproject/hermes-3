#include <bout/fv_ops.hxx>

#include "../include/diamagnetic_drift.hxx"

using bout::globals::mesh;

DiamagneticDrift::DiamagneticDrift(std::string name, Options &alloptions, Solver *UNUSED(solver)) {

  // Get options for this component
  auto& options = alloptions[name];
  
  bndry_flux = options["bndry_flux"]
                   .doc("Allow fluxes through boundary?")
                   .withDefault<bool>(true);

  // Read curvature vector
  try {
    Curlb_B.covariant = false; // Contravariant
    mesh->get(Curlb_B, "bxcv");
  
  } catch (BoutException &e) {
    try {
      // May be 2D, reading as 3D
      Vector2D curv2d;
      curv2d.covariant = false;
      mesh->get(curv2d, "bxcv");
      Curlb_B = curv2d;
    } catch (BoutException &e) {
      // Need curvature
      throw;
    }
  }

  if (Options::root()["mesh"]["paralleltransform"].withDefault<std::string>(
          "none") == "shifted") {
    Field2D I;
    mesh->get(I, "sinty");
    Curlb_B.z += I * Curlb_B.x;
  }

  // Normalise

  // Get the units
  const auto& units = alloptions["units"];
  BoutReal Bnorm = get<BoutReal>(units["Tesla"]);
  BoutReal Lnorm = get<BoutReal>(units["meters"]);
  
  Curlb_B.x /= Bnorm;
  Curlb_B.y *= SQ(Lnorm);
  Curlb_B.z *= SQ(Lnorm);

  Curlb_B *= 2. / mesh->getCoordinates()->Bxy;
}

void DiamagneticDrift::transform(Options &state) {
  // Iterate through all subsections
  for (auto& kv : state.getChildren()) {
    Options& species = state[kv.first]; // Note: Need non-const

    if (!(species.isSet("charge") and species.isSet("temperature")))
      continue; // Skip, go to next species

    // Calculate diamagnetic drift velocity for this species
    auto q = get<BoutReal>(species["charge"]);
    auto T = get<Field3D>(species["temperature"]);

    // Diamagnetic drift velocity
    Vector3D vD = (T / q) * Curlb_B;

    if (species.isSet("density")) {
      auto N = get<Field3D>(species["density"]);
      
      subtract(species["density_source"],
          FV::Div_f_v(N, vD, bndry_flux));
    }

    if (species.isSet("pressure")) {
      auto P = get<Field3D>(species["pressure"]);
      subtract(species["energy_source"],
               (5. / 2) * FV::Div_f_v(P, vD, bndry_flux));
    }

    if (species.isSet("momentum")) {
      auto NV = get<Field3D>(species["momentum"]);
      subtract(species["momentum_source"],
               FV::Div_f_v(NV, vD, bndry_flux));
    }
  }
}
