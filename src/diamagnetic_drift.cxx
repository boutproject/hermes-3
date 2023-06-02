#include <bout/fv_ops.hxx>
#include <bout/vecops.hxx>

#include "../include/diamagnetic_drift.hxx"

using bout::globals::mesh;

DiamagneticDrift::DiamagneticDrift(std::string name, Options& alloptions,
                                   Solver* UNUSED(solver)) {

  // Get options for this component
  auto& options = alloptions[name];

  bndry_flux =
      options["bndry_flux"].doc("Allow fluxes through boundary?").withDefault<bool>(true);

  diamag_form = options["diamag_form"]
    .doc("Form of diamagnetic drift: 0 = gradient; 1 = divergence")
    .withDefault(Field3D(1.0));

  // Read curvature vector
  Curlb_B.covariant = false; // Contravariant
  if (mesh->get(Curlb_B, "bxcv")) {
    Curlb_B.x = Curlb_B.y = Curlb_B.z = 0.0;
  }

  Options& paralleltransform = Options::root()["mesh"]["paralleltransform"];
  if (paralleltransform.isSet("type") and
      paralleltransform["type"].as<std::string>() == "shifted") {
    Field2D I;
    if (mesh->get(I, "sinty")) {
      I = 0.0;
    }
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

  // Set drift to zero through sheath boundaries.
  // Flux through those cell faces should be set by sheath.
  for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
    for (int jz = 0; jz < mesh->LocalNz; jz++) {
      Curlb_B.y(r.ind, mesh->ystart - 1, jz) = -Curlb_B.y(r.ind, mesh->ystart, jz);
    }
  }
  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
    for (int jz = 0; jz < mesh->LocalNz; jz++) {
      Curlb_B.y(r.ind, mesh->yend + 1, jz) = -Curlb_B.y(r.ind, mesh->yend, jz);
    }
  }
}

void DiamagneticDrift::transform(Options& state) {
  // Iterate through all subsections
  Options& allspecies = state["species"];

  for (auto& kv : allspecies.getChildren()) {
    Options& species = allspecies[kv.first]; // Note: Need non-const

    if (!(species.isSet("charge") and species.isSet("temperature")))
      continue; // Skip, go to next species

    // Calculate diamagnetic drift velocity for this species
    auto q = get<BoutReal>(species["charge"]);
    auto T = GET_VALUE(Field3D, species["temperature"]);

    // Diamagnetic drift velocity
    Vector3D vD = (T / q) * Curlb_B;

    if (IS_SET(species["density"])) {
      auto N = GET_VALUE(Field3D, species["density"]);

      // Divergence form: Div(n v_D)
      Field3D div_form = FV::Div_f_v(N, vD, bndry_flux);
      // Gradient form: Curlb_B dot Grad(N T / q)
      Field3D grad_form = Curlb_B * Grad(N * T / q);

      subtract(species["density_source"], diamag_form * div_form + (1. - diamag_form) * grad_form);
    }

    if (IS_SET(species["pressure"])) {
      auto P = get<Field3D>(species["pressure"]);

      Field3D div_form = FV::Div_f_v(P, vD, bndry_flux);
      Field3D grad_form = Curlb_B * Grad(P * T / q);
      subtract(species["energy_source"], (5. / 2) * (diamag_form * div_form + (1. - diamag_form) * grad_form));
    }

    if (IS_SET(species["momentum"])) {
      auto NV = get<Field3D>(species["momentum"]);
      Field3D div_form = FV::Div_f_v(NV, vD, bndry_flux);
      Field3D grad_form = Curlb_B * Grad(NV * T / q);
      subtract(species["momentum_source"], diamag_form * div_form + (1. - diamag_form) * grad_form);
    }
  }
}
