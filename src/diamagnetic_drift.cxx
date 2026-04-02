#include <bout/fv_ops.hxx>
#include <bout/vecops.hxx>

#include "../include/diamagnetic_drift.hxx"

using bout::globals::mesh;

DiamagneticDrift::DiamagneticDrift(std::string name, Options& alloptions,
                                   Solver* UNUSED(solver))
    : Component({readIfSet("species:{all_species}:{input}"),
                 readWrite("species:{all_species}:{output}")}) {

  // Get options for this component
  auto& options = alloptions[name];

  bndry_flux = options["bndry_flux"]
                   .doc("Allow fluxes through boundary?")
                   .withDefault<bool>(false);

  divergence_form = options["divergence_form"]
                        .doc("Use divergence form of diamagnetic drift?")
                        .withDefault<bool>(true);

  average_core =
      options["average_core"].doc("Average around core boundary?").withDefault<bool>(true)
      and
      // Only average if on core boundary
      mesh->periodicY(mesh->xstart) and mesh->firstX();

  if (average_core) {
    const auto* coords = mesh->getCoordinates();
    this->cell_volume = coords->dx * coords->dy * coords->dz * coords->J;
    BoutReal local_core_volume = 0.0;
    for (int jy = mesh->ystart; jy <= mesh->yend; ++jy) {
      local_core_volume += this->cell_volume(mesh->xstart, jy);
    }
    // Sum over processors on core boundary
    MPI_Allreduce(&local_core_volume, &this->core_ring_volume, 1, MPI_DOUBLE, MPI_SUM,
                  mesh->getYcomm(0));
  }

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
    Curlb_B.y(r.ind, mesh->ystart - 1) = -Curlb_B.y(r.ind, mesh->ystart);
  }
  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
    Curlb_B.y(r.ind, mesh->yend + 1) = -Curlb_B.y(r.ind, mesh->yend);
  }

  // FIXME: density, pressure, and momentum will not be read even if
  // they are defined if charge and temperature were not defined for
  // that species.
  substitutePermissions("input",
                        {"charge", "temperature", "density", "pressure", "momentum"});
  // FIXME: These will actually only be written if density, pressure,
  // and momentum are set, respectively. They also require charge and
  // temperature to have been set.
  substitutePermissions("output", {"density_source", "energy_source", "momentum_source"});
}

void DiamagneticDrift::coreAverage(Field3D& f) {
  BoutReal local_sum = 0.0;
  for (int jy = mesh->ystart; jy <= mesh->yend; ++jy) {
    BoutReal zavg = 0.0;
    for (int jz = mesh->zstart; jz <= mesh->zend; ++jz) {
      zavg += f(mesh->xstart, jy, jz);
    }
    zavg /= mesh->zend - mesh->zstart + 1;
    local_sum += this->cell_volume(mesh->xstart, jy) * zavg;
  }
  // Sum over processors on core boundary
  BoutReal global_sum{0.0};
  MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, mesh->getYcomm(0));

  const BoutReal average = global_sum / core_ring_volume;
  for (int jy = mesh->ystart; jy <= mesh->yend; ++jy) {
    for (int jz = mesh->zstart; jz <= mesh->zend; ++jz) {
      f(mesh->xstart, jy, jz) = average;
    }
  }
}

void DiamagneticDrift::transform_impl(GuardedOptions& state) {
  // Iterate through all subsections
  GuardedOptions allspecies = state["species"];

  for (auto& kv : allspecies.getChildren()) {
    GuardedOptions species = allspecies[kv.first]; // Note: Need non-const

    if (!(species.isSet("charge") and species.isSet("temperature")))
      continue; // Skip, go to next species

    // Calculate diamagnetic drift velocity for this species
    auto q = get<BoutReal>(species["charge"]);
    if (fabs(q) < 1e-5) {
      continue;
    }
    auto T = GET_VALUE(Field3D, species["temperature"]);

    // Diamagnetic drift velocity
    Vector3D vD = (T / q) * Curlb_B;

    if (IS_SET(species["density"])) {
      auto N = GET_VALUE(Field3D, species["density"]);

      Field3D sink = (divergence_form) ?
                                       // Divergence form: Div(n v_D)
                         FV::Div_f_v(N, vD, bndry_flux)
                                       :
                                       // Gradient form: Curlb_B dot Grad(N T / q)
                         Curlb_B * Grad(N * T / q);

      if (average_core) {
        coreAverage(sink);
      }

      subtract(species["density_source"], sink);
    }

    if (IS_SET(species["pressure"])) {
      auto P = get<Field3D>(species["pressure"]);

      Field3D sink = (divergence_form) ? (5. / 2) * FV::Div_f_v(P, vD, bndry_flux)
                                       : (5. / 2) * Curlb_B * Grad(P * T / q);

      if (average_core) {
        coreAverage(sink);
      }

      subtract(species["energy_source"], sink);
    }

    if (IS_SET(species["momentum"])) {
      auto NV = get<Field3D>(species["momentum"]);
      Field3D sink = (divergence_form) ? FV::Div_f_v(NV, vD, bndry_flux)
                                       : Curlb_B * Grad(NV * T / q);

      if (average_core) {
        coreAverage(sink);
      }

      subtract(species["momentum_source"], sink);
    }
  }
}
