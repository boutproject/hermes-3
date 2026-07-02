#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/field2d.hxx>
#include <bout/field3d.hxx>
#include <bout/fv_ops.hxx>
#include <bout/globals.hxx>
#include <bout/options.hxx>
#include <bout/output.hxx>
#include <bout/solver.hxx>
#include <bout/sys/range.hxx>
#include <bout/utils.hxx>
#include <bout/vecops.hxx>

#include "../include/component.hxx"
#include "../include/diamagnetic_drift.hxx"
#include "../include/guarded_options.hxx"
#include "../include/permissions.hxx"

#include <string>

using bout::globals::mesh;

DiamagneticDrift::DiamagneticDrift(std::string name, Options& alloptions,
                                   [[maybe_unused]] Solver* solver)
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
  if (mesh->get(Curlb_B, "bxcv") != 0) {
    output_warn.write("Diamagnetic drift: Couldn't read curvature vector 'bxcv'");
    Curlb_B.x = Curlb_B.y = Curlb_B.z = 0.0;
  }

  Options& paralleltransform = Options::root()["mesh"]["paralleltransform"];
  if (paralleltransform.isSet("type")
      and paralleltransform["type"].as<std::string>() == "shifted") {
    Field2D I;
    if (mesh->get(I, "sinty") != 0) {
      I = 0.0;
    }
    Curlb_B.z += I * Curlb_B.x;
  }

  // Normalise

  // Get the units
  const auto& units = alloptions["units"];
  const BoutReal Bnorm = get<BoutReal>(units["Tesla"]);
  const BoutReal Lnorm = get<BoutReal>(units["meters"]);

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

Field3D DiamagneticDrift::calculateDivergenceForm(const Field3D& quantity,
                                                  const Field3D& temperature,
                                                  BoutReal charge, BoutReal factor) {
  if (fabs(charge) < 1e-5) {
    throw BoutException(
        "DiamagneticDrift::calculateDivergenceForm requires a non-zero charge");
  }

  const Vector3D vD = (temperature / charge) * Curlb_B;
  return factor * FV::Div_f_v(quantity, vD, bndry_flux);
}

Field3D DiamagneticDrift::calculateGradientForm(const Field3D& quantity,
                                                const Field3D& temperature,
                                                BoutReal charge, BoutReal factor) {
  if (fabs(charge) < 1e-5) {
    throw BoutException(
        "DiamagneticDrift::calculateGradientForm requires a non-zero charge");
  }

  return factor * Curlb_B * Grad(quantity * temperature / charge);
}

void DiamagneticDrift::addDiamagneticSources(GuardedOptions& species) {
  if (!(species.isSet("charge") and species.isSet("temperature"))) {
    return;
  }

  const auto charge = get<BoutReal>(species["charge"]);
  if (fabs(charge) < 1e-5) {
    return;
  }
  const auto temperature = GET_VALUE(Field3D, species["temperature"]);

  if (IS_SET(species["density"])) {
    Field3D sink = divergence_form
                       ? calculateDivergenceForm(GET_VALUE(Field3D, species["density"]),
                                                 temperature, charge)
                       : calculateGradientForm(GET_VALUE(Field3D, species["density"]),
                                               temperature, charge);

    if (average_core) {
      coreAverage(sink);
    }

    subtract(species["density_source"], sink);
  }

  if (IS_SET(species["pressure"])) {
    Field3D sink = divergence_form
                       ? calculateDivergenceForm(get<Field3D>(species["pressure"]),
                                                 temperature, charge, 5. / 2)
                       : calculateGradientForm(get<Field3D>(species["pressure"]),
                                               temperature, charge, 5. / 2);

    if (average_core) {
      coreAverage(sink);
    }

    subtract(species["energy_source"], sink);
  }

  if (IS_SET(species["momentum"])) {
    Field3D sink = divergence_form
                       ? calculateDivergenceForm(get<Field3D>(species["momentum"]),
                                                 temperature, charge)
                       : calculateGradientForm(get<Field3D>(species["momentum"]),
                                               temperature, charge);

    if (average_core) {
      coreAverage(sink);
    }

    subtract(species["momentum_source"], sink);
  }
}

void DiamagneticDrift::coreAverage(Field3D& f) {
  if (!average_core) {
    throw BoutException("DiamagneticDrift::coreAverage requires average_core=true");
  }

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
    addDiamagneticSources(species);
  }
}
