#include <bout/bout_types.hxx>
#include <bout/boutcomm.hxx>
#include <bout/boutexception.hxx>
#include <bout/coordinates.hxx>
#include <bout/derivs.hxx>
#include <bout/field.hxx>
#include <bout/field2d.hxx>
#include <bout/field3d.hxx>
#include <bout/fv_ops.hxx>
#include <bout/globals.hxx>
#include <bout/mpi_wrapper.hxx>
#include <bout/options.hxx>
#include <bout/output.hxx>
#include <bout/paralleltransform.hxx>
#include <bout/solver.hxx>
#include <bout/sys/range.hxx>
#include <bout/utils.hxx>
#include <bout/vecops.hxx>

#include "../include/component.hxx"
#include "../include/diamagnetic_drift.hxx"
#include "../include/guarded_options.hxx"
#include "../include/permissions.hxx"

#include <string>
#include <utility>
#include <vector>

using bout::globals::mesh;

DiamagneticDrift::DiamagneticDrift(std::string name, Options& alloptions,
                                   [[maybe_unused]] Solver* solver)
    : NamedComponent(name, {readIfSet("species:{all_species}:{input}"),
                            readWrite("species:{all_species}:{output}")}) {

  // Get options for this component
  auto& options = alloptions[name];

  if (options.isSet("diamag_form")) {
    throw BoutException(
        "Diamagnetic drift: diamag_form has been removed. Replace diamag_form = 1 "
        "with form = divergence, or diamag_form = 0 with form = gradient. "
        "Spatially varying blends were removed because they are not conservative.");
  }

  if (options.isSet("bndry_flux")) {
    if (options.isSet("radial_boundary")) {
      throw BoutException(
          "Diamagnetic drift: bndry_flux and radial_boundary cannot both be set. "
          "Use radial_boundary only.");
    }
    const auto boundary = options["bndry_flux"].as<bool>() ? DiamagneticBoundary::open
                                                           : DiamagneticBoundary::closed;
    options["radial_boundary"].assign(toString(boundary), "bndry_flux");
    if (BoutComm::rank() == 0) {
      output_warn.write("Diamagnetic drift: bndry_flux is deprecated; use "
                        "radial_boundary = {} instead.\n",
                        toString(boundary));
    }
  }

  form = options["form"]
             .doc("Form of the diamagnetic sink: divergence or gradient. "
                  "radial_boundary applies to x boundaries only; sheath (y) "
                  "faces are always closed.")
             .withDefault(DiamagneticForm::divergence);

  radial_boundary = options["radial_boundary"]
                        .doc("Treatment of x boundaries only: open, closed, or "
                             "closed_average (average the innermost core ring). "
                             "Sheath (y) faces are always closed.")
                        .withDefault(DiamagneticBoundary::closed_average);

  if (form == DiamagneticForm::gradient
      and !mesh->getCoordinates()->getParallelTransform().canToFromFieldAligned()) {
    throw BoutException(
        "Diamagnetic drift: gradient form is not supported by parallel transforms "
        "that cannot convert to and from field-aligned coordinates (such as FCI).");
  }

  average_core_active = radial_boundary == DiamagneticBoundary::closed_average
                        and mesh->periodicY(mesh->xstart) and mesh->firstX();

  if (mesh->firstX()) {
    // BoutMesh builds Y communicators per x-processor column and separates
    // core and private-flux ranks. All members must take the same collectives.
    const int local_active = average_core_active ? 1 : 0;
    int min_active{0}, max_active{0};
    const auto comm = mesh->getYcomm(mesh->xstart);
    bout::globals::mpi->MPI_Allreduce(&local_active, &min_active, 1, MPI_INT, MPI_MIN,
                                      comm);
    bout::globals::mpi->MPI_Allreduce(&local_active, &max_active, 1, MPI_INT, MPI_MAX,
                                      comm);
    if (min_active != max_active) {
      throw BoutException(
          "Diamagnetic drift: core averaging must be consistent across all "
          "processors in getYcomm(xstart); average_core_active flags differ.");
    }
  }

  if (BoutComm::rank() == 0) {
    output_info.write("Diamagnetic drift: form = {}, radial_boundary = {}.\n",
                      toString(form), toString(radial_boundary));
    if (radial_boundary == DiamagneticBoundary::closed_average and !average_core_active) {
      output_info.write("Diamagnetic drift: closed_average acts as closed here "
                        "because this processor is not on the core ring.\n");
    }
  }

  if (average_core_active) {
    const auto* coords = mesh->getCoordinates();
    this->cell_volume = coords->dx() * coords->dy() * coords->dz() * coords->J();
    BoutReal local_core_volume = 0.0;
    for (int jy = mesh->ystart; jy <= mesh->yend; ++jy) {
      for (int jz = mesh->zstart; jz <= mesh->zend; ++jz) {
        local_core_volume += this->cell_volume(mesh->xstart, jy, jz);
      }
    }
    // Sum over processors on core boundary
    bout::globals::mpi->MPI_Allreduce(&local_core_volume, &this->core_ring_volume, 1,
                                      MPI_DOUBLE, MPI_SUM, mesh->getYcomm(mesh->xstart));
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
    Coordinates::FieldMetric I;
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

  Curlb_B *= 2. / mesh->getCoordinates()->Bxy();

  // Set drift to zero through sheath boundaries.
  // Flux through those cell faces should be set by sheath.
  for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
    for (int k = 0; k < Curlb_B.y.getNz(); ++k) {
      Curlb_B.y(r.ind, mesh->ystart - 1, k) = -Curlb_B.y(r.ind, mesh->ystart, k);
    }
  }
  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
    for (int k = 0; k < Curlb_B.y.getNz(); ++k) {
      Curlb_B.y(r.ind, mesh->yend + 1, k) = -Curlb_B.y(r.ind, mesh->yend, k);
    }
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

  // Copy in the aligned frame so the sign-flipped curvature cancels the face
  // velocity even when the input guard temperature differs or varies in z.
  Field3D temperature_al = toFieldAligned(copy(temperature), "RGN_ALL");
  for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
    for (int jz = mesh->zstart; jz <= mesh->zend; ++jz) {
      temperature_al(r.ind, mesh->ystart - 1, jz) =
          temperature_al(r.ind, mesh->ystart, jz);
    }
  }
  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
    for (int jz = mesh->zstart; jz <= mesh->zend; ++jz) {
      temperature_al(r.ind, mesh->yend + 1, jz) = temperature_al(r.ind, mesh->yend, jz);
    }
  }
  // Include x and y guards: vD is also used at the radial faces.
  const Field3D sheath_temperature = fromFieldAligned(temperature_al, "RGN_ALL");
  const Vector3D vD = (sheath_temperature / charge) * Curlb_B;
  return factor * FV::Div_f_v(quantity, vD, radial_boundary == DiamagneticBoundary::open);
}

Field3D DiamagneticDrift::calculateGradientForm(const Field3D& quantity,
                                                const Field3D& temperature,
                                                BoutReal charge, BoutReal factor) {
  if (fabs(charge) < 1e-5) {
    throw BoutException(
        "DiamagneticDrift::calculateGradientForm requires a non-zero charge");
  }

  const Field3D g = quantity * temperature / charge;
  Field3D sink =
      factor
      * (Curlb_B.x * DDX(g, CELL_DEFAULT, "C2") + Curlb_B.y * DDY(g, CELL_DEFAULT, "C2")
         + Curlb_B.z * DDZ(g, CELL_DEFAULT, "C2"));

  if (radial_boundary != DiamagneticBoundary::open) {
    const auto* coords = quantity.getCoordinates();
    const auto& J = coords->J();
    const auto& dx = coords->dx();
    const Coordinates::FieldMetric a = J * Curlb_B.x;

    // C2 has the implied face flux F = (a_left*g_right + a_right*g_left)/2.
    // Remove the boundary faces while retaining the cell-centred residual.
    if (mesh->firstX()) {
      const int jx = mesh->xstart;
      for (int jy = mesh->ystart; jy <= mesh->yend; ++jy) {
        for (int jz = mesh->zstart; jz <= mesh->zend; ++jz) {
          const BoutReal flux = 0.5
                                * ((a(jx - 1, jy, jz) * g(jx, jy, jz))
                                   + (a(jx, jy, jz) * g(jx - 1, jy, jz)));
          sink(jx, jy, jz) += factor * flux / (J(jx, jy, jz) * dx(jx, jy, jz));
        }
      }
    }
    if (mesh->lastX()) {
      const int jx = mesh->xend;
      for (int jy = mesh->ystart; jy <= mesh->yend; ++jy) {
        for (int jz = mesh->zstart; jz <= mesh->zend; ++jz) {
          const BoutReal flux = 0.5
                                * ((a(jx, jy, jz) * g(jx + 1, jy, jz))
                                   + (a(jx + 1, jy, jz) * g(jx, jy, jz)));
          sink(jx, jy, jz) -= factor * flux / (J(jx, jy, jz) * dx(jx, jy, jz));
        }
      }
    }
  }

  const auto* coords = quantity.getCoordinates();
  const Field3D g_al = toFieldAligned(g, "RGN_NOX");
  const auto a_al =
      toFieldAligned(Coordinates::FieldMetric{coords->J() * Curlb_B.y}, "RGN_NOX");
  const auto Jdy_al =
      toFieldAligned(Coordinates::FieldMetric{coords->J() * coords->dy()}, "RGN_NOX");
  Field3D correction = zeroFrom(g_al);

  // The gradient form needs the unflipped continuation of J*C.y in each
  // target guard. Extrapolate from interior cells rather than reading the
  // curvature guard used to close the divergence-form velocity.
  for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
    const int jy = mesh->ystart;
    for (int jz = mesh->zstart; jz <= mesh->zend; ++jz) {
      const BoutReal a_guard = (2.0 * a_al(r.ind, jy, jz)) - a_al(r.ind, jy + 1, jz);
      const BoutReal flux = 0.5
                            * ((a_guard * g_al(r.ind, jy, jz))
                               + (a_al(r.ind, jy, jz) * g_al(r.ind, jy - 1, jz)));
      correction(r.ind, jy, jz) += factor * flux / Jdy_al(r.ind, jy, jz);
    }
  }
  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
    const int jy = mesh->yend;
    for (int jz = mesh->zstart; jz <= mesh->zend; ++jz) {
      const BoutReal a_guard = (2.0 * a_al(r.ind, jy, jz)) - a_al(r.ind, jy - 1, jz);
      const BoutReal flux = 0.5
                            * ((a_al(r.ind, jy, jz) * g_al(r.ind, jy + 1, jz))
                               + (a_guard * g_al(r.ind, jy, jz)));
      correction(r.ind, jy, jz) -= factor * flux / Jdy_al(r.ind, jy, jz);
    }
  }
  sink += fromFieldAligned(correction, "RGN_NOBNDRY");

  return sink;
}

BoutReal DiamagneticDrift::coreRingSum(const Field3D& f) const {
  BoutReal local_sum = 0.0;
  for (int jy = mesh->ystart; jy <= mesh->yend; ++jy) {
    for (int jz = mesh->zstart; jz <= mesh->zend; ++jz) {
      local_sum += this->cell_volume(mesh->xstart, jy, jz) * f(mesh->xstart, jy, jz);
    }
  }
  return local_sum;
}

void DiamagneticDrift::setCoreAverage(Field3D& f, BoutReal global_sum) const {
  const BoutReal average = global_sum / core_ring_volume;
  for (int jy = mesh->ystart; jy <= mesh->yend; ++jy) {
    for (int jz = mesh->zstart; jz <= mesh->zend; ++jz) {
      f(mesh->xstart, jy, jz) = average;
    }
  }
}

void DiamagneticDrift::transform_impl(GuardedOptions& state) {
  struct Source {
    GuardedOptions output;
    Field3D sink;
  };
  std::vector<Source> sources;
  std::vector<BoutReal> ring_sums;
  GuardedOptions allspecies = state["species"];

  // First compute every sink and collect its local core-ring sum in the
  // same species/quantity order on all members of the Y communicator.
  for (auto& kv : allspecies.getChildren()) {
    GuardedOptions species = allspecies[kv.first];
    if (!(species.isSet("charge") and species.isSet("temperature"))) {
      continue;
    }
    const auto charge = get<BoutReal>(species["charge"]);
    if (fabs(charge) < 1e-5) {
      continue;
    }
    const auto temperature = GET_VALUE(Field3D, species["temperature"]);

    const auto collect = [&](const char* input, const char* output, BoutReal factor) {
      if (!IS_SET(species[input])) {
        return;
      }
      const auto quantity = GET_VALUE(Field3D, species[input]);
      Field3D sink = form == DiamagneticForm::divergence
                         ? calculateDivergenceForm(quantity, temperature, charge, factor)
                         : calculateGradientForm(quantity, temperature, charge, factor);
      if (average_core_active) {
        ring_sums.push_back(coreRingSum(sink));
      }
      sources.push_back({species[output], std::move(sink)});
    };
    collect("density", "density_source", 1.0);
    collect("pressure", "energy_source", 5. / 2);
    collect("momentum", "momentum_source", 1.0);
  }

  if (average_core_active) {
    // Even an empty batch participates, keeping one collective per transform.
    BoutReal empty_buffer{0.0};
    bout::globals::mpi->MPI_Allreduce(
        MPI_IN_PLACE, ring_sums.empty() ? &empty_buffer : ring_sums.data(),
        static_cast<int>(ring_sums.size()), MPI_DOUBLE, MPI_SUM,
        mesh->getYcomm(mesh->xstart));
  }

  // Apply the global means and write the source terms in a second pass.
  for (std::size_t i = 0; i < sources.size(); ++i) {
    if (average_core_active) {
      setCoreAverage(sources[i].sink, ring_sums[i]);
    }
    subtract(sources[i].output.getWritable(), sources[i].sink);
  }
}
