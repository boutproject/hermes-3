#include "../include/upstream_density_feedback.hxx"

#include <bout/mesh.hxx>
using bout::globals::mesh;

void UpstreamDensityFeedback::transform(Options& state) {
  Options& species = state["species"][name];

  // Doesn't need all boundaries to be set
  Field3D N = getNoBoundary<Field3D>(species["density"]);

  auto time = get<BoutReal>(state["time"]);

  for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
    int jz = 0;
    error = density_upstream - N(r.ind, mesh->ystart, jz);

    // PI controller, using crude integral of the error
    if (density_error_lasttime < 0.0) {
      // First time this has run
      density_error_lasttime = time;
      density_error_last = error;
    }

    // Integrate using Trapezium rule
    if (time > density_error_lasttime) { // Since time can decrease
      density_error_integral += (time - density_error_lasttime) * 0.5 *
        (error + density_error_last);
    }

    if ((density_error_integral < 0.0) && density_integral_positive) {
      // Limit density_error_integral to be >= 0
      density_error_integral = 0.0;
    }

    // Calculate source from combination of error and integral
    integral_term = density_controller_i * density_error_integral;
    proportional_term = density_controller_p * error;
    source_multiplier = proportional_term + integral_term;

    if ((source_multiplier < 0.0) && density_source_positive) {
      source_multiplier = 0.0; // Don't remove particles
    }

    density_error_last = error;
    density_error_lasttime = time;
    break;
  }

  // The upstream density is only present in processor 0, so for other
  // processors the above loop is skipped. Need to broadcast the values from
  // processor 0 to the other processors
  MPI_Bcast(&source_multiplier, 1, MPI_DOUBLE, 0, BoutComm::get());
  MPI_Bcast(&proportional_term, 1, MPI_DOUBLE, 0, BoutComm::get());
  MPI_Bcast(&integral_term, 1, MPI_DOUBLE, 0, BoutComm::get());
  ASSERT2(std::isfinite(source_multiplier));

  // Scale the source and add to the species density source
  add(species["density_source"], source_multiplier * density_source_shape);

  // Adding particles to a flowing plasma reduces its kinetic energy
  // -> Increase internal energy
  if (IS_SET_NOBOUNDARY(species["velocity"]) and IS_SET(species["AA"])) {
    const Field3D V = GET_NOBOUNDARY(Field3D, species["velocity"]);
    const BoutReal Mi = get<BoutReal>(species["AA"]);
    // Internal energy source
    add(species["energy_source"], 0.5 * Mi * SQ(V) * source_multiplier * density_source_shape);
  }
}
