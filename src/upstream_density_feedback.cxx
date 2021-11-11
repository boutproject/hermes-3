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
    BoutReal error = density_upstream - N(r.ind, mesh->ystart, jz);

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
    source_multiplier = density_controller_p * error +
      density_controller_i * density_error_integral;

    if ((source_multiplier < 0.0) && density_source_positive) {
      source_multiplier = 0.0; // Don't remove particles
    }

    density_error_last = error;
    density_error_lasttime = time;
    break;
  }

  // Broadcast the value of source from processor 0
  MPI_Bcast(&source_multiplier, 1, MPI_DOUBLE, 0, BoutComm::get());
  ASSERT2(std::isfinite(source_multiplier));

  // Scale the source and add to the species density source
  add(species["density_source"], source_multiplier * density_source_shape);
}
