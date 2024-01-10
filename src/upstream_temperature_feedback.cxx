#include "../include/upstream_temperature_feedback.hxx"

#include <bout/mesh.hxx>
using bout::globals::mesh;

void UpstreamTemperatureFeedback::transform(Options& state) {
  Options& species = state["species"][name];

  // Doesn't need all boundaries to be set
  Field3D T = getNoBoundary<Field3D>(species["temperature"]);

  auto time = get<BoutReal>(state["time"]);

  for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
    int jz = 0;

    if (control_target_temperature) {
      error = temperature_setpoint - T(r.ind, mesh->yend, jz);
    } else {
      error = temperature_setpoint - T(r.ind, mesh->ystart, jz);
    }

    // PI controller, using crude integral of the error
    if (temperature_error_lasttime < 0.0) {
      // First time this has run
      temperature_error_lasttime = time;
      temperature_error_last = error;
    }

    // Integrate using Trapezium rule
    if (time > temperature_error_lasttime) { // Since time can decrease
      temperature_error_integral += (time - temperature_error_lasttime) * 0.5 *
        (error + temperature_error_last);
    }

    if ((temperature_error_integral < 0.0) && temperature_integral_positive) {
      // Limit temperature_error_integral to be >= 0
      temperature_error_integral = 0.0;
    }

    // Calculate source from combination of error and integral
    integral_term = temperature_controller_i * temperature_error_integral;
    proportional_term = temperature_controller_p * error;
    source_multiplier = proportional_term + integral_term;

    if ((source_multiplier < 0.0) && temperature_source_positive) {
      source_multiplier = 0.0; // Don't remove particles
    }

    temperature_error_last = error;
    temperature_error_lasttime = time;
    break;
  }

  // The upstream temperature is only present in processor 0, so for other
  // processors the above loop is skipped. Need to broadcast the values from
  // processor 0 to the other processors
  MPI_Bcast(&source_multiplier, 1, MPI_DOUBLE, 0, BoutComm::get());
  MPI_Bcast(&proportional_term, 1, MPI_DOUBLE, 0, BoutComm::get());
  MPI_Bcast(&integral_term, 1, MPI_DOUBLE, 0, BoutComm::get());
  ASSERT2(std::isfinite(source_multiplier));

  // Scale the source and add to the species temperature source
  add(species["pressure_source"], source_multiplier * pressure_source_shape);
}
