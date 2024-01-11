#include "../include/temperature_feedback.hxx"

#include <bout/mesh.hxx>
using bout::globals::mesh;

void TemperatureFeedback::transform(Options& state) {
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

  // std::list<std::string>::iterator species_it = species_list.begin();
  // std::list<std::string>::iterator scaling_factor_it = scaling_factors_list.begin();

  // while (species_it != species_list.end() && scaling_factor_it != scaling_factors_list.end()) {

  //     std::string sourced_species = trim(*species_it, " \t\r()"); // The species name in the list
  //     BoutReal scaling_factor = stringToReal(trim(*scaling_factor_it, " \t\r()"));

  //     if (sourced_species.empty())
  //       continue; // Missing
      
  //     add(state["species"][sourced_species]["energy_source"], scaling_factor * source_multiplier * source_shape);

  //     ++species_it;
  //     ++scaling_factor_it;
  // }

  auto species_it = species_list.begin();
  auto scaling_factor_it = scaling_factors_list.begin();

  while (species_it != species_list.end() && scaling_factor_it != scaling_factors_list.end()) {
      std::string trimmed_species = trim(*species_it);
      std::string trimmed_scaling_factor = trim(*scaling_factor_it);

      if (trimmed_species.empty() || trimmed_scaling_factor.empty()) {
          ++species_it;
          ++scaling_factor_it;
          continue; // Skip this iteration if either trimmed string is empty
      }

      BoutReal scaling_factor = stringToReal(trimmed_scaling_factor);
      add(state["species"][trimmed_species]["energy_source"], scaling_factor * source_multiplier * source_shape);

      ++species_it;
      ++scaling_factor_it;
  }

  // Scale the source and add to the species temperature source
  // add(species["energy_source"], source_multiplier * source_shape);
  // add(state["species"]["d+"]["energy_source"], source_multiplier * source_shape);
}
