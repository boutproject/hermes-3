#include "../include/detachment_controller.hxx"

#include <bout/mesh.hxx>
using bout::globals::mesh;
#include <algorithm>  // Include for std::max
#include <iostream>
#include <iomanip>

void DetachmentController::transform(Options& state) {

    std::cout << std::fixed << std::setprecision(15);

    Field3D neutral_density = getNoBoundary<Field3D>(state["species"][neutral_species]["density"]);
    Field3D electron_density = getNoBoundary<Field3D>(state["species"]["e"]["density"]);
    Coordinates *coord = mesh->getCoordinates();

    bool first_time = (previous_time < 0.0);
    auto time = get<BoutReal>(state["time"]);
    if (first_time) {
        previous_time = time;
    }
    bool state_changed = (fabs(time - previous_time) > min_time_for_change);

    // Set the initial value so that if no point has Nn > Ne, the detachment front is
    // at the target.
    detachment_front_location = connection_length;

    BoutReal distance_from_upstream = -(coord->dy(0, mesh->ystart - 1, 0) + coord->dy(0, mesh->ystart, 0)) / 4.0;
    // Looking at https://github.com/boutproject/xhermes/blob/main/xhermes/accessors.py#L58
    // dy from the first two cells cancels out
    for (int j = mesh->ystart; j <= mesh->yend; ++j) {

        distance_from_upstream = distance_from_upstream + 0.5 * coord->dy(0, j-1, 0) + 0.5 * coord->dy(0, j, 0);
        BoutReal Nn = neutral_density(0, j, 0);
        BoutReal Ne = electron_density(0, j, 0);

        if (
            Nn > Ne
        ) {
            detachment_front_location = distance_from_upstream;
            break;
        }

    }
    if (error_on_log_distance) {
        // Equivalent to taking log10(detachment_front_desired_location/detachment_front_location)
        error = log10(std::max(connection_length - detachment_front_location, log_floor)) - 
                log10(std::max(connection_length - detachment_front_desired_location, log_floor));
    } else {
        error = detachment_front_desired_location - detachment_front_location;
    }

    if (first_time) {
        previous_error = error;
    }

    bool state_changed = state_changed && (fabs(error - previous_error) > min_error_for_change);

    if (state_changed) {
        std::cout << "Time: " << time << std::endl;
        std::cout << "Time change: " << (time - previous_time) << std::endl;
        std::cout << "Error: " << error << std::endl;
        std::cout << "Error change: " << (error - previous_error) << std::endl;
    }

    // Integrate using Trapezium rule
    // Don't add to the integral when error is larger than integral_threshold. This prevents excessive
    // integral windup.
    if (state_changed && (fabs(error) < integral_threshold)) {
        error_integral = previous_error_integral + (time - previous_time) * 0.5 * (error + previous_error);
        std::cout << "Error integral: " << error_integral << std::endl;
    } else {
        error_integral = previous_error_integral;
    }

    if (fabs(error) < derivative_threshold) {
        error_derivative = 0.0;
    } else if (state_changed) {
        error_derivative = (error - previous_error) / (time - previous_time);
        std::cout << "Error derivative: " << error_derivative << std::endl;
    } else {
        error_derivative = previous_error_derivative;
    }
    

    // Calculate source from combination of error and integral
    if (exponential_control) {
        proportional_term = controller_p * pow(10, error);
        integral_term = controller_i * pow(10, error_integral);
        derivative_term = -1.0 * controller_d * pow(10, error_derivative);
        source_multiplier = proportional_term + integral_term + derivative_term;
        detachment_source_feedback = source_shape * source_multiplier;
    } else {
        proportional_term = controller_p * error;
        integral_term = controller_i * error_integral;
        derivative_term = -1.0 * controller_d * error_derivative;
        source_multiplier = proportional_term + integral_term + derivative_term;
        detachment_source_feedback = source_shape * source_multiplier;
    }

    // Only update the error when the time has changed
    if (state_changed) {
        previous_error = error;
        previous_time = time;
        previous_error_integral = error_integral;
        previous_error_derivative = error_derivative;
    }

    // Need to think about how to do this in parallel runs.

    ASSERT2(std::isfinite(source_multiplier));

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

        if (control_power) {
            add(state["species"][trimmed_species]["energy_source"], scaling_factor * detachment_source_feedback);
        } else {
            add(state["species"][trimmed_species]["density_source"], scaling_factor * detachment_source_feedback);
        }

        ++species_it;
        ++scaling_factor_it;
    }
}
