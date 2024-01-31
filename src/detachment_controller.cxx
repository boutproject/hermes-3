#include "../include/detachment_controller.hxx"

#include <bout/mesh.hxx>
using bout::globals::mesh;
#include <algorithm>  // Include for std::max
#include <iostream>
#include <iomanip>

void DetachmentController::transform(Options& state) {

    std::cout << std::fixed << std::setprecision(15);

    // Part 1: compute the detachment front location
    Field3D neutral_density = getNoBoundary<Field3D>(state["species"][neutral_species]["density"]);
    Field3D electron_density = getNoBoundary<Field3D>(state["species"]["e"]["density"]);
    Coordinates *coord = mesh->getCoordinates();

    // Set the initial value so that if no point has Nn > Ne, the detachment front is
    // at the target.

    bool detachment_front_found = false;
    BoutReal distance_from_upstream = -(coord->dy(0, mesh->ystart - 1, 0) + coord->dy(0, mesh->ystart, 0)) / 4.0;
    // Looking at https://github.com/boutproject/xhermes/blob/main/xhermes/accessors.py#L58
    // dy from the first two cells cancels out
    for (int j = mesh->ystart; j <= mesh->yend; ++j) {
        BoutReal x1 = distance_from_upstream; //y position of previous point
        BoutReal a1 = neutral_density(0, j-1, 0); //Nn at previous point
        BoutReal b1 = electron_density(0, j-1, 0); //Ne at previous point

        distance_from_upstream = distance_from_upstream + 0.5 * coord->dy(0, j-1, 0) + 0.5 * coord->dy(0, j, 0);
        
        BoutReal x2 = distance_from_upstream; //y position of current point
        BoutReal a2 = neutral_density(0, j, 0); //Nn at current point
        BoutReal b2 = electron_density(0, j, 0); //Ne at current point
        
        // Find the first point where Nn > Ne, when iterating from upstream to target
        if (a2 > b2) {
            
            // Compute the x-value of the intersection between
            // (x1, a1)->(x2, a2) and (x1, b1)->(x2, b2)
            // https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
            // 
            // This gives a linear approximation of the detachment front position
            distance_from_upstream = ((x1*b2 - b1*x2) - (x1*a2 - a1*x2)) / ((a1 - a2) - (b1 - b2));

            // Make sure that the linear interpolation returns a point between the two sample
            // points
            ASSERT2(((x1 < distance_from_upstream) && (distance_from_upstream < x2)));

            if (debug >= 2) {
                std::cout << std::endl;
                std::cout << "x1: " << x1 << std::endl;
                std::cout << "x2: " << x2 << std::endl;
                std::cout << "xp: " << distance_from_upstream << std::endl;
                std::cout << "a1: " << a1 << std::endl;
                std::cout << "a2: " << a2 << std::endl;
                std::cout << "b1: " << b1 << std::endl;
                std::cout << "b2: " << b2 << std::endl;
                std::cout << "detachment_front_location: " << (connection_length - distance_from_upstream) << std::endl;
                std::cout << std::endl;
            }

            detachment_front_found = true;
            break;
        }
    }
    // Apply corrections in case something funky happened in the linear interpolation
    distance_from_upstream = std::max(distance_from_upstream, 0.0);
    distance_from_upstream = std::min(distance_from_upstream, connection_length);

    detachment_front_location = connection_length - distance_from_upstream;

    // Part 2: compute the response
    if (first_step) {
        previous_control = initial_control;
        first_step = false;
    }
    
    // Get the time in real units
    time = get<BoutReal>(state["time"]) * time_normalisation;
    // Compute the error
    error = (1 - alpha_e) * (detachment_front_setpoint - detachment_front_location) + alpha_e * previous_error;

    if (((time - previous_time) > min_time_for_change) && (fabs(error - previous_error) > min_error_for_change)) {
        change_in_time = time - previous_time;

        change_in_error = (1.0 - alpha_de) * (error - previous_error) + alpha_de * previous_change_in_error;
        derivative = change_in_error / change_in_time;

        change_in_derivative = (1.0 - alpha_d2e) * (derivative - previous_derivative) - alpha_d2e * previous_change_in_derivative;

        change_in_control = response_sign * controller_gain * (
            change_in_error
            + (change_in_time / integral_time) * error
            + (derivative_time / change_in_time) * change_in_derivative
        );

        control = (1.0 - alpha_c) * (previous_control + change_in_control) + alpha_c * previous_change_in_control;
        control = std::max(control, minval_for_source_multiplier);
        control = std::min(control, maxval_for_source_multiplier);

        if (debug >= 1) {
            std::cout << std::endl;
            std::cout << "detachment_front_location: " << detachment_front_location << std::endl;
            std::cout << "previous_time:             " << previous_time << std::endl;
            std::cout << "change_in_time:            " << change_in_time << std::endl;
            std::cout << "time:                      " << time << std::endl;
            std::cout << "previous_error:            " << previous_error << std::endl;
            std::cout << "change_in_error:           " << change_in_error << std::endl;
            std::cout << "error:                     " << error << std::endl;
            std::cout << "previous_derivative:       " << previous_derivative << std::endl;
            std::cout << "change_in_derivative:      " << change_in_derivative << std::endl;
            std::cout << "derivative:                " << derivative << std::endl;
            std::cout << "previous_control:          " << previous_control << std::endl;
            std::cout << "change_in_control:         " << change_in_control << std::endl;
            std::cout << "control:                   " << control << std::endl;
            std::cout << std::endl;
        }

        previous_time = time;
        previous_error = error;
        previous_derivative = derivative;
        previous_control = control;
        previous_change_in_error = change_in_error;
        previous_change_in_derivative = change_in_derivative;
        previous_change_in_control = change_in_control;
            
    } else {
        control = previous_control;
    }

    // Part 3: Apply the source
    detachment_source_feedback = control * source_shape;

    ASSERT2(std::isfinite(control));

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

        add(state["species"][trimmed_species]["energy_source"], scaling_factor * detachment_source_feedback);

        ++species_it;
        ++scaling_factor_it;
    }
}
