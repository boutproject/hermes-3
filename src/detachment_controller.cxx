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

            if (debug >= 3) {
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

            // Make sure that the linear interpolation returns a point between the two sample
            // points
            ASSERT2(((x1 < distance_from_upstream) && (distance_from_upstream < x2)));
            ASSERT2(((0.0 < distance_from_upstream) && (distance_from_upstream < connection_length)));

            detachment_front_found = true;
            break;
        }
    }
    // Apply corrections in case something funky happened in the linear interpolation
    distance_from_upstream = std::max(distance_from_upstream, 0.0);
    distance_from_upstream = std::min(distance_from_upstream, connection_length);

    detachment_front_location = connection_length - distance_from_upstream;

    // Part 2: compute the response
    
    // Get the time in real units
    time = get<BoutReal>(state["time"]) * time_normalisation;
    // Compute the error
    error = detachment_front_setpoint - detachment_front_location;

    if ((((time - previous_time) > min_time_for_change) and (debug >= 2)) or (debug >= 3)) {
        output << endl;
        output << "detachment_front_location: " << detachment_front_location << endl;
        output << "time:                      " << time << endl;
        output << "time - previous_time       " << (time - previous_time) << endl;
        output << "error - previous_error     " << (error - previous_error) << endl;
        output << "time threshold met?        " << ((time - previous_time) >= min_time_for_change) << endl;
        output << "error threshold met?       " << (fabs(error - previous_error) >= min_error_for_change) << endl;
        output << "passed threshold time?     " << (time > settling_time) << endl;
        output << "reevaluate control?        " << ((time > settling_time) && ((time - previous_time) >= min_time_for_change) && (fabs(error - previous_error) >= min_error_for_change)) << endl;
        output << "control:                   " << control << endl;
        output << endl;
    }

    if ((time > settling_time) && ((time - previous_time) >= min_time_for_change) && (fabs(error - previous_error) >= min_error_for_change)) {
        change_in_time = time - previous_time;
        change_in_error = first_step ? 0.0 : error - previous_error;
            
        derivative = ((1.0 - alpha_de) * change_in_error + alpha_de * previous_change_in_error) / change_in_time;
        change_in_derivative = (1.0 - alpha_d2e) * (derivative - previous_derivative) + alpha_d2e * previous_change_in_derivative;
        
        error_integral = first_step ? 0.0 : error_integral + change_in_time * 0.5 * (error + previous_error);

        if (velocity_form) {
            change_in_control = response_sign * controller_gain * (
                change_in_error
                + (change_in_time / integral_time) * error
                + derivative_time * change_in_derivative
            );

            control = previous_control + change_in_control;
        } else {
            control = control_offset + response_sign * controller_gain * (
                error
                + error_integral / integral_time
                + derivative_time * derivative
            );

            change_in_control = control - previous_control;
        }
        
        control = std::max(control, minval_for_source_multiplier);
        control = std::min(control, maxval_for_source_multiplier);

        if (debug >= 1) {
            output << endl;
            output << "detachment_front_location: " << detachment_front_location << endl;
            output << "previous_time:             " << previous_time << endl;
            output << "change_in_time:            " << change_in_time << endl;
            output << "time:                      " << time << endl;
            output << "previous_error:            " << previous_error << endl;
            output << "change_in_error:           " << change_in_error << endl;
            output << "error:                     " << error << endl;
            output << "error_integral:            " << error_integral << endl;
            output << "previous_derivative:       " << previous_derivative << endl;
            output << "change_in_derivative:      " << change_in_derivative << endl;
            output << "derivative:                " << derivative << endl;
            output << "previous_control:          " << previous_control << endl;
            output << "change_in_control:         " << change_in_control << endl;
            output << "control:                   " << control << endl;
            output << endl;
        }

        previous_time = time;
        previous_error = error;
        
        previous_derivative = derivative;
        previous_control = control;
        previous_change_in_error = change_in_error;
        previous_change_in_derivative = change_in_derivative;
    
        first_step = false;
    }
    ASSERT2(std::isfinite(control));

    // Part 3: Apply the source
    detachment_source_feedback = control * source_shape;
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

        if (control_mode == control_power) {
            add(state["species"][trimmed_species]["energy_source"], scaling_factor * detachment_source_feedback);
        } else if (control_mode == control_particles) {
            add(state["species"][trimmed_species]["density_source"], scaling_factor * detachment_source_feedback);
        } else {
            ASSERT2(false);
        }

        ++species_it;
        ++scaling_factor_it;
    }
}
