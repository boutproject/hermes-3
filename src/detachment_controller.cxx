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
        distance_from_upstream = distance_from_upstream + 0.5 * coord->dy(0, j-1, 0) + 0.5 * coord->dy(0, j, 0);
        if (neutral_density(0, j, 0) > electron_density(0, j, 0)) {
            detachment_front_found = true;
            break;
        }
    }
    // Part 2: compute the response
    if (first_step) {
        previous_control = initial_control;
        first_step = false;
    }
    detachment_front_location = connection_length - distance_from_upstream;
        
    BoutReal time = get<BoutReal>(state["time"]);
    error = detachment_front_setpoint - detachment_front_location;

    change_in_time = time - previous_time;
    change_in_error = error - previous_error;
    derivative = change_in_error / change_in_time;
    change_in_derivative = derivative - previous_derivative;

    change_in_control = response_sign * controller_gain * (
        change_in_error
        + (change_in_time / integral_time) * error
        + (derivative_time / change_in_time) * change_in_derivative
    );
    
    if (debug >= 2) {
        std::cout << std::endl;
        std::cout << "detachment_front_location: " << detachment_front_location << std::endl;
        std::cout << "time: " << time << std::endl;
        std::cout << "error: " << error << std::endl;
        std::cout << "change_in_time: " << change_in_time << std::endl;
        std::cout << "change_in_error: " << change_in_error << std::endl;
        std::cout << "control: " << error << std::endl;
        std::cout << std::endl;
    }

    if ((change_in_time > min_time_for_change) && (fabs(change_in_error) > min_error_for_change)) {
        control = previous_control + change_in_control;
        control = std::max(control, minval_for_source_multiplier);
        control = std::min(control, maxval_for_source_multiplier);

        previous_time = time;
        previous_error = error;
        previous_derivative = derivative;
        previous_control = control;

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
