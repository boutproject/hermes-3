#include "../include/detachment_controller.hxx"

#include <bout/mesh.hxx>
using bout::globals::mesh;

void DetachmentController::transform(Options& state) {

    Field3D neutral_density = getNoBoundary<Field3D>(state["species"][neutral_species]["density"]);
    Field3D electron_density = getNoBoundary<Field3D>(state["species"]["e"]["density"]);
    Coordinates* coord = mesh.getCoordinates();

    auto time = get<BoutReal>(state["time"]);

    // Iterate over all y points, starting at the upstream end of the simulation.
    BoutReal actual_front_location = connection_length;

    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        int jz = 0;

        // Find the point where neutral density exceeds electron density, which is usually a pretty
        // good indicator of the detachment front.
        for (int jy = mesh->ystart; jy <= mesh->yend; jy++) {
            if neutral_density(r.ind, jy, jz) > electron_density(r.ind, jy, jz) {
                actual_front_location = coord->y(r.ind, jy, jz)
                break;
            }
        }

        error = detachment_front_location - actual_front_location

        // PI controller, using crude integral of the error
        if (detachment_control_error_lasttime < 0.0) {
        // First time this has run
        detachment_control_error_lasttime = time;
        detachment_control_error_last = error;
        }

        // Integrate using Trapezium rule
        if (time > detachment_control_error_lasttime) { // Since time can decrease
        detachment_control_error_integral += (time - detachment_control_error_lasttime) * 0.5 *
            (error + detachment_control_error_last);
        }

        if ((detachment_control_error_integral < 0.0) && detachment_control_integral_positive) {
        // Limit detachment_control_error_integral to be >= 0
        detachment_control_error_integral = 0.0;
        }

        // Calculate source from combination of error and integral
        integral_term = detachment_control_controller_i * detachment_control_error_integral;
        proportional_term = detachment_control_controller_p * error;
        source_multiplier = proportional_term + integral_term;

        if ((source_multiplier < 0.0) && detachment_control_source_positive) {
        source_multiplier = 0.0; // Don't remove particles
        }

        detachment_control_error_last = error;
        detachment_control_error_lasttime = time;
        break;

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

        if control_power {
            add(state["species"][trimmed_species]["energy_source"], scaling_factor * source_multiplier * source_shape);
        else {
            add(state["species"][trimmed_species]["density_source"], scaling_factor * source_multiplier * source_shape);
        }

        ++species_it;
        ++scaling_factor_it;
    }
}
