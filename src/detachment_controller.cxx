#include "../include/detachment_controller.hxx"

#include <bout/mesh.hxx>
using bout::globals::mesh;

void DetachmentController::transform(Options& state) {

    Field3D neutral_density = getNoBoundary<Field3D>(state["species"][neutral_species]["density"]);
    Field3D electron_density = getNoBoundary<Field3D>(state["species"]["e"]["density"]);
    Coordinates *coord = mesh->getCoordinates();

    auto time = get<BoutReal>(state["time"]);

    BoutReal closest = INFINITY;
    
    // Set the initial value so that if no point has Nn > Ne, the detachment front is
    // at the target.
    detachment_front_index = mesh->yend;

    // Iterate over all y points
    BOUT_FOR_SERIAL(i, electron_density.getRegion("RGN_NOBNDRY")) {
        BoutReal N_neutral = neutral_density[i];
        BoutReal N_electron = electron_density[i];
        BoutReal test = fabs(N_neutral - N_electron);

        if (
             ((test < closest) && (N_neutral > N_electron))
        ) {
            closest = test;
            detachment_front_index = i.y();
        }
    }

    detachment_front_location = - (coord->dy(0, mesh->ystart - 1, 0) + coord->dy(0, mesh->ystart, 0)) / 4.0;
    // Looking at https://github.com/boutproject/xhermes/blob/main/xhermes/accessors.py#L58
    // dy from the first two cells cancels out
    for (int j = mesh->ystart; j <= detachment_front_index; ++j) {
        detachment_front_location = detachment_front_location + 0.5 * coord->dy(0, j-1, 0) + 0.5 * coord->dy(0, j, 0);
    }

    // Control on the sqrt of the distance
    // if (detachment_front_desired_location >= detachment_front_location) {
    //     error = sqrtf(detachment_front_desired_location - detachment_front_location);
    // } else {
    //     error = -1.0 * sqrtf(detachment_front_location - detachment_front_desired_location);
    // }
    error = detachment_front_desired_location - detachment_front_location;

    // PI controller, using crude integral of the error
    if (error_lasttime < 0.0) {
    // First time this has run
    error_lasttime = time;
    error_last = error;
    }

    // Integrate using Trapezium rule
    if (time > error_lasttime) { // Since time can decrease
    error_integral += (time - error_lasttime) * 0.5 *
        (error + error_last);
    }

    if ((error_integral < 0.0) && force_integral_positive) {
    // Limit error_integral to be >= 0
    error_integral = 0.0;
    }

    // Calculate source from combination of error and integral
    integral_term = controller_i * error_integral;
    proportional_term = controller_p * error;
    source_multiplier = proportional_term + integral_term;

    if (force_source_positive) {
        if ((source_multiplier > 0.0) && (control_power)) {
            source_multiplier = 0.0;
        }
        if ((source_multiplier < 0.0) && (not control_power)) {
            source_multiplier = 0.0;
        }
    }

    error_last = error;
    error_lasttime = time;

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
            add(state["species"][trimmed_species]["energy_source"], -1.0 * scaling_factor * source_multiplier * source_shape);
        } else {
            add(state["species"][trimmed_species]["density_source"], scaling_factor * source_multiplier * source_shape);
        }

        ++species_it;
        ++scaling_factor_it;
    }
}
