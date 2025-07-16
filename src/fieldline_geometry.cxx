#include "../include/fieldline_geometry.hxx"
#include <bout/constants.hxx>
#include <bout/coordinates.hxx>
#include <bout/mesh.hxx>
#include <bout/field_factory.hxx>
#include <cmath>

using bout::globals::mesh;

/**
 * @brief Calculates the parallel length lpar along the field line.
 * @details
 * The parallel length is calculated by integrating along the poloidal direction.
 * This function calculates the parallel distance from the upstream boundary
 * for each cell in the simulation domain.
 * @return Field3D containing the parallel length.
 */
Field3D calculate_Lpar() {
    // Get lpar
    const int MYPE = BoutComm::rank();   // Current rank of the processor
    const int NPES = BoutComm::size();    // Total number of processors
    const int NYPE = NPES / mesh->NXPE;    // Number of processors in the Y direction
    Coordinates *coord = mesh->getCoordinates(); // Get the coordinates object from the mesh
    Field3D lpar{0.0}; // Initialize the parallel length field with 0.0

    BoutReal offset = 0;   // Offset to ensure ylow domain boundary starts at 0. This is needed
                           // because the integration starts at the first cell center, not the boundary.
    auto dy = coord->dy;  // Get the poloidal cell width
    lpar(0,0,0) = 0.5 * dy(0,0,0);  // Initialize lpar at the first cell center
                                   // The parallel length at the first cell center is half of the cell width.
    for (int id = 0; id <= NYPE-1; ++id) {   // Iterate through each processor in the poloidal direction
        for (int j = 0; j <= mesh->LocalNy; ++j) {  // Iterate through each cell in the local poloidal domain
            if (j > 0) {
                // Calculate lpar at the current cell center by adding the average of the current
                // and previous cell widths to the previous lpar value. This is a simple trapezoidal rule.
                lpar(0, j, 0) = lpar(0, j-1, 0) + 0.5 * dy(0, j-1, 0) + 0.5 * dy(0, j, 0);
            }
        }
        mesh->communicate(lpar);  // Communicate lpar values across guard cells.
                                 // This is necessary to ensure that the lpar calculation
                                 // is consistent across processors, as the integration is
                                 // performed in a loop over processors.
        if (MYPE == 0) {
            offset = (lpar(0,1,0) + lpar(0,2,0)) / 2;  // Calculate the offset at the first processor (MYPE == 0).
                                                     // The offset is the average of lpar at the second and third cell centers.
        }
    }
    MPI_Bcast(&offset, 1, MPI_DOUBLE, 0, BoutComm::get());  // Broadcast the offset to all processors.
                                                         // This ensures that all processors have the same offset value.
    lpar -= offset;  // Subtract the offset from lpar. This ensures that lpar is 0 at the ylow boundary.

    return lpar;  // Return the calculated parallel length field.
}

/**
 * @brief Constructor for the FieldlineGeometry component.
 * @param name The name of the component.
 * @param options Options for the component.
 * @param solver Pointer to the solver.
 * @details
 * This constructor initializes the FieldlineGeometry component.
 * It calculates and stores various geometric quantities such as the parallel length,
 * magnetic field components, flux expansion, and cell dimensions.
 */
FieldlineGeometry::FieldlineGeometry(std::string, Options& options, Solver*) {
    Options& geo_options = options["fieldline_geometry"];  // Get options specific to fieldline_geometry
    const Options& mesh_options = options["mesh"];        // Get mesh options
    const Options& units = options["units"];              // Get unit options
    BoutReal Lnorm = get<BoutReal>(units["meters"]);     // Get the length normalization factor from the units options
    BoutReal Bnorm = get<BoutReal>(units["Tesla"]);      // Get the magnetic field normalization factor from the units options

    Coordinates *coord = mesh->getCoordinates();          // Get the coordinates object from the mesh
    lpar = calculate_Lpar() / Lnorm;                     // Calculate the parallel length and normalize it

    // Get the string expressions for lambda_int, fieldline_radius and poloidal_magnetic_field from the options
    std::string lambda_int_str = geo_options["lambda_int"]
        .doc("Function for the integral heat flux width lambda_int = lambda_q + 1.64 S [m].")
        .as<std::string>();
    std::string fieldline_radius_str = geo_options["fieldline_radius"]
        .doc("Function for the fieldline major radius R [m].")
        .as<std::string>();
    std::string poloidal_magnetic_field_str = geo_options["poloidal_magnetic_field"]
        .doc("Function for the poloidal magnetic field strength Bpol [T].")
        .as<std::string>();

    // Create FieldGenerator objects for lambda_int, fieldline_radius and poloidal_magnetic_field
    // These FieldGenerators will be used to evaluate the string expressions
    FieldGeneratorPtr lambda_int_function = FieldFactory::get()->parse(lambda_int_str, &geo_options);
    FieldGeneratorPtr poloidal_magnetic_field_function = FieldFactory::get()->parse(poloidal_magnetic_field_str, &geo_options);
    FieldGeneratorPtr fieldline_radius_function = FieldFactory::get()->parse(fieldline_radius_str, &geo_options);

    // Allocate memory for the fields
    lambda_int.allocate();
    fieldline_radius.allocate();
    poloidal_magnetic_field.allocate();
    toroidal_magnetic_field.allocate();
    total_magnetic_field.allocate();
    pitch_angle.allocate();

    // Generate the field data for lambda_int, fieldline_radius and poloidal_magnetic_field
    // using the FieldGenerator objects and normalize the fields
    BOUT_FOR(i, lpar.getRegion("RGN_ALL")) {
        lambda_int[i] = lambda_int_function->generate(bout::generator::Context().set("lpar", lpar[i] * Lnorm)) / Lnorm;
        fieldline_radius[i] = fieldline_radius_function->generate(bout::generator::Context().set("lpar", lpar[i] * Lnorm)) / Lnorm;
        poloidal_magnetic_field[i] = poloidal_magnetic_field_function->generate(bout::generator::Context().set("lpar", lpar[i] * Lnorm)) / Bnorm;
    }

    // Determine whether to compute Btor from R or use a provided function
    bool compute_B_from_R = geo_options["compute_Btor_from_R"]
        .doc("Compute Btor = B_tor,upstream * R_upstream / R if true, or else use a function for toroidal_magnetic_field")
        .withDefault<bool>(true);

    // If compute_B_from_R is true, compute Btor from R
    if (compute_B_from_R) {
        geo_options["toroidal_magnetic_field"].setConditionallyUsed();  // Ensure toroidal_magnetic_field is not used
        BoutReal upstream_toroidal_magnetic_field = geo_options["upstream_toroidal_magnetic_field"]
            .doc("Upstream toroidal magnetic field strength Btor,u [T]")
            .as<BoutReal>();
        // Compute toroidal_magnetic_field using the formula: Btor = B_tor,upstream * R_upstream / R
        toroidal_magnetic_field = (upstream_toroidal_magnetic_field / Bnorm) * fieldline_radius(0, mesh->ystart, 0) / fieldline_radius;
    } else {
        // Otherwise, use the provided function for toroidal_magnetic_field
        geo_options["upstream_toroidal_magnetic_field"].setConditionallyUsed(); // Ensure upstream_toroidal_magnetic_field is not used
        std::string toroidal_magnetic_field_str = geo_options["toroidal_magnetic_field"]
            .doc("Function for the toroidal magnetic field strength Btor [T].")
            .as<std::string>();
        FieldGeneratorPtr toroidal_magnetic_field_function = FieldFactory::get()->parse(toroidal_magnetic_field_str, &geo_options);

        // Generate the field data for toroidal_magnetic_field using the FieldGenerator object
        // and normalize the field
        BOUT_FOR(i, lpar.getRegion("RGN_ALL")) {
            toroidal_magnetic_field[i] = toroidal_magnetic_field_function->generate(bout::generator::Context().set("lpar", lpar[i] * Lnorm)) / Bnorm;
        }
    }

    // Compute the total magnetic field strength
    total_magnetic_field = sqrt(toroidal_magnetic_field*toroidal_magnetic_field + poloidal_magnetic_field*poloidal_magnetic_field);
    // Compute the pitch angle
    pitch_angle = poloidal_magnetic_field / total_magnetic_field;
    // Compute the flux tube broadening factor due to cross-field transport
    transport_broadening = lambda_int / lambda_int(0, mesh->ystart, 0);
    // Compute the flux expansion factor
    flux_expansion = pitch_angle(0, mesh->ystart, 0) / pitch_angle;

    // Compute the effective magnetic field strength, which is the actual magnetic field strength divided by the transport broadening.
    // This is used in the Jacobian calculation.
    Field3D effective_magnetic_field_strength = total_magnetic_field / transport_broadening;

    // Bxy is no longer consistent with the Jacobian.
    // Set equal to NaN, to prevent anyone from using it.
    BOUT_FOR(i, coord->Bxy.getRegion("RGN_ALL")) {
        coord->Bxy[i] = std::numeric_limits<BoutReal>::quiet_NaN();
    }

    // Calculate the Jacobian
    for (int j = mesh->ystart; j <= mesh->yend; ++j) {
        // N.b.
        // The Jacobian has units of [m / radian T], which is why we need an extra factor of Lnorm.
        coord->J(0, j) = 1 / effective_magnetic_field_strength(0, j, 0) / Lnorm;
    }

    // Parallel length of cell
    Field3D dlpar = coord->dy / Lnorm;
    // Width of flux tube in the radial direction
    flux_tube_width = lambda_int * flux_expansion;
    // Length of the cell in the poloidal direction
    cell_poloidal_length = dlpar * pitch_angle;
    // Poloidal area of the cell (poloidal length times circumference)
    cell_side_area = cell_poloidal_length * 2.0 * PI * fieldline_radius;
    // Volume of the cell (poloidal area times flux tube width)
    cell_volume = cell_side_area * flux_tube_width;

    // Determine whether to output additional diagnostics
    diagnose = geo_options["diagnose"]
                    .doc("Output additional diagnostics?")
                    .withDefault<bool>(false);
}

void FieldlineGeometry::transform(Options& state) {
    // This method is intentionally left empty.
    // If you want the geometry to evolve during the simulation (for instance, increasing
    // the cross-field broadening based on divertor conditions) you can implement it here.
}

/**
 * @brief Output variables to the state.
 * @param state Options containing the state.
 * @details
 * This method adds the calculated geometric variables to the output state,
 * so that they can be used by other components or output to files.
 */
void FieldlineGeometry::outputVars(Options& state) {
    AUTO_TRACE();  // Automatically add tracing information for debugging purposes
    auto Lnorm = get<BoutReal>(state["rho_s0"]);  // Get the length normalization factor from the state
    auto Bnorm = get<BoutReal>(state["Bnorm"]);   // Get the magnetic field normalization factor from the state

    // Only output the variables if the diagnose option is true
    if (diagnose) {

        set_with_attrs(
            state[std::string("fieldline_geometry_lpar")], lpar,
            {
                {"units", "m"},
                {"conversion", Lnorm},
                {"long_name", "Parallel distance from upstream"},
                {"source", "fieldline_geometry"}
            }
        );

        set_with_attrs(
            state[std::string("fieldline_geometry_lambda_int")], lambda_int,
            {
                {"units", "m"},
                {"conversion", Lnorm},
                {"long_name", "Flux tube radial width mapped upstream (lambda_q + 1.64S)"},
                {"source", "fieldline_geometry"}
            }
        );

        set_with_attrs(
            state[std::string("fieldline_geometry_pitch_angle")], pitch_angle,
            {
                {"units", ""},
                {"long_name", "sin(theta) = Bpol/B"},
                {"source", "fieldline_geometry"}
            }
        );

        set_with_attrs(
            state[std::string("fieldline_geometry_fieldline_radius")], fieldline_radius,
            {
                {"units", "m"},
                {"conversion", Lnorm},
                {"long_name", "Major radius of fieldline"},
                {"source", "fieldline_geometry"}
            }
        );

        set_with_attrs(
            state[std::string("fieldline_geometry_poloidal_magnetic_field")], poloidal_magnetic_field,
            {
                {"units", "T"},
                {"conversion", Bnorm},
                {"long_name", "Poloidal magnetic field strength along fieldline"},
                {"source", "fieldline_geometry"}
            }
        );

        set_with_attrs(
            state[std::string("fieldline_geometry_toroidal_magnetic_field")], toroidal_magnetic_field,
            {
                {"units", "T"},
                {"conversion", Bnorm},
                {"long_name", "Toroidal magnetic field strength along fieldline"},
                {"source", "fieldline_geometry"}
            }
        );

        set_with_attrs(
            state[std::string("fieldline_geometry_total_magnetic_field")], total_magnetic_field,
            {
                {"units", "T"},
                {"conversion", Bnorm},
                {"long_name", "Magnetic field strength along fieldline"},
                {"source", "fieldline_geometry"}
            }
        );

        set_with_attrs(
            state[std::string("fieldline_geometry_transport_broadening")], transport_broadening,
            {
                {"units", ""},
                {"long_name", "Flux tube broadening factor due to cross-field transport"},
                {"source", "fieldline_geometry"}
            }
        );

        set_with_attrs(
            state[std::string("fieldline_geometry_flux_expansion")], flux_expansion,
            {
                {"units", ""},
                {"long_name", "Flux expansion"},
                {"source", "fieldline_geometry"}
            }
        );

        set_with_attrs(
            state[std::string("fieldline_geometry_flux_tube_width")], flux_tube_width,
            {
                {"units", "m"},
                {"conversion", Lnorm},
                {"long_name", "Local flux tube radial width (lambda_q + 1.64S) * flux_expansion"},
                {"source", "fieldline_geometry"}
            }
        );

        set_with_attrs(
            state[std::string("fieldline_geometry_cell_poloidal_length")], cell_poloidal_length,
            {
                {"units", "m"},
                {"conversion", Lnorm},
                {"long_name", "Poloidal length of cell (dl * Bpol / B)"},
                {"source", "fieldline_geometry"}
            }
        );

        set_with_attrs(
            state[std::string("fieldline_geometry_cell_side_area")], cell_side_area,
            {
                {"units", "m^2"},
                {"conversion", Lnorm*Lnorm},
                {"long_name", "Poloidal length of cell times circumference"},
                {"source", "fieldline_geometry"}
            }
        );

        set_with_attrs(
            state[std::string("fieldline_geometry_cell_volume")], cell_volume,
            {
                {"units", "m^3"},
                {"conversion", Lnorm*Lnorm*Lnorm},
                {"long_name", "Poloidal length of cell times circumference times flux tube radial width"},
                {"source", "fieldline_geometry"}
            }
        );
    }
}