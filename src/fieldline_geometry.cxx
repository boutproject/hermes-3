#include "../include/fieldline_geometry.hxx"
#include <bout/constants.hxx>
#include <bout/coordinates.hxx>
#include <bout/mesh.hxx>
#include <bout/field_factory.hxx>
#include <cmath>

using bout::globals::mesh;

Field3D calculate_Lpar() {
    // Get lpar
    const int MYPE = BoutComm::rank();   // Current rank
    const int NPES = BoutComm::size();    // Number of procs
    const int NYPE = NPES / mesh->NXPE;    // Number of procs in Y
    Coordinates *coord = mesh->getCoordinates();
    Field3D lpar{0.0};

    BoutReal offset = 0;   // Offset to ensure ylow domain boundary starts at 0
    auto dy = coord->dy;
    lpar(0,0,0) = 0.5 * dy(0,0,0);
    for (int id = 0; id <= NYPE-1; ++id) {   // Iterate through each proc
        for (int j = 0; j <= mesh->LocalNy; ++j) {
            if (j > 0) {
                lpar(0, j, 0) = lpar(0, j-1, 0) + 0.5 * dy(0, j-1, 0) + 0.5 * dy(0, j, 0);
            }
        }
        mesh->communicate(lpar);  // Communicate across guard cells so other procs keep counting
        if (MYPE == 0) {
        offset = (lpar(0,1,0) + lpar(0,2,0)) / 2;  // Offset lives on first proc
        }
    }
    MPI_Bcast(&offset, 1, MPI_DOUBLE, 0, BoutComm::get());  // Ensure all procs get offset
    lpar -= offset;

    return lpar;
}

FieldlineGeometry::FieldlineGeometry(std::string, Options& options, Solver*) {
    Options& geo_options = options["fieldline_geometry"];
    const Options& mesh_options = options["mesh"];
    const Options& units = options["units"];
    BoutReal Lnorm = get<BoutReal>(units["meters"]);
    BoutReal Bnorm = get<BoutReal>(units["Tesla"]);

    Coordinates *coord = mesh->getCoordinates();

    lpar = calculate_Lpar() / Lnorm;

    std::string lambda_int_str = geo_options["lambda_int"]
        .doc("Function for the integral heat flux width lambda_int = lambda_q + 1.64 S [m].")
        .as<std::string>();
    
    std::string fieldline_radius_str = geo_options["fieldline_radius"]
        .doc("Function for the fieldline major radius R [m].")
        .as<std::string>();
    
    std::string poloidal_magnetic_field_str = geo_options["poloidal_magnetic_field"]
        .doc("Function for the poloidal magnetic field strength Bpol [T].")
        .as<std::string>();
    
    FieldGeneratorPtr lambda_int_function = FieldFactory::get()->parse(lambda_int_str, &geo_options);
    FieldGeneratorPtr poloidal_magnetic_field_function = FieldFactory::get()->parse(poloidal_magnetic_field_str, &geo_options);
    FieldGeneratorPtr fieldline_radius_function = FieldFactory::get()->parse(fieldline_radius_str, &geo_options);

    lambda_int.allocate();
    fieldline_radius.allocate();
    poloidal_magnetic_field.allocate();
    toroidal_magnetic_field.allocate();
    total_magnetic_field.allocate();
    pitch_angle.allocate();

    BOUT_FOR(i, lpar.getRegion("RGN_ALL")) {
        lambda_int[i] = lambda_int_function->generate(bout::generator::Context().set("lpar", lpar[i] * Lnorm)) / Lnorm;
        fieldline_radius[i] = fieldline_radius_function->generate(bout::generator::Context().set("lpar", lpar[i] * Lnorm)) / Lnorm;
        poloidal_magnetic_field[i] = poloidal_magnetic_field_function->generate(bout::generator::Context().set("lpar", lpar[i] * Lnorm)) / Bnorm;
    }

    bool compute_B_from_R = geo_options["compute_Btor_from_R"]
        .doc("Compute Btor = B_tor,upstream * R_upstream / R if true, or else use a function for toroidal_magnetic_field")
        .withDefault<bool>(true);
    
    if (compute_B_from_R) {
        geo_options["toroidal_magnetic_field"].setConditionallyUsed();

        BoutReal upstream_toroidal_magnetic_field = geo_options["upstream_toroidal_magnetic_field"]
            .doc("Upstream toroidal magnetic field strength Btor,u [T]")
            .as<BoutReal>();
        
        toroidal_magnetic_field = (upstream_toroidal_magnetic_field / Bnorm) * fieldline_radius(0, mesh->ystart, 0) / fieldline_radius;
    } else {
        geo_options["upstream_toroidal_magnetic_field"].setConditionallyUsed();
        
        std::string toroidal_magnetic_field_str = geo_options["toroidal_magnetic_field"]
            .doc("Function for the toroidal magnetic field strength Btor [T].")
            .as<std::string>();
        FieldGeneratorPtr toroidal_magnetic_field_function = FieldFactory::get()->parse(toroidal_magnetic_field_str, &geo_options);

        BOUT_FOR(i, lpar.getRegion("RGN_ALL")) {
            toroidal_magnetic_field[i] = toroidal_magnetic_field_function->generate(bout::generator::Context().set("lpar", lpar[i] * Lnorm)) / Bnorm;
        }
    }
    total_magnetic_field = sqrt(toroidal_magnetic_field*toroidal_magnetic_field + poloidal_magnetic_field*poloidal_magnetic_field);
    pitch_angle = poloidal_magnetic_field / total_magnetic_field;

    transport_broadening = lambda_int / lambda_int(0, mesh->ystart, 0);
    flux_expansion = pitch_angle(0, mesh->ystart, 0) / pitch_angle;

    // Compute the effective magnetic field strength, which is the actual magnetic field strength divided by the transport broadening.
    Field3D effective_magnetic_field_strength = total_magnetic_field / transport_broadening;
    // Bxy is no longer consistent with the Jacobian. Set equal to NaN, to prevent anyone from using it.
    BOUT_FOR(i, coord->Bxy.getRegion("RGN_ALL")) {
        coord->Bxy[i] = std::numeric_limits<BoutReal>::quiet_NaN();
    }
    for (int j = mesh->ystart; j <= mesh->yend; ++j) {
        // N.b. The Jacobian has units of [m / radian T], which is why we need an extra factor of Lnorm.
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

    diagnose = geo_options["diagnose"]
                    .doc("Output additional diagnostics?")
                    .withDefault<bool>(false);
}

void FieldlineGeometry::transform(Options& state) {

}

void FieldlineGeometry::outputVars(Options& state) {
    AUTO_TRACE();
    auto Lnorm = get<BoutReal>(state["rho_s0"]);
    auto Bnorm = get<BoutReal>(state["Bnorm"]);

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

}}