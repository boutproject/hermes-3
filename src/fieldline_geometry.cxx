#include "../include/fieldline_geometry.hxx"
#include <bout/constants.hxx>
#include <bout/coordinates.hxx>
#include <bout/mesh.hxx>
#include <bout/field_factory.hxx>
#include <iostream> // For outputting to log file
#include <cmath>

using bout::globals::mesh;

Field3D calculateLpar() {
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
    Options& mesh_options = options["mesh"];
    Coordinates *coord = mesh->getCoordinates();

    const auto& units = options["units"];
    Lnorm = get<BoutReal>(units["meters"]);

    lpar = calculateLpar();

    auto geometric_broadening_str = geo_options["geometric_broadening"]
      .doc("Function for R(lpar) / R(lpar=0) (for lpar in [m]).")
      .as<std::string>();
    FieldGeneratorPtr geometric_broadening_function = FieldFactory::get()->parse(geometric_broadening_str, &geo_options);

    auto transport_broadening_factor_str = geo_options["transport_broadening"]
      .doc("Function for lambda_INT(lpar) / lambda_q(lpar=0) (for lpar in [m]).")
      .as<std::string>();
    FieldGeneratorPtr transport_broadening_factor_function = FieldFactory::get()->parse(transport_broadening_factor_str, &geo_options);

    auto flux_expansion_str = geo_options["flux_expansion"]
      .doc("Function for (B_pol/B_tot(lpar=0)) / (B_pol/B_tot(lpar)) (for lpar in [m]).")
      .as<std::string>();
    FieldGeneratorPtr flux_expansion_function = FieldFactory::get()->parse(flux_expansion_str, &geo_options);

    geometric_broadening.allocate();
    transport_broadening.allocate();
    flux_expansion.allocate();

    BOUT_FOR(i, lpar.getRegion("RGN_ALL")) {
        geometric_broadening[i] = geometric_broadening_function->generate(bout::generator::Context().set("lpar", lpar[i]));
        transport_broadening[i] = transport_broadening_factor_function->generate(bout::generator::Context().set("lpar", lpar[i]));
        flux_expansion[i] = flux_expansion_function->generate(bout::generator::Context().set("lpar", lpar[i]));
    }

    // Make sure that all the broadening factors are normalized to their upstream values
    geometric_broadening /= geometric_broadening(0, mesh->ystart, 0);
    transport_broadening /= transport_broadening(0, mesh->ystart, 0);
    flux_expansion /= flux_expansion(0, mesh->ystart, 0);

    flux_tube_broadening = geometric_broadening * transport_broadening * flux_expansion;
    // Bxy is no longer consistent with the Jacobian. Set equal to NaN, to prevent anyone from using it.
    BOUT_FOR(i, coord->Bxy.getRegion("RGN_ALL")) {
        coord->Bxy[i] = std::numeric_limits<BoutReal>::quiet_NaN();
    }
    for (int j = mesh->ystart; j <= mesh->yend; ++j) {
        coord->J(0, j) = flux_tube_broadening(0, j, 0);
    }

    diagnose = geo_options["diagnose"]
                    .doc("Output additional diagnostics?")
                    .withDefault<bool>(false);

}

void FieldlineGeometry::transform(Options& state) {

}

void FieldlineGeometry::outputVars(Options& state) {
    AUTO_TRACE();
    if (diagnose) {

        set_with_attrs(
            state[std::string("fieldline_geometry_parallel_length")], lpar,
            {{"units", "m"},
            {"long_name", "Parallel length"},
            {"source", "fieldline_geometry"}});

        set_with_attrs(
            state[std::string("fieldline_geometry_geometric_broadening")], geometric_broadening,
            {{"units", ""},
            {"long_name", "R(lpar) / R(lpar=0) for lpar in [m]"},
            {"source", "fieldline_geometry"}});
        
        set_with_attrs(
            state[std::string("fieldline_geometry_transport_broadening")], transport_broadening,
            {{"units", ""},
            {"long_name", "lambda_INT(lpar) / lambda_q(lpar=0) for lpar in [m]"},
            {"source", "fieldline_geometry"}});

        set_with_attrs(
            state[std::string("fieldline_geometry_flux_expansion")], flux_expansion,
            {{"units", ""},
            {"long_name", "(B_pol/B_tot(lpar=0)) / (B_pol/B_tot(lpar)) for lpar in [m]"},
            {"source", "fieldline_geometry"}});
        
        set_with_attrs(
            state[std::string("fieldline_geometry_flux_tube_broadening")], flux_tube_broadening,
            {{"units", ""},
            {"long_name", "Total factor increase for effective flux tube cross-sectional area."},
            {"source", "fieldline_geometry"}});

}}