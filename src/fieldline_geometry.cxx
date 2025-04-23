#include "../include/fieldline_geometry.hxx"
#include <bout/constants.hxx>
#include <bout/coordinates.hxx>
#include <bout/mesh.hxx>
#include <bout/field_factory.hxx>
#include <iostream> // For outputting to log file

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
    Coordinates *coord = mesh->getCoordinates();

    const auto& units = options["units"];
    Lnorm = get<BoutReal>(units["meters"]);

    lpar = calculateLpar();

    lambda_q_omp = geo_options["lambda_q_omp"]
      .doc("Heat-flux decay length at outboard midplane [m]")
      .as<BoutReal>();

    auto B_poloidal_str = geo_options["B_poloidal"]
      .doc("Function for B_poloidal(lpar) (for B_poloidal in [T] and lpar in [m]).")
      .as<std::string>();
    FieldGeneratorPtr B_poloidal_function = FieldFactory::get()->parse(B_poloidal_str, &geo_options);

    auto B_toroidal_str = geo_options["B_toroidal"]
      .doc("Function for B_toroidal(lpar) (for B_toroidal in [T] and lpar in [m]).")
      .as<std::string>();
    FieldGeneratorPtr B_toroidal_function = FieldFactory::get()->parse(B_toroidal_str, &geo_options);

    auto major_radius_str = geo_options["major_radius"]
      .doc("Function for major_radius(lpar) (for major_radius in [m] and lpar in [m]).")
      .as<std::string>();
    FieldGeneratorPtr major_radius_function = FieldFactory::get()->parse(major_radius_str, &geo_options);

    auto effective_flux_exp_str = geo_options["effective_flux_exp"]
      .doc("Function for effective_flux_exp(lpar) (for effective_flux_exp in [~] and lpar in [m]).")
      .as<std::string>();
    FieldGeneratorPtr effective_flux_exp_function = FieldFactory::get()->parse(effective_flux_exp_str, &geo_options);

    B_poloidal.allocate();
    B_toroidal.allocate();
    major_radius.allocate();
    effective_flux_exp.allocate();

    BOUT_FOR(i, lpar.getRegion("RGN_ALL")) {
        B_poloidal[i] = B_poloidal_function->generate(bout::generator::Context().set("lpar", lpar[i]));
        B_toroidal[i] = B_toroidal_function->generate(bout::generator::Context().set("lpar", lpar[i]));
        major_radius[i] = major_radius_function->generate(bout::generator::Context().set("lpar", lpar[i]));
        effective_flux_exp[i] = effective_flux_exp_function->generate(bout::generator::Context().set("lpar", lpar[i]));
    }

    B_total = sqrt(B_toroidal*B_toroidal + B_poloidal*B_poloidal);
    pitch_angle = B_poloidal / B_total;

    auto dy = coord->dy;
    area_external = pitch_angle * dy * 2.0 * PI * major_radius;
    
    BoutReal upstream_pitch_angle = pitch_angle(0, mesh->ystart, 0);
    BoutReal upstream_effective_flux_exp = effective_flux_exp(0, mesh->ystart, 0);

    // "lambda_q" of equation 6 of Derks 2024
    // N.b. B_trans is the inverse of effective_flux_exp.
    poloidal_SOL_width = lambda_q_omp * (upstream_pitch_angle / pitch_angle) * (effective_flux_exp / upstream_effective_flux_exp);

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
            state[std::string("fieldline_geometry_B_poloidal")], B_poloidal,
            {{"units", "T"},
            {"long_name", "B poloidal"},
            {"source", "fieldline_geometry"}});
        
        set_with_attrs(
            state[std::string("fieldline_geometry_B_toroidal")], B_toroidal,
            {{"units", "T"},
            {"long_name", "B toroidal"},
            {"source", "fieldline_geometry"}});

        set_with_attrs(
            state[std::string("fieldline_geometry_major_radius")], major_radius,
            {{"units", "m"},
            {"long_name", "major radius"},
            {"source", "fieldline_geometry"}});
        
        set_with_attrs(
            state[std::string("fieldline_geometry_effective_flux_exp")], effective_flux_exp,
            {{"units", ""},
            {"long_name", "effective flux expansion due to cross-field transport"},
            {"source", "fieldline_geometry"}});
        
        set_with_attrs(
            state[std::string("fieldline_geometry_pitch_angle")], pitch_angle,
            {{"units", ""},
            {"long_name", "B-poloidal / B-total"},
            {"source", "fieldline_geometry"}});
        
        set_with_attrs(
            state[std::string("fieldline_geometry_area_external")], area_external,
            {{"units", "m^2"},
            {"long_name", "area perpendicular to poloidal SOL width"},
            {"source", "fieldline_geometry"}});
        
        set_with_attrs(
            state[std::string("fieldline_geometry_poloidal_SOL_width")], poloidal_SOL_width,
            {{"units", "m"},
            {"long_name", "poloidal scrape-off-layer width"},
            {"source", "fieldline_geometry"}});

}}