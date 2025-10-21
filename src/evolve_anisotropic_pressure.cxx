
#include <bout/constants.hxx>
#include <bout/derivs.hxx>
#include <bout/difops.hxx>
#include <bout/field_factory.hxx>
#include <bout/fv_ops.hxx>
#include <bout/initialprofiles.hxx>
#include <bout/invert_pardiv.hxx>
#include <bout/output_bout_types.hxx>

#include "../include/div_ops.hxx"
#include "../include/evolve_anisotropic_pressure.hxx"
#include "../include/hermes_build_config.hxx"
#include "../include/hermes_utils.hxx"

using bout::globals::mesh;

EvolveAnisotropicPressure::EvolveAnisotropicPressure(std::string name,
                                                     Options& alloptions, Solver* solver)
    : name(name) {
  AUTO_TRACE();

  auto& options = alloptions[name];

  adiabatic_index =
      options["adiabatic_index"]
          .doc("Ratio of specific heats γ = Cp/Cv [5/3 for monatomic ideal gas]")
          .withDefault(5. / 3);
  Cv = 1. / (adiabatic_index - 1.);
  
  density_floor = options["density_floor"].doc("Minimum density floor").withDefault(1e-7);

  temperature_floor = options["temperature_floor"]
                          .doc("Low temperature scale for low_T_diffuse_perp")
                          .withDefault<BoutReal>(0.1)
                      / get<BoutReal>(alloptions["units"]["eV"]);

  pressure_floor = density_floor * temperature_floor;

  // Evolve the energy
  solver->add(E, std::string("E") + name);
  // Evolve the pressure anisotropy Ppar - Pperp
  solver->add(PA, std::string("PA") + name);

  bndry_flux = options["bndry_flux"]
                   .doc("Allow flows through radial boundaries")
                   .withDefault<bool>(true);

  poloidal_flows =
      options["poloidal_flows"].doc("Include poloidal ExB flow").withDefault<bool>(true);

  diagnose = options["diagnose"]
                 .doc("Save additional output diagnostics")
                 .withDefault<bool>(false);

  const Options& units = alloptions["units"];
  const BoutReal Nnorm = units["inv_meters_cubed"];
  const BoutReal Tnorm = units["eV"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();

  auto& p_options = alloptions[std::string("P") + name];
  source_normalisation =
      SI::qe * Nnorm * Tnorm * Omega_ci; // [Pa/s] or [W/m^3] if converted to energy
  time_normalisation = 1. / Omega_ci;    // [s]

  // Try to read the pressure source from the mesh
  // Units of Pascals per second
  source = 0.0;
  mesh->get(source, std::string("P") + name + "_src");
  // Allow the user to override the source
  source = p_options["source"]
               .doc(std::string("Source term in ddt(P") + name
                    + std::string("). Units [Pa/s], note P = 2/3 E"))
               .withDefault(source)
           / (source_normalisation);

  source_time_dependent = p_options["source_time_dependent"]
                              .doc("Use a time-dependent source?")
                              .withDefault<bool>(false);

  // If time dependent, parse the function with respect to time from the input file
  if (source_time_dependent) {
    auto str = p_options["source_prefactor"]
                   .doc("Time-dependent function of multiplier on ddt(P" + name
                        + std::string(") source."))
                   .as<std::string>();
    source_prefactor_function = FieldFactory::get()->parse(str, &p_options);
  }

  if (p_options["source_only_in_core"]
          .doc("Zero the source outside the closed field-line region?")
          .withDefault<bool>(false)) {
    for (int x = mesh->xstart; x <= mesh->xend; x++) {
      if (!mesh->periodicY(x)) {
        // Not periodic, so not in core
        for (int y = mesh->ystart; y <= mesh->yend; y++) {
          for (int z = mesh->zstart; z <= mesh->zend; z++) {
            source(x, y, z) = 0.0;
          }
        }
      }
    }
  }

  neumann_boundary_average_z = p_options["neumann_boundary_average_z"]
                                   .doc("Apply neumann boundary with Z average?")
                                   .withDefault<bool>(false);
}

void EvolveAnisotropicPressure::transform(Options& state) {
  AUTO_TRACE();

  mesh->communicate(E, PA);

  auto& species = state["species"][name];
  N = getNoBoundary<Field3D>(species["density"]);
  const Field3D V = getNoBoundary<Field3D>(species["velocity"]);
  const BoutReal AA = get<BoutReal>(species["AA"]);

  // Calculate pressure
  // E = Cv * P + (1/2) m n v^2
  // P = (2 P_perp + P_par) / 3
  // PA = P_par - P_perp
  // -> P = ( 2 (P_par - PA) + P_par ) / 3 = P_par - (2/3) PA
  //    P_par = P + (2/3) PA
  //
  // -> P = ( 2 P_perp + (PA + P_perp) ) / 3 = P_perp + (1/3) PA
  //    P_perp = P - (1/3) PA
  //
  // => Limit -(3/2) P < PA < 3 P
  P.allocate();
  P_par.allocate();
  P_perp.allocate();
  BOUT_FOR(i, P.getRegion("RGN_ALL")) {
    P[i] = (E[i] - 0.5 * AA * N[i] * SQ(V[i])) / Cv;
    if (P[i] < 0.0) {
      P[i] = 0.0;
    }

    // Limited anisotropy, ensuring that P_perp, P_par >= 0
    BoutReal pa = BOUTMAX(-(3. / 2) * P[i], BOUTMIN(PA[i], 3. * P[i]));

    P_par[i] = P[i] + (2. / 3) * pa;
    P_perp[i] = P[i] - (1. / 3) * pa;
  }
  P.applyBoundary("neumann");
  P_par.applyBoundary("neumann");
  P_perp.applyBoundary("neumann");

  if (neumann_boundary_average_z) {
    // Take Z (usually toroidal) average and apply as X (radial) boundary condition

    auto applyBoundaryAverageZ = [](Field3D& f) {
      if (mesh->firstX()) {
        for (int j = mesh->ystart; j <= mesh->yend; j++) {
          BoutReal favg = 0.0; // Average P in Z
          for (int k = 0; k < mesh->LocalNz; k++) {
            favg += f(mesh->xstart, j, k);
          }
          favg /= mesh->LocalNz;

          // Apply boundary condition
          for (int k = 0; k < mesh->LocalNz; k++) {
            f(mesh->xstart - 1, j, k) = 2. * favg - f(mesh->xstart, j, k);
            f(mesh->xstart - 2, j, k) = f(mesh->xstart - 1, j, k);
          }
        }
      }

      if (mesh->lastX()) {
        for (int j = mesh->ystart; j <= mesh->yend; j++) {
          BoutReal favg = 0.0; // Average P in Z
          for (int k = 0; k < mesh->LocalNz; k++) {
            favg += f(mesh->xend, j, k);
          }
          favg /= mesh->LocalNz;

          for (int k = 0; k < mesh->LocalNz; k++) {
            f(mesh->xend + 1, j, k) = 2. * favg - f(mesh->xend, j, k);
            f(mesh->xend + 2, j, k) = f(mesh->xend + 1, j, k);
          }
        }
      }
    };
    applyBoundaryAverageZ(P);
    applyBoundaryAverageZ(P_par);
    applyBoundaryAverageZ(P_perp);
  }

  // Calculate temperature
  // Not using density boundary condition
  N = getNoBoundary<Field3D>(species["density"]);

  T = P / floor(N, density_floor);
  Field3D T_par = P_par / floor(N, density_floor);
  Field3D T_perp = P_perp / floor(N, density_floor);

  set(species["pressure"], P);
  set(species["pressure_perp"], P_perp);
  set(species["pressure_par"], P_par);
  set(species["temperature"], T);
  set(species["temperature_par"], T_par);
  set(species["temperature_perp"], T_perp);
}

void EvolveAnisotropicPressure::finally(const Options& state) {
  AUTO_TRACE();

  /// Get the section containing this species
  const auto& species = state["species"][name];

  // Get updated pressure and temperature with boundary conditions
  // Note: Retain pressures which fall below zero
  P.clearParallelSlices();
  P.setBoundaryTo(get<Field3D>(species["pressure"]));
  Field3D Pfloor = floor(P, 0.0); // Restricted to never go below zero

  T = get<Field3D>(species["temperature"]);
  N = get<Field3D>(species["density"]);

  if (species.isSet("charge") and (fabs(get<BoutReal>(species["charge"])) > 1e-5)
      and state.isSection("fields") and state["fields"].isSet("phi")) {
    // Electrostatic potential set and species is charged -> include ExB flow

    Field3D phi = get<Field3D>(state["fields"]["phi"]);

    ddt(E) = -Div_n_bxGrad_f_B_XPPM(E, phi, bndry_flux, poloidal_flows, true);

    ddt(PA) = -Div_n_bxGrad_f_B_XPPM(PA, phi, bndry_flux, poloidal_flows, true);
  } else {
    ddt(E) = 0.0;
    ddt(PA) = 0.0;
  }

  if (species.isSet("velocity")) {
    Field3D V = get<Field3D>(species["velocity"]);

    // Typical wave speed used for numerical diffusion
    Field3D fastest_wave;
    if (state.isSet("fastest_wave")) {
      fastest_wave = get<Field3D>(state["fastest_wave"]);
    } else {
      BoutReal AA = get<BoutReal>(species["AA"]);
      fastest_wave = sqrt(T / AA);
    }

    Field3D flow_ylow;
    ddt(E) -= FV::Div_par_mod<hermes::Limiter>(E + P_par, V, fastest_wave, flow_ylow);

    ddt(PA) -= FV::Div_par_mod<hermes::Limiter>(PA, V, fastest_wave, flow_ylow)
               + (2. * P_par + P_perp) * Grad_par(V) - P_perp * Div_par(V);
  }

  //////////////////////
  // Other sources

  if (source_time_dependent) {
    // Evaluate the source_prefactor function at the current time in seconds and scale
    // source with it
    BoutReal time = get<BoutReal>(state["time"]);
    BoutReal source_prefactor =
        source_prefactor_function->generate(bout::generator::Context().set(
            "x", 0, "y", 0, "z", 0, "t", time * time_normalisation));
    final_source = source * source_prefactor;
  } else {
    final_source = source;
  }

  Sp = final_source;
  if (species.isSet("energy_source")) {
    Sp += (2. / 3) * get<Field3D>(species["energy_source"]); // For diagnostic output
  }

  ddt(E) += (3. / 2) * Sp;
}

void EvolveAnisotropicPressure::outputVars(Options& state) {
  AUTO_TRACE();
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto rho_s0 = get<BoutReal>(state["rho_s0"]);

  BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

  state[std::string("E") + name].setAttributes(
      {{"time_dimension", "t"},
       {"units", "J/m^3"},
       {"conversion", Pnorm},
       {"standard_name", "energy density"},
       {"long_name", name + " energy density"},
       {"species", name},
       {"source", "evolve_anisotropic_pressure"}});

  state[std::string("PA") + name].setAttributes(
      {{"time_dimension", "t"},
       {"units", "Pa"},
       {"conversion", Pnorm},
       {"standard_name", "pressure anisotropy"},
       {"long_name", name + " pressure anisotropy"},
       {"species", name},
       {"source", "evolve_anisotropic_pressure"}});

  if (diagnose) {
    set_with_attrs(state[std::string("T") + name], T,
                   {{"time_dimension", "t"},
                    {"units", "eV"},
                    {"conversion", Tnorm},
                    {"standard_name", "temperature"},
                    {"long_name", name + " temperature"},
                    {"species", name},
                    {"source", "evolve_anisotropic_pressure"}});

    set_with_attrs(state[std::string("ddt(E") + name + std::string(")")], ddt(E),
                   {{"time_dimension", "t"},
                    {"units", "J m^-3 s^-1"},
                    {"conversion", Pnorm * Omega_ci},
                    {"long_name", std::string("Rate of change of ") + name + " energy"},
                    {"species", name},
                    {"source", "evolve_anisotropic_pressure"}});

    set_with_attrs(state[std::string("SP") + name], Sp,
                   {{"time_dimension", "t"},
                    {"units", "Pa s^-1"},
                    {"conversion", Pnorm * Omega_ci},
                    {"standard_name", "pressure source"},
                    {"long_name", name + " pressure source"},
                    {"species", name},
                    {"source", "evolve_anisotropic_pressure"}});

    set_with_attrs(state[std::string("P") + name + std::string("_src")], final_source,
                   {{"time_dimension", "t"},
                    {"units", "Pa s^-1"},
                    {"conversion", Pnorm * Omega_ci},
                    {"standard_name", "pressure source"},
                    {"long_name", name + " pressure source"},
                    {"species", name},
                    {"source", "evolve_anisotropic_pressure"}});
  }
}
