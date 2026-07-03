
#include <bout/assert.hxx>
#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/constants.hxx>
#include <bout/derivs.hxx>
#include <bout/difops.hxx>
#include <bout/field.hxx>
#include <bout/field3d.hxx>
#include <bout/fv_ops.hxx>
#include <bout/globals.hxx>
#include <bout/output.hxx>
#include <bout/output_bout_types.hxx>
#include <bout/solver.hxx>

#include "../include/component.hxx"
#include "../include/div_ops.hxx"
#include "../include/guarded_options.hxx"
#include "../include/hermes_build_config.hxx"
#include "../include/hermes_utils.hxx"
#include "../include/neutral_mixed.hxx"
#include "../include/permissions.hxx"

#include <algorithm>
#include <string>

using bout::globals::mesh;

using ParLimiter = hermes::Limiter;

NeutralMixed::NeutralMixed(const std::string& name, Options& alloptions, Solver* solver)
    : Component({readWrite("species:{name}:{outputs}")}), name(name) {

  // Normalisations
  const Options& units = alloptions["units"];
  const BoutReal meters = units["meters"];
  const BoutReal seconds = units["seconds"];
  const BoutReal Nnorm = units["inv_meters_cubed"];
  const BoutReal Tnorm = units["eV"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();
  const BoutReal Cs0 = sqrt(SI::qe * Tnorm / SI::Mp);

  // Need to take derivatives in X for cross-field diffusion terms
  ASSERT0(mesh->xstart > 0);

  auto& options = alloptions[name];

  // Evolving variables e.g name is "h" or "h+"
  solver->add(Nn, std::string("N") + name);
  solver->add(Pn, std::string("P") + name);

  evolve_momentum = options["evolve_momentum"]
                        .doc("Evolve parallel neutral momentum?")
                        .withDefault<bool>(true);

  if (evolve_momentum) {
    solver->add(NVn, std::string("NV") + name);
  } else {
    output_warn.write(
        "WARNING: Not evolving neutral parallel momentum. NVn and Vn set to zero\n");
    NVn = 0.0;
    Vn = 0.0;
  }

  nonorthogonal_operators =
      options["nonorthogonal_operators"]
          .doc("Use nonorthogonal operators for radial transport? NOTE: may be broken")
          .withDefault<bool>(false);

  sheath_ydown = options["sheath_ydown"]
                     .doc("Enable wall boundary conditions at ydown")
                     .withDefault<bool>(true);

  sheath_yup = options["sheath_yup"]
                   .doc("Enable wall boundary conditions at yup")
                   .withDefault<bool>(true);

  density_floor = options["density_floor"]
                      .doc("A minimum density used when dividing NVn by Nn. "
                           "Normalised units.")
                      .withDefault(1e-8);

  freeze_low_density = options["freeze_low_density"]
                           .doc("Freeze evolution in low density regions?")
                           .withDefault<bool>(false);

  temperature_floor = options["temperature_floor"]
                          .doc("Low temperature scale for low_T_diffuse_perp")
                          .withDefault<BoutReal>(0.1)
                      / get<BoutReal>(alloptions["units"]["eV"]);

  pressure_floor = density_floor * temperature_floor;

  precondition = options["precondition"]
                     .doc("Enable preconditioning in neutral model?")
                     .withDefault<bool>(false);

  lax_flux =
      options["lax_flux"].doc("Enable stabilising lax flux?").withDefault<bool>(true);

  neutral_lmax = options["neutral_lmax"]
                     .doc("Maximum neutral mean free path in [m]. Used to calculate "
                          "collisionality floor.")
                     .withDefault(0.1)
                 / meters;

  limiter_gradient_floor = options["limiter_gradient_floor"]
                               .doc("Floor for |grad log Pn| in the D limiter "
                                    "denominator. Normalised inverse-length units.")
                               .withDefault(1.0e-3);

  flux_limit =
      options["flux_limit"]
          .doc("Limit diffusive fluxes to fraction of free-streaming flux. <0 means off.")
          .withDefault(0.5);

  flux_limiter_sharpness = options["flux_limiter_sharpness"]
                               .doc("Sharpness parameter for flux limiter. Higher values "
                                    "make the limiter more abrupt. Must be >0")
                               .withDefault(1.0);

  if (flux_limiter_sharpness <= 0.0) {
    throw BoutException("flux_limiter_sharpness must be > 0.0, got {:g}",
                        flux_limiter_sharpness);
  }

  diffusion_limit = options["diffusion_limit"]
                        .doc("Upper limit on diffusion coefficient [m^2/s]. <0 means off")
                        .withDefault(-1.0)
                    / (meters * meters / seconds); // Normalise

  neutral_viscosity = options["neutral_viscosity"]
                          .doc("Include neutral gas viscosity?")
                          .withDefault<bool>(true);

  neutral_conduction = options["neutral_conduction"]
                           .doc("Include neutral gas heat conduction?")
                           .withDefault<bool>(true);

  collisionality_override =
      options["collisionality_override"]
          .doc(
              "Paramter for overriding the neutral collision frequency in Dn for testing")
          .withDefault(-1.0);

  normalise_sources = options["normalise_sources"]
                          .doc("Normalise input sources?")
                          .withDefault<bool>(true);
  diffusion_collisions_mode = options["diffusion_collisions_mode"]
                                  .doc("Can be multispecies: all enabled collisions "
                                       "excl. IZ, or afn: CX, IZ and NN collisions")
                                  .withDefault<std::string>("multispecies");

  if (precondition) {
    inv = Laplacian::create(&options["precon_laplace"]);

    inv->setInnerBoundaryFlags(INVERT_DC_GRAD | INVERT_AC_GRAD);
    inv->setOuterBoundaryFlags(INVERT_DC_GRAD | INVERT_AC_GRAD);
  }

  // Optionally output time derivatives
  output_ddt =
      options["output_ddt"].doc("Save derivatives to output?").withDefault<bool>(false);

  diagnose =
      options["diagnose"].doc("Save additional diagnostics?").withDefault<bool>(false);

  AA = options["AA"].doc("Particle atomic mass. Proton = 1").withDefault(1.0);

  if (normalise_sources) {
    density_norm = Nnorm * Omega_ci;
    pressure_norm = SI::qe * Nnorm * Tnorm * Omega_ci;
    momentum_norm = SI::Mp * Nnorm * Cs0 * Omega_ci;
  } else {
    density_norm = 1.0;
    pressure_norm = 1.0;
    momentum_norm = 1.0;
  }
  // Try to read the density source from the mesh
  // Units of particles per cubic meter per second
  density_source = 0.0;
  mesh->get(density_source, std::string("N") + name + "_src");
  // Allow the user to override the source
  density_source =
      alloptions[std::string("N") + name]["source"]
          .doc("Source term in ddt(N" + name + std::string("). Units [m^-3/s]"))
          .withDefault(density_source)
      / density_norm;

  // Try to read the pressure source from the mesh
  // Units of Pascals per second
  pressure_source = 0.0;
  mesh->get(pressure_source, std::string("P") + name + "_src");
  // Allow the user to override the source
  pressure_source = alloptions[std::string("P") + name]["source"]
                        .doc(std::string("Source term in ddt(P") + name
                             + std::string("). Units [N/m^2/s]"))
                        .withDefault(pressure_source)
                    / pressure_norm;
  // Try to read the momentum source from the mesh
  momentum_source = 0.0;
  mesh->get(momentum_source, fmt::format("NV{}_src", name));
  // Allow the user to override the source
  momentum_source =
      alloptions[fmt::format("NV{}", name)]["source"]
          .doc(fmt::format("Source term in ddt(NV{}). Units [kg m^-2 s^-2]", name))
          .withDefault(momentum_source)
      / momentum_norm;
  // need some normalisation convention here

  // Set boundary condition defaults: Neumann for all but the diffusivity.
  // The dirichlet on diffusivity ensures no radial flux.
  // NV and V are ignored as they are hardcoded in the parallel BC code.
  alloptions[std::string("Dnn") + name]["bndry_all"] =
      alloptions[std::string("Dnn") + name]["bndry_all"].withDefault("dirichlet");
  alloptions[std::string("T") + name]["bndry_all"] =
      alloptions[std::string("T") + name]["bndry_all"].withDefault("neumann");
  alloptions[std::string("P") + name]["bndry_all"] =
      alloptions[std::string("P") + name]["bndry_all"].withDefault("neumann");
  alloptions[std::string("N") + name]["bndry_all"] =
      alloptions[std::string("N") + name]["bndry_all"].withDefault("neumann");

  // Pick up BCs from input file
  Dnn.setBoundary(std::string("Dnn") + name);
  Tn.setBoundary(std::string("T") + name);
  Pn.setBoundary(std::string("P") + name);
  Nn.setBoundary(std::string("N") + name);

  // All floored versions of variables get the same boundary as the original
  Pnlim.setBoundary(std::string("P") + name);
  logPnlim.setBoundary(std::string("P") + name);
  Nnlim.setBoundary(std::string("N") + name);

  // Product of Dnn and another parameter has same BC as Dnn - see eqns to see why this is
  // necessary
  DnnNn.setBoundary(std::string("Dnn") + name);
  DnnPn.setBoundary(std::string("Dnn") + name);
  DnnNVn.setBoundary(std::string("Dnn") + name);

  substitutePermissions("name", {name});
  substitutePermissions(
      "outputs", {"AA", "density", "pressure", "temperature", "momentum", "velocity"});
}

void NeutralMixed::transform_impl(GuardedOptions& state) {

  mesh->communicate(Nn, Pn, NVn);

  Nn.clearParallelSlices();
  Pn.clearParallelSlices();
  NVn.clearParallelSlices();

  Nn = floor(Nn, 0.0);
  Pn = floor(Pn, 0.0);

  // Nnlim Used where division by neutral density is needed
  Nnlim = softFloor(Nn, density_floor);
  Tn = Pn / Nnlim;
  Tn.applyBoundary();

  Vn = NVn / (AA * Nnlim);
  Vn.applyBoundary("neumann");

  /////////////////////////////////////////////////////
  // Parallel boundary conditions
  TRACE("Neutral boundary conditions");

  if (sheath_ydown) {
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        // Free boundary (constant gradient) density
        const BoutReal nnwall = std::max(
            0.5 * (3. * Nn(r.ind, mesh->ystart, jz) - Nn(r.ind, mesh->ystart + 1, jz)),
            0.0);

        const BoutReal tnwall = Tn(r.ind, mesh->ystart, jz);

        Nn(r.ind, mesh->ystart - 1, jz) = 2 * nnwall - Nn(r.ind, mesh->ystart, jz);

        // Zero gradient temperature, heat flux added later
        Tn(r.ind, mesh->ystart - 1, jz) = tnwall;

        // Set pressure consistent at the boundary
        // Pn(r.ind, mesh->ystart - 1, jz) =
        //     2. * nnwall * tnwall - Pn(r.ind, mesh->ystart, jz);

        // Zero-gradient pressure
        Pn(r.ind, mesh->ystart - 1, jz) = Pn(r.ind, mesh->ystart, jz);

        // No flow into wall
        Vn(r.ind, mesh->ystart - 1, jz) = -Vn(r.ind, mesh->ystart, jz);
        NVn(r.ind, mesh->ystart - 1, jz) = -NVn(r.ind, mesh->ystart, jz);
      }
    }
  }

  if (sheath_yup) {
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        // Free boundary (constant gradient) density
        const BoutReal nnwall = std::max(
            0.5 * (3. * Nn(r.ind, mesh->yend, jz) - Nn(r.ind, mesh->yend - 1, jz)), 0.0);

        const BoutReal tnwall = Tn(r.ind, mesh->yend, jz);

        Nn(r.ind, mesh->yend + 1, jz) = 2 * nnwall - Nn(r.ind, mesh->yend, jz);

        // Zero gradient temperature, heat flux added later
        Tn(r.ind, mesh->yend + 1, jz) = tnwall;

        // Zero-gradient pressure
        Pn(r.ind, mesh->yend + 1, jz) = Pn(r.ind, mesh->yend, jz);

        // No flow into wall
        Vn(r.ind, mesh->yend + 1, jz) = -Vn(r.ind, mesh->yend, jz);
        NVn(r.ind, mesh->yend + 1, jz) = -NVn(r.ind, mesh->yend, jz);
      }
    }
  }

  // Set values in the state
  auto localstate = state["species"][name];
  set(localstate["density"], Nn);
  set(localstate["AA"], AA); // Atomic mass
  set(localstate["pressure"], Pn);
  set(localstate["momentum"], NVn);
  set(localstate["velocity"], Vn);
  set(localstate["temperature"], Tn);
}

void NeutralMixed::finally(const Options& state) {
  const auto& localstate = state["species"][name];

  // extract auxiliary variables derived from
  // Nn, Pn, NVn, from the local state
  // and set boundary conditions on evolved quantities
  Tn = get<Field3D>(localstate["temperature"]);
  Vn = get<Field3D>(localstate["velocity"]);
  Pn = get<Field3D>(localstate["pressure"]);
  Nn = get<Field3D>(localstate["density"]);
  NVn = get<Field3D>(localstate["momentum"]);

  // Logarithms used to calculate perpendicular velocity
  // V_perp = -Dnn * ( Grad_perp(Nn)/Nn + Grad_perp(Tn)/Tn )
  //
  // Grad(Pn) / Pn = Grad(Tn)/Tn + Grad(Nn)/Nn
  //               = Grad(logTn + logNn)
  // Field3D logNn = log(Nn);
  // Field3D logTn = log(Tn);

  // Nnlim Used where division by neutral density is needed
  Nnlim = softFloor(Nn, density_floor);
  // Tnlim used where positivity of Tn is required
  const Field3D Tnlim = softFloor(Tn, temperature_floor);
  // Pnlim used where positivity of Pn is required
  Pnlim = softFloor(Pn, pressure_floor);
  logPnlim = log(Pnlim);
  logPnlim.applyBoundary();
  ///////////////////////////////////////////////////////
  // Calculate cross-field diffusion from collision frequency
  //
  //
  // Pseudo-collisionality. A collision frequency calculated from user-set
  // maximum neutral mean free path to use as a floor for neutral collisionality.
  const Field3D nu_pseudo_mfp = sqrt(Tnlim / AA) / neutral_lmax;

  if (collisionality_override > 0.0) {
    // user has set an override for collision frequency
    Dnn_unlimited = (Tn / AA) / collisionality_override;
  } else {
    if (localstate.isSet("collision_frequency")) {
      // Collisionality
      // Braginskii mode: plasma - self collisions and ei, neutrals - CX, IZ
      if (collision_names.empty()) { // Calculate only once - at the beginning

        if (diffusion_collisions_mode == "afn") {
          for (const auto& collision :
               localstate["collision_frequencies"].getChildren()) {

            const std::string collision_name = collision.second.name();

            if ( // Charge exchange
                (collisionSpeciesMatch(collision_name, name, "+", "cx", "partial")) or
                // Ionisation
                (collisionSpeciesMatch(collision_name, name, "+", "iz", "partial")) or
                // Neutral-neutral collisions
                (collisionSpeciesMatch(collision_name, name, name, "coll", "exact"))) {
              collision_names.push_back(collision_name);
            }
          }
          // Multispecies mode: all collisions and CX are included
        } else if (diffusion_collisions_mode == "multispecies") {
          for (const auto& collision :
               localstate["collision_frequencies"].getChildren()) {

            const std::string collision_name = collision.second.name();

            if ( // Charge exchange
                (collisionSpeciesMatch(collision_name, name, "", "cx", "partial")) or
                // Any collision (en, in, ee, ii, nn)
                (collisionSpeciesMatch(collision_name, name, "", "coll", "partial"))) {
              collision_names.push_back(collision_name);
            }
          }

        } else {
          throw BoutException("\ndiffusion_collisions_mode for {:s} must be either "
                              "multispecies or braginskii",
                              name);
        }

        if (collision_names.empty()) {
          throw BoutException("\tNo collisions found for {:s} in neutral_mixed for "
                              "selected collisions mode",
                              name);
        }

        // Write chosen collisions to log file
        output_info.write("\t{:s} neutral collisionality mode: '{:s}' using ", name,
                          diffusion_collisions_mode);
        for (const auto& collision : collision_names) {
          output_info.write("{:s} ", collision);
        }
        output_info.write("\n");
      }

      // Collect the collisionalities based on list of names
      nu = 0;
      for (const auto& collision_name : collision_names) {
        nu += GET_VALUE(Field3D, localstate["collision_frequencies"][collision_name]);
      }

      // Floor collisionality to avoid unphysically large Dnn at low plasma densities
      Field3D nu_floored = nu + nu_pseudo_mfp * exp(-nu / nu_pseudo_mfp);

      // Dnn = Vth^2 / sigma
      Dnn_unlimited = (Tnlim / AA) / nu_floored;
    } else {
      Dnn_unlimited = (Tnlim / AA) / nu_pseudo_mfp;
    }
  }

  // Start from the unlimited coefficient. copy() forces independent storage, so
  // that later writes to Dnn / Dmax do not alias back onto Dnn_unlimited.
  Dnn = copy(Dnn_unlimited);
  Dmax = copy(Dnn_unlimited);

  // Flux limit: cap diffusion at a fraction of the free-streaming particle flux,
  // set through the ceiling coefficient Dmax.
  if (flux_limit > 0.0) {

    // Thermal speed in a non-drifting Maxwellian [Stangeby eq. 2.21, p.67]
    const Field3D Vn_th = sqrt(8.0 * Tnlim / (PI * AA));

    // Particle flux is 0.25 * Nn * Vth [Stangeby, under eq. 2.24, p.67]; the Nn
    // factor enters at the operator. The denominator uses a gradient floor:
    // smoothly differentiable and avoids division by zero at shallow gradients.
    const Vector3D grad_logPnlim = Grad_perp(logPnlim);
    const Field3D square_gradient = grad_logPnlim * grad_logPnlim;
    Dmax = flux_limit * 0.25 * Vn_th / sqrt(square_gradient + SQ(limiter_gradient_floor));
  }

  // Hard upper limit on the diffusion coefficient. Clamp Dmax down to whichever
  // cap is tighter, so both limiters compose through the single blend below.
  if (diffusion_limit > 0.0) {
    BOUT_FOR(i, Dmax.getRegion("RGN_ALL")) {
      Dmax[i] = BOUTMIN(Dmax[i], diffusion_limit);
    }
  }

  // Blend the unlimited coefficient with the ceiling (harmonic mean). Skipped
  // when no limiter is active, leaving Dnn = Dnn_unlimited from the copy above.
  if (flux_limit > 0.0 || diffusion_limit > 0.0) {

    if (flux_limiter_sharpness == 1.0) {
      // Avoid expensive pow() if not needed
      BOUT_FOR(i, Dnn.getRegion("RGN_NOBNDRY")) {
        Dnn[i] = Dnn_unlimited[i] * Dmax[i] / (Dnn_unlimited[i] + Dmax[i]);
      }
    } else {
      BOUT_FOR(i, Dnn.getRegion("RGN_NOBNDRY")) {
        Dnn[i] = Dnn_unlimited[i]
                 * pow(1.0 + pow(Dnn_unlimited[i] / Dmax[i], flux_limiter_sharpness),
                       -1.0 / flux_limiter_sharpness);
      }
    }
  }

  mesh->communicate(Dnn);
  Dnn.clearParallelSlices();
  Dnn.applyBoundary();

  // Neutral diffusion parameters have the same boundary condition as Dnn
  DnnNn = Dnn * Nnlim;
  DnnPn = Dnn * Pnlim;
  DnnNVn = Dnn * NVn;

  DnnPn.applyBoundary();
  DnnNn.applyBoundary();
  DnnNVn.applyBoundary();

  if (sheath_ydown) {
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        Dnn(r.ind, mesh->ystart - 1, jz) = -Dnn(r.ind, mesh->ystart, jz);
        DnnNn(r.ind, mesh->ystart - 1, jz) = -DnnNn(r.ind, mesh->ystart, jz);
        DnnPn(r.ind, mesh->ystart - 1, jz) = -DnnPn(r.ind, mesh->ystart, jz);
        DnnNVn(r.ind, mesh->ystart - 1, jz) = -DnnNVn(r.ind, mesh->ystart, jz);
      }
    }
  }

  if (sheath_yup) {
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        Dnn(r.ind, mesh->yend + 1, jz) = -Dnn(r.ind, mesh->yend, jz);
        DnnNn(r.ind, mesh->yend + 1, jz) = -DnnNn(r.ind, mesh->yend, jz);
        DnnPn(r.ind, mesh->yend + 1, jz) = -DnnPn(r.ind, mesh->yend, jz);
        DnnNVn(r.ind, mesh->yend + 1, jz) = -DnnNVn(r.ind, mesh->yend, jz);
      }
    }
  }

  // Sound speed appearing in Lax flux for advection terms
  sound_speed = 0;
  if (lax_flux) {
    if (state.isSet("fastest_wave")) {
      sound_speed = get<Field3D>(state["fastest_wave"]);
    } else {
      sound_speed = sqrt(Tn * (5. / 3) / AA);
    }
  }

  // Heat conductivity
  // Note: This is kappa_n = (5/2) * Pn / (m * nu)
  //       where nu is the collision frequency used in Dnn
  kappa_n = (5. / 2) * DnnNn;

  // Viscosity
  // Relationship between heat conduction and viscosity for neutral
  // gas Chapman, Cowling "The Mathematical Theory of Non-Uniform
  // Gases", CUP 1952 Ferziger, Kaper "Mathematical Theory of
  // Transport Processes in Gases", 1972
  // eta_n = (2. / 5) * m_n * kappa_n;
  //
  eta_n = AA * (2. / 5) * kappa_n;

  /////////////////////////////////////////////////////
  // Neutral density
  TRACE("Neutral density");
  ddt(Nn) = -FV::Div_par_mod<ParLimiter>(Nn, Vn, sound_speed,
                                         pf_adv_par_ylow); // Parallel advection

  // Perpendicular diffusion
  if (nonorthogonal_operators) {
    ddt(Nn) +=
        Div_a_Grad_perp_nonorthog(DnnNn, logPnlim, pf_adv_perp_xlow, pf_adv_perp_ylow);
  } else {
    ddt(Nn) += Div_a_Grad_perp_flows(DnnNn, logPnlim, pf_adv_perp_xlow, pf_adv_perp_ylow);
  }

  Sn = density_source; // Save for possible output
  if (localstate.isSet("density_source")) {
    Sn += get<Field3D>(localstate["density_source"]);
  }
  ddt(Nn) += Sn; // Always add density_source

  /////////////////////////////////////////////////////
  // Neutral pressure
  TRACE("Neutral pressure");

  ddt(Pn) = -(5. / 3)
                * FV::Div_par_mod<ParLimiter>( // Parallel advection
                    Pn, Vn, sound_speed, ef_adv_par_ylow)
            + (2. / 3) * Vn * Grad_par(Pn); // Work done

  // Perpendicular advection of pressure
  if (nonorthogonal_operators) {
    ddt(Pn) +=
        (5. / 3)
        * Div_a_Grad_perp_nonorthog(DnnPn, logPnlim, ef_adv_perp_xlow, ef_adv_perp_ylow);
  } else {
    ddt(Pn) +=
        (5. / 3)
        * Div_a_Grad_perp_flows(DnnPn, logPnlim, ef_adv_perp_xlow, ef_adv_perp_ylow);
  }

  // The factor here is 5/2 as we're advecting internal energy and pressure.
  ef_adv_par_ylow *= 5. / 2;
  ef_adv_perp_xlow *= 5. / 2;
  ef_adv_perp_ylow *= 5. / 2;

  if (neutral_conduction) {
    ddt(Pn) += (2. / 3)
               * Div_par_K_Grad_par_mod(kappa_n, Tn, // Parallel conduction
                                        ef_cond_par_ylow,
                                        false); // No conduction through target boundary

    // Perpendicular conduction
    if (nonorthogonal_operators) {
      ddt(Pn) +=
          (2. / 3)
          * Div_a_Grad_perp_nonorthog(kappa_n, Tn, ef_cond_perp_xlow, ef_cond_perp_ylow);
    } else {
      ddt(Pn) +=
          (2. / 3)
          * Div_a_Grad_perp_flows(kappa_n, Tn, ef_cond_perp_xlow, ef_cond_perp_ylow);
    }

    // The factor here is likely 3/2 as this is pure energy flow, but needs checking.
    ef_cond_perp_xlow *= 3. / 2;
    ef_cond_perp_ylow *= 3. / 2;
    ef_cond_par_ylow *= 3. / 2;
  }

  Sp = pressure_source;
  if (localstate.isSet("energy_source")) {
    Sp += (2. / 3) * get<Field3D>(localstate["energy_source"]);
  }
  ddt(Pn) += Sp;

  if (evolve_momentum) {

    /////////////////////////////////////////////////////
    // Neutral momentum
    TRACE("Neutral momentum");

    ddt(NVn) = -AA
                   * FV::Div_par_fvv<ParLimiter>( // Momentum flow
                       Nnlim, Vn, sound_speed)

               - Grad_par(Pn); // Pressure gradient

    // Perpendicular advection of momentum
    if (nonorthogonal_operators) {
      ddt(NVn) +=
          Div_a_Grad_perp_nonorthog(DnnNVn, logPnlim, mf_adv_perp_xlow, mf_adv_perp_ylow);
    } else {
      ddt(NVn) +=
          Div_a_Grad_perp_flows(DnnNVn, logPnlim, mf_adv_perp_xlow, mf_adv_perp_ylow);
    }

    if (neutral_viscosity) {
      // NOTE: The following viscosity terms are not (yet) balanced
      //       by a viscous heating term
      // Relationship between heat conduction and viscosity for neutral
      // gas Chapman, Cowling "The Mathematical Theory of Non-Uniform
      // Gases", CUP 1952 Ferziger, Kaper "Mathematical Theory of
      // Transport Processes in Gases", 1972
      // eta_n = (2. / 5) * kappa_n;

      Field3D viscosity_source = Div_par_K_Grad_par_mod( // Parallel viscosity
          eta_n, Vn, mf_visc_par_ylow,
          false) // No viscosity through target boundary
          ;

      // Perpendicular viscosity
      if (nonorthogonal_operators) {
        viscosity_source +=
            Div_a_Grad_perp_nonorthog(eta_n, Vn, mf_visc_perp_xlow, mf_visc_perp_ylow);
      } else {
        viscosity_source +=
            Div_a_Grad_perp_flows(eta_n, Vn, mf_visc_perp_xlow, mf_visc_perp_ylow);
      }

      ddt(NVn) += viscosity_source;
      ddt(Pn) += -(2. / 3) * Vn * viscosity_source;
    }
    Snv = momentum_source;
    if (localstate.isSet("momentum_source")) {
      Snv += get<Field3D>(localstate["momentum_source"]);
    }
    ddt(NVn) += Snv;

  } else {
    ddt(NVn) = 0;
    Snv = 0;
  }

  // Scale time derivatives
  if (state.isSet("scale_timederivs")) {
    Field3D scale_timederivs = get<Field3D>(state["scale_timederivs"]);
    ddt(Nn) *= scale_timederivs;
    ddt(Pn) *= scale_timederivs;
    ddt(NVn) *= scale_timederivs;
  }

  if (freeze_low_density) {
    // Apply a factor to time derivatives in low density regions.
    // Keep the sources and sinks, so that temperature and flow
    // equilibriates with the plasma through collisions.

    Field3D Nn_s, Pn_s, NVn_s;
    if (localstate.isSet("density_source")) {
      Nn_s = get<Field3D>(localstate["density_source"]);
    } else {
      Nn_s = 0.0;
    }
    if (localstate.isSet("energy_source")) {
      Pn_s = (2. / 3) * get<Field3D>(localstate["energy_source"]);
    } else {
      Pn_s = 0.0;
    }
    if (localstate.isSet("momentum_source")) {
      NVn_s = get<Field3D>(localstate["momentum_source"]);
    } else {
      NVn_s = 0.0;
    }

    for (auto& i : Nn.getRegion("RGN_NOBNDRY")) {
      // Local average density.
      // The purpose is to turn on evolution when nearby cells contain significant
      // density.
      const BoutReal meanNn =
          (1. / 6) * (2 * Nn[i] + Nn[i.xp()] + Nn[i.xm()] + Nn[i.yp()] + Nn[i.ym()]);
      const BoutReal factor = exp(-density_floor / meanNn);
      ddt(Nn)[i] = factor * ddt(Nn)[i] + (1. - factor) * Nn_s[i];
      ddt(Pn)[i] = factor * ddt(Pn)[i] + (1. - factor) * Pn_s[i];
      ddt(NVn)[i] = factor * ddt(NVn)[i] + (1. - factor) * NVn_s[i];
    }
  }

#if CHECKLEVEL >= 1
  for (auto& i : Nn.getRegion("RGN_NOBNDRY")) {
    if (!std::isfinite(ddt(Nn)[i])) {
      throw BoutException("ddt(N{}) non-finite at {}\n", name, i);
    }
    if (!std::isfinite(ddt(Pn)[i])) {
      throw BoutException("ddt(P{}) non-finite at {}\n", name, i);
    }
    if (!std::isfinite(ddt(NVn)[i])) {
      throw BoutException("ddt(NV{}) non-finite at {}\n", name, i);
    }
  }
#endif
}

void NeutralMixed::outputVars(Options& state) {
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto Cs0 = get<BoutReal>(state["Cs0"]);
  auto rho_s0 = get<BoutReal>(state["rho_s0"]);
  const BoutReal Pnorm = SI::qe * Tnorm * Nnorm;

  state[std::string("N") + name].setAttributes({{"time_dimension", "t"},
                                                {"units", "m^-3"},
                                                {"conversion", Nnorm},
                                                {"standard_name", "density"},
                                                {"long_name", name + " number density"},
                                                {"species", name},
                                                {"source", "neutral_mixed"}});

  state[std::string("P") + name].setAttributes({{"time_dimension", "t"},
                                                {"units", "Pa"},
                                                {"conversion", Pnorm},
                                                {"standard_name", "pressure"},
                                                {"long_name", name + " pressure"},
                                                {"species", name},
                                                {"source", "neutral_mixed"}});

  state[std::string("NV") + name].setAttributes(
      {{"time_dimension", "t"},
       {"units", "kg / m^2 / s"},
       {"conversion", SI::Mp * Nnorm * Cs0},
       {"standard_name", "momentum"},
       {"long_name", name + " parallel momentum"},
       {"species", name},
       {"source", "neutral_mixed"}});

  set_with_attrs(state[std::string("V") + name], Vn,
                 {{"time_dimension", "t"},
                  {"units", "m / s"},
                  {"conversion", Cs0},
                  {"standard_name", "velocity"},
                  {"long_name", name + " parallel velocity"},
                  {"source", "neutral_mixed"}});

  if (output_ddt) {
    set_with_attrs(
        state[std::string("ddt(N") + name + std::string(")")], ddt(Nn),
        {{"time_dimension", "t"},
         {"units", "m^-3 s^-1"},
         {"conversion", Nnorm * Omega_ci},
         {"long_name", std::string("Rate of change of ") + name + " number density"},
         {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("ddt(P") + name + std::string(")")], ddt(Pn),
                   {{"time_dimension", "t"},
                    {"units", "Pa s^-1"},
                    {"conversion", Pnorm * Omega_ci},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("ddt(NV") + name + std::string(")")], ddt(NVn),
                   {{"time_dimension", "t"},
                    {"units", "kg m^-2 s^-2"},
                    {"conversion", SI::Mp * Nnorm * Cs0 * Omega_ci},
                    {"source", "neutral_mixed"}});
  }
  if (diagnose) {
    set_with_attrs(state[std::string("T") + name], Tn,
                   {{"time_dimension", "t"},
                    {"units", "eV"},
                    {"conversion", Tnorm},
                    {"standard_name", "temperature"},
                    {"long_name", name + " temperature"},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("Dnn") + name], Dnn,
                   {{"time_dimension", "t"},
                    {"units", "m^2/s"},
                    {"conversion", Cs0 * Cs0 / Omega_ci},
                    {"standard_name", "diffusion coefficient"},
                    {"long_name", name + " diffusion coefficient"},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[fmt::format("Dnn{}_unlimited", name)], Dnn_unlimited,
                   {{"time_dimension", "t"},
                    {"units", "m^2/s"},
                    {"conversion", Cs0 * Cs0 / Omega_ci},
                    {"standard_name", "diffusion coefficient"},
                    {"long_name", name + " unlimited diffusion coefficient"},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[fmt::format("Dnn{}_max", name)], Dmax,
                   {{"time_dimension", "t"},
                    {"units", "m^2/s"},
                    {"conversion", Cs0 * Cs0 / Omega_ci},
                    {"standard_name", "diffusion coefficient"},
                    {"long_name", name + " maximum diffusion coefficient"},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("SN") + name], Sn,
                   {{"time_dimension", "t"},
                    {"units", "m^-3 s^-1"},
                    {"conversion", Nnorm * Omega_ci},
                    {"standard_name", "density source"},
                    {"long_name", name + " number density source"},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("SP") + name], Sp,
                   {{"time_dimension", "t"},
                    {"units", "Pa s^-1"},
                    {"conversion", SI::qe * Tnorm * Nnorm * Omega_ci},
                    {"standard_name", "pressure source"},
                    {"long_name", name + " pressure source"},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("SNV") + name], Snv,
                   {{"time_dimension", "t"},
                    {"units", "kg m^-2 s^-2"},
                    {"conversion", SI::Mp * Nnorm * Cs0 * Omega_ci},
                    {"standard_name", "momentum source"},
                    {"long_name", name + " momentum source"},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("S") + name + std::string("_src")], density_source,
                   {{"time_dimension", "t"},
                    {"units", "m^-3 s^-1"},
                    {"conversion", Nnorm * Omega_ci},
                    {"standard_name", "density source"},
                    {"long_name", name + " number density source"},
                    {"species", name},
                    {"source", "neutral_mixed"}});
    set_with_attrs(state[std::string("P") + name + std::string("_src")], pressure_source,
                   {{"time_dimension", "t"},
                    {"units", "Pa s^-1"},
                    {"conversion", Pnorm * Omega_ci},
                    {"standard_name", "pressure source"},
                    {"long_name", name + " pressure source"},
                    {"species", name},
                    {"source", "neutral_mixed"}});

    ///////////////////////////////////////////////////
    // Parallel flow diagnostics

    // Particle flows due to advection
    if (pf_adv_perp_xlow.isAllocated()) {
      set_with_attrs(
          state[fmt::format("pf{}_adv_perp_xlow", name)], pf_adv_perp_xlow,
          {{"time_dimension", "t"},
           {"units", "s^-1"},
           {"conversion", rho_s0 * SQ(rho_s0) * Nnorm * Omega_ci},
           {"standard_name", "particle flow"},
           {"long_name", name + " radial component of perpendicular advection flow."},
           {"species", name},
           {"source", "neutral_mixed"}});
    }
    if (pf_adv_perp_ylow.isAllocated()) {
      set_with_attrs(
          state[fmt::format("pf{}_adv_perp_ylow", name)], pf_adv_perp_ylow,
          {{"time_dimension", "t"},
           {"units", "s^-1"},
           {"conversion", rho_s0 * SQ(rho_s0) * Nnorm * Omega_ci},
           {"standard_name", "particle flow"},
           {"long_name", name + " poloidal component of perpendicular advection flow."},
           {"species", name},
           {"source", "evolve_density"}});
    }
    if (pf_adv_par_ylow.isAllocated()) {
      set_with_attrs(state[fmt::format("pf{}_adv_par_ylow", name)], pf_adv_par_ylow,
                     {{"time_dimension", "t"},
                      {"units", "s^-1"},
                      {"conversion", rho_s0 * SQ(rho_s0) * Nnorm * Omega_ci},
                      {"standard_name", "particle flow"},
                      {"long_name", name + " parallel advection flow."},
                      {"species", name},
                      {"source", "evolve_density"}});
    }

    // Momentum flows due to advection
    if (mf_adv_perp_xlow.isAllocated()) {
      set_with_attrs(
          state[fmt::format("mf{}_adv_perp_xlow", name)], mf_adv_perp_xlow,
          {{"time_dimension", "t"},
           {"units", "N"},
           {"conversion", rho_s0 * SQ(rho_s0) * SI::Mp * Nnorm * Cs0 * Omega_ci},
           {"standard_name", "momentum flow"},
           {"long_name",
            name + " radial component of perpendicular momentum advection flow."},
           {"species", name},
           {"source", "evolve_momentum"}});
    }
    if (mf_adv_perp_ylow.isAllocated()) {
      set_with_attrs(
          state[fmt::format("mf{}_adv_perp_ylow", name)], mf_adv_perp_ylow,
          {{"time_dimension", "t"},
           {"units", "N"},
           {"conversion", rho_s0 * SQ(rho_s0) * SI::Mp * Nnorm * Cs0 * Omega_ci},
           {"standard_name", "momentum flow"},
           {"long_name",
            name + " poloidal component of perpendicular momentum advection flow."},
           {"species", name},
           {"source", "evolve_momentum"}});
    }
    // This one is awaiting flow implementation into Div_par_fvv

    // if (mf_adv_par_ylow.isAllocated()) {
    //   set_with_attrs(state[fmt::format("mf{}_adv_par_ylow", name)], mf_adv_par_ylow,
    //                {{"time_dimension", "t"},
    //                 {"units", "N"},
    //                 {"conversion", rho_s0 * SQ(rho_s0) * SI::Mp * Nnorm * Cs0 *
    //                 Omega_ci},
    //                 {"standard_name", "momentum flow"},
    //                 {"long_name", name + " parallel momentum advection flow. Note: May
    //                 be incomplete."},
    //                 {"species", name},
    //                 {"source", "evolve_momentum"}});
    // }

    // Momentum flows due to viscosity
    if (mf_visc_perp_ylow.isAllocated()) {
      set_with_attrs(
          state[fmt::format("mf{}_visc_perp_ylow", name)], mf_visc_perp_ylow,
          {{"time_dimension", "t"},
           {"units", "N"},
           {"conversion", rho_s0 * SQ(rho_s0) * SI::Mp * Nnorm * Cs0 * Omega_ci},
           {"standard_name", "momentum flow"},
           {"long_name", name + " poloidal component of perpendicular viscosity."},
           {"species", name},
           {"source", "evolve_momentum"}});
    }
    if (mf_visc_perp_xlow.isAllocated()) {
      set_with_attrs(
          state[fmt::format("mf{}_visc_perp_xlow", name)], mf_visc_perp_xlow,
          {{"time_dimension", "t"},
           {"units", "N"},
           {"conversion", rho_s0 * SQ(rho_s0) * SI::Mp * Nnorm * Cs0 * Omega_ci},
           {"standard_name", "momentum flow"},
           {"long_name", name + " radial component of perpendicular viscosity."},
           {"species", name},
           {"source", "evolve_momentum"}});
    }
    if (mf_visc_par_ylow.isAllocated()) {
      set_with_attrs(
          state[fmt::format("mf{}_visc_par_ylow", name)], mf_visc_par_ylow,
          {{"time_dimension", "t"},
           {"units", "N"},
           {"conversion", rho_s0 * SQ(rho_s0) * SI::Mp * Nnorm * Cs0 * Omega_ci},
           {"standard_name", "momentum flow"},
           {"long_name", name + " parallel viscosity."},
           {"species", name},
           {"source", "evolve_momentum"}});
    }

    // Energy flows due to advection
    if (ef_adv_perp_xlow.isAllocated()) {
      set_with_attrs(
          state[fmt::format("ef{}_adv_perp_xlow", name)], ef_adv_perp_xlow,
          {{"time_dimension", "t"},
           {"units", "W"},
           {"conversion", rho_s0 * SQ(rho_s0) * Pnorm * Omega_ci},
           {"standard_name", "power"},
           {"long_name", name + " radial component of perpendicular energy advection."},
           {"species", name},
           {"source", "evolve_pressure"}});
    }
    if (ef_adv_perp_ylow.isAllocated()) {
      set_with_attrs(
          state[fmt::format("ef{}_adv_perp_ylow", name)], ef_adv_perp_ylow,
          {{"time_dimension", "t"},
           {"units", "W"},
           {"conversion", rho_s0 * SQ(rho_s0) * Pnorm * Omega_ci},
           {"standard_name", "power"},
           {"long_name", name + " poloidal component of perpendicular energy advection."},
           {"species", name},
           {"source", "evolve_pressure"}});
    }
    if (ef_adv_par_ylow.isAllocated()) {
      set_with_attrs(state[fmt::format("ef{}_adv_par_ylow", name)], ef_adv_par_ylow,
                     {{"time_dimension", "t"},
                      {"units", "W"},
                      {"conversion", rho_s0 * SQ(rho_s0) * Pnorm * Omega_ci},
                      {"standard_name", "power"},
                      {"long_name", name + " parallel energy advection."},
                      {"species", name},
                      {"source", "evolve_pressure"}});
    }

    // Energy flows due to conduction
    if (ef_cond_perp_xlow.isAllocated()) {
      set_with_attrs(
          state[fmt::format("ef{}_cond_perp_xlow", name)], ef_cond_perp_xlow,
          {{"time_dimension", "t"},
           {"units", "W"},
           {"conversion", rho_s0 * SQ(rho_s0) * Pnorm * Omega_ci},
           {"standard_name", "power"},
           {"long_name", name + " radial component of perpendicular conduction."},
           {"species", name},
           {"source", "evolve_pressure"}});
    }
    if (ef_cond_perp_ylow.isAllocated()) {
      set_with_attrs(
          state[fmt::format("ef{}_cond_perp_ylow", name)], ef_cond_perp_ylow,
          {{"time_dimension", "t"},
           {"units", "W"},
           {"conversion", rho_s0 * SQ(rho_s0) * Pnorm * Omega_ci},
           {"standard_name", "power"},
           {"long_name", name + " poloidal component of perpendicular conduction."},
           {"species", name},
           {"source", "evolve_pressure"}});
    }
    if (ef_cond_par_ylow.isAllocated()) {
      set_with_attrs(state[fmt::format("ef{}_cond_par_ylow", name)], ef_cond_par_ylow,
                     {{"time_dimension", "t"},
                      {"units", "W"},
                      {"conversion", rho_s0 * SQ(rho_s0) * Pnorm * Omega_ci},
                      {"standard_name", "power"},
                      {"long_name", name + " parallel conduction."},
                      {"species", name},
                      {"source", "evolve_pressure"}});
    }
  }
}

void NeutralMixed::precon([[maybe_unused]] const Options& state, BoutReal gamma) {
  if (!precondition) {
    return;
  }

  // First matrix
  //   ( I   0)
  //   (-LE  I)

  Field3D DTdtN = Dnn * Tn * ddt(Nn);
  mesh->communicate(DTdtN);
  DTdtN.applyBoundary("dirichlet");

  ddt(Pn) -= (gamma * 5. / 3) * FV::Div_a_Grad_perp(DTdtN, logPnlim);

  // Second matrix: Invert Pshur
  //   (E^-1   0  )
  //   ( 0    P^-1)
  //
  // d Laplace_perp(x) + a x + (1/c1)Grad(c2) dot Grad_perp(x) = b
  inv->setCoefA(1 - gamma * FV::Div_a_Grad_perp(Dnn, logPnlim));
  inv->setCoefC1(-1. / ((gamma * 5. / 3) * Dnn));
  inv->setCoefC2(logPnlim);
  inv->setCoefD((-gamma * 5. / 3) * Dnn);

  // inv->setInnerBoundaryFlags(INVERT_DC_GRAD);
  // inv->setOuterBoundaryFlags(INVERT_DC_GRAD);

  ddt(Pn) = inv->solve(ddt(Pn));
  mesh->communicate(ddt(Pn));
  ddt(Pn).applyBoundary("dirichlet");

  // Third matrix: update Nn and NVn equations
  // ( I   E^-1U )
  // ( 0     I   )

  ddt(Nn) -= gamma * FV::Div_a_Grad_perp(DnnNn / Pnlim, ddt(Pn));

  if (evolve_momentum) {
    ddt(NVn) -= gamma * FV::Div_a_Grad_perp(DnnNVn / Pnlim, ddt(Pn));
  }

  for (auto& i : Nn.getRegion("RGN_NOBNDRY")) {
    if (!std::isfinite(ddt(Nn)[i])) {
      throw BoutException("Precon ddt(N{}) non-finite at {}\n", name, i);
    }
    if (!std::isfinite(ddt(Pn)[i])) {
      throw BoutException("Precon ddt(P{}) non-finite at {}\n", name, i);
    }
    if (!std::isfinite(ddt(NVn)[i])) {
      throw BoutException("Precon ddt(NV{}) non-finite at {}\n", name, i);
    }
  }
}
