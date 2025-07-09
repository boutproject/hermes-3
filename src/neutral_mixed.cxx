
#include <bout/constants.hxx>
#include <bout/derivs.hxx>
#include <bout/difops.hxx>
#include <bout/fv_ops.hxx>
#include <bout/output_bout_types.hxx>

#include "../include/hermes_utils.hxx"
#include "../include/div_ops.hxx"
#include "../include/hermes_utils.hxx"
#include "../include/hermes_build_config.hxx"
#include "../include/neutral_mixed.hxx"

using bout::globals::mesh;

using ParLimiter = hermes::Limiter;

NeutralMixed::NeutralMixed(const std::string& name, Options& alloptions, Solver* solver)
    : name(name) {
  AUTO_TRACE();

  // Normalisations
  const Options& units = alloptions["units"];
  const BoutReal meters = units["meters"];
  const BoutReal seconds = units["seconds"];
  const BoutReal Nnorm = units["inv_meters_cubed"];
  const BoutReal Tnorm = units["eV"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();

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

  temperature_floor = options["temperature_floor"].doc("Low temperature scale for low_T_diffuse_perp")
    .withDefault<BoutReal>(0.1) / get<BoutReal>(alloptions["units"]["eV"]);

  pressure_floor = density_floor * temperature_floor;

  precondition = options["precondition"]
                     .doc("Enable preconditioning in neutral model?")
                     .withDefault<bool>(false);

  lax_flux = options["lax_flux"]
                     .doc("Enable stabilising lax flux?")
                     .withDefault<bool>(true);

  flux_limit =
      options["flux_limit"]
          .doc("Limit diffusive fluxes to fraction of thermal speed. <0 means off.")
          .withDefault(0.2);

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

  diffusion_collisions_mode = options["diffusion_collisions_mode"]
      .doc("Can be multispecies: all enabled collisions excl. IZ, or afn: CX, IZ and NN collisions")
      .withDefault<std::string>("multispecies");

  if (precondition) {
    inv = std::unique_ptr<Laplacian>(Laplacian::create(&options["precon_laplace"]));

    inv->setInnerBoundaryFlags(INVERT_DC_GRAD | INVERT_AC_GRAD);
    inv->setOuterBoundaryFlags(INVERT_DC_GRAD | INVERT_AC_GRAD);
  }

  // Optionally output time derivatives
  output_ddt =
      options["output_ddt"].doc("Save derivatives to output?").withDefault<bool>(false);

  diagnose =
      options["diagnose"].doc("Save additional diagnostics?").withDefault<bool>(false);

  AA = options["AA"].doc("Particle atomic mass. Proton = 1").withDefault(1.0);

  // Try to read the density source from the mesh
  // Units of particles per cubic meter per second
  density_source = 0.0;
  mesh->get(density_source, std::string("N") + name + "_src");
  // Allow the user to override the source
  density_source =
      alloptions[std::string("N") + name]["source"]
          .doc("Source term in ddt(N" + name + std::string("). Units [m^-3/s]"))
          .withDefault(density_source)
      / (Nnorm * Omega_ci);

  // Try to read the pressure source from the mesh
  // Units of Pascals per second
  pressure_source = 0.0;
  mesh->get(pressure_source, std::string("P") + name + "_src");
  // Allow the user to override the source
  pressure_source = alloptions[std::string("P") + name]["source"]
                        .doc(std::string("Source term in ddt(P") + name
                             + std::string("). Units [N/m^2/s]"))
                        .withDefault(pressure_source)
                    / (SI::qe * Nnorm * Tnorm * Omega_ci);

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
}

void NeutralMixed::transform(Options& state) {
  AUTO_TRACE();

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

  Pnlim = softFloor(Pn, pressure_floor);
  Pnlim.applyBoundary();

  /////////////////////////////////////////////////////
  // Parallel boundary conditions
  TRACE("Neutral boundary conditions");

  if (sheath_ydown) {
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        // Free boundary (constant gradient) density
        BoutReal nnwall =
            0.5 * (3. * Nn(r.ind, mesh->ystart, jz) - Nn(r.ind, mesh->ystart + 1, jz));
        if (nnwall < 0.0)
          nnwall = 0.0;

        BoutReal tnwall = Tn(r.ind, mesh->ystart, jz);

        Nn(r.ind, mesh->ystart - 1, jz) = 2 * nnwall - Nn(r.ind, mesh->ystart, jz);

        // Zero gradient temperature, heat flux added later
        Tn(r.ind, mesh->ystart - 1, jz) = tnwall;

        // Set pressure consistent at the boundary
        // Pn(r.ind, mesh->ystart - 1, jz) =
        //     2. * nnwall * tnwall - Pn(r.ind, mesh->ystart, jz);

        // Zero-gradient pressure
        Pn(r.ind, mesh->ystart - 1, jz) = Pn(r.ind, mesh->ystart, jz);
        Pnlim(r.ind, mesh->ystart - 1, jz) = Pnlim(r.ind, mesh->ystart, jz);

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
        BoutReal nnwall =
            0.5 * (3. * Nn(r.ind, mesh->yend, jz) - Nn(r.ind, mesh->yend - 1, jz));
        if (nnwall < 0.0)
          nnwall = 0.0;

        BoutReal tnwall = Tn(r.ind, mesh->yend, jz);

        Nn(r.ind, mesh->yend + 1, jz) = 2 * nnwall - Nn(r.ind, mesh->yend, jz);

        // Zero gradient temperature, heat flux added later
        Tn(r.ind, mesh->yend + 1, jz) = tnwall;

        // Zero-gradient pressure
        Pn(r.ind, mesh->yend + 1, jz) = Pn(r.ind, mesh->yend, jz);
        Pnlim(r.ind, mesh->yend + 1, jz) = Pnlim(r.ind, mesh->yend, jz);

        // No flow into wall
        Vn(r.ind, mesh->yend + 1, jz) = -Vn(r.ind, mesh->yend, jz);
        NVn(r.ind, mesh->yend + 1, jz) = -NVn(r.ind, mesh->yend, jz);
      }
    }
  }

  // Set values in the state
  auto& localstate = state["species"][name];
  set(localstate["density"], Nn);
  set(localstate["AA"], AA); // Atomic mass
  set(localstate["pressure"], Pn);
  set(localstate["momentum"], NVn);
  set(localstate["velocity"], Vn);
  set(localstate["temperature"], Tn);
}

void NeutralMixed::finally(const Options& state) {
  AUTO_TRACE();
  auto& localstate = state["species"][name];

  // Logarithms used to calculate perpendicular velocity
  // V_perp = -Dnn * ( Grad_perp(Nn)/Nn + Grad_perp(Tn)/Tn )
  //
  // Grad(Pn) / Pn = Grad(Tn)/Tn + Grad(Nn)/Nn
  //               = Grad(logTn + logNn)
  // Field3D logNn = log(Nn);
  // Field3D logTn = log(Tn);

  logPnlim = log(Pnlim);
  logPnlim.applyBoundary();

  ///////////////////////////////////////////////////////
  // Calculate cross-field diffusion from collision frequency
  //
  //

  Field3D Tnlim = softFloor(Tn, temperature_floor);

  BoutReal neutral_lmax =
    0.1 / get<BoutReal>(state["units"]["meters"]); // Normalised length

  Field3D Rnn =
    sqrt(Tnlim / AA) / neutral_lmax; // Neutral-neutral collisions [normalised frequency]

  if (localstate.isSet("collision_frequency")) {

    // Collisionality
    // Braginskii mode: plasma - self collisions and ei, neutrals - CX, IZ
    if (collision_names.empty()) {     /// Calculate only once - at the beginning

      if (diffusion_collisions_mode == "afn") {
        for (const auto& collision : localstate["collision_frequencies"].getChildren()) {

          std::string collision_name = collision.second.name();

          if (/// Charge exchange
              (collisionSpeciesMatch(    
                collision_name, name, "+", "cx", "partial")) or
              /// Ionisation
              (collisionSpeciesMatch(    
                collision_name, name, "+", "iz", "partial")) or
              /// Neutral-neutral collisions
              (collisionSpeciesMatch(    
                collision_name, name, name, "coll", "exact"))) {
                  collision_names.push_back(collision_name);
                }
        }
      // Multispecies mode: all collisions and CX are included
      } else if (diffusion_collisions_mode == "multispecies") {
        for (const auto& collision : localstate["collision_frequencies"].getChildren()) {

          std::string collision_name = collision.second.name();

          if (/// Charge exchange
              (collisionSpeciesMatch(    
                collision_name, name, "", "cx", "partial")) or
              /// Any collision (en, in, ee, ii, nn)
              (collisionSpeciesMatch(    
                collision_name, name, "", "coll", "partial"))) {
                  collision_names.push_back(collision_name);
                }
        }
        
      } else {
        throw BoutException("\ndiffusion_collisions_mode for {:s} must be either multispecies or braginskii", name);
      }

      if (collision_names.empty()) {
        throw BoutException("\tNo collisions found for {:s} in neutral_mixed for selected collisions mode", name);
      }

      /// Write chosen collisions to log file
      output_info.write("\t{:s} neutral collisionality mode: '{:s}' using ",
                      name, diffusion_collisions_mode);
      for (const auto& collision : collision_names) {        
        output_info.write("{:s} ", collision);
      }
      output_info.write("\n");
      }

    /// Collect the collisionalities based on list of names
    nu = 0;
    for (const auto& collision_name : collision_names) {
      nu += GET_VALUE(Field3D, localstate["collision_frequencies"][collision_name]);
    }


    // Dnn = Vth^2 / sigma
    Dnn = (Tnlim / AA) / (nu + Rnn);
  } else {
    Dnn = (Tnlim / AA) / Rnn;
  }

  if (flux_limit > 0.0) {
    // Apply flux limit to diffusion,
    // using the local thermal speed and pressure gradient magnitude
    Field3D Dmax = flux_limit * sqrt(Tnlim / AA) / (abs(Grad(logPnlim)) + 1. / neutral_lmax);
    BOUT_FOR(i, Dnn.getRegion("RGN_NOBNDRY")) {
      Dnn[i] = Dnn[i] * Dmax[i] / (Dnn[i] + Dmax[i]);
    }
  }

  if (diffusion_limit > 0.0) {
    // Impose an upper limit on the diffusion coefficient
    BOUT_FOR(i, Dnn.getRegion("RGN_NOBNDRY")) {
      Dnn[i] = Dnn[i] * diffusion_limit / (Dnn[i] + diffusion_limit);
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

  ddt(Nn) =
    - FV::Div_par_mod<ParLimiter>(Nn, Vn, sound_speed, pf_adv_par_ylow); // Parallel advection

  ddt(Nn) += Div_a_Grad_perp_flows(DnnNn, logPnlim,
                                   pf_adv_perp_xlow,
                                   pf_adv_perp_ylow);    // Perpendicular advection

  Sn = density_source; // Save for possible output
  if (localstate.isSet("density_source")) {
    Sn += get<Field3D>(localstate["density_source"]);
  }
  ddt(Nn) += Sn; // Always add density_source

  /////////////////////////////////////////////////////
  // Neutral pressure
  TRACE("Neutral pressure");

  ddt(Pn) = - (5. / 3) * FV::Div_par_mod<ParLimiter>(       // Parallel advection
                    Pn, Vn, sound_speed, ef_adv_par_ylow)
            + (2. / 3) * Vn * Grad_par(Pn)                  // Work done
            + (5. / 3) * Div_a_Grad_perp_flows(             // Perpendicular advection
                    DnnPn, logPnlim,
                    ef_adv_perp_xlow, ef_adv_perp_ylow)  
     ;

  // The factor here is 5/2 as we're advecting internal energy and pressure.
  ef_adv_par_ylow  *= 5/2;
  ef_adv_perp_xlow *= 5/2; 
  ef_adv_perp_ylow *= 5/2;

  if (neutral_conduction) {
    ddt(Pn) += (2. / 3) * Div_a_Grad_perp_flows(
                    kappa_n, Tn,                            // Perpendicular conduction
                    ef_cond_perp_xlow, ef_cond_perp_ylow)

            + (2. / 3) * Div_par_K_Grad_par_mod(kappa_n, Tn,           // Parallel conduction 
                      ef_cond_par_ylow,        
                      false)  // No conduction through target boundary
      ;
    // The factor here is likely 3/2 as this is pure energy flow, but needs checking.
    ef_cond_perp_xlow *= 3/2;
    ef_cond_perp_ylow *= 3/2;
    ef_cond_par_ylow *= 3/2;
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

    ddt(NVn) =
        -AA * FV::Div_par_fvv<ParLimiter>(             // Momentum flow
              Nnlim, Vn, sound_speed)                  

        - Grad_par(Pn)                                 // Pressure gradient
        
        + Div_a_Grad_perp_flows(DnnNVn, logPnlim,
                                     mf_adv_perp_xlow,
                                     mf_adv_perp_ylow) // Perpendicular advection
      ;

    if (neutral_viscosity) {
      // NOTE: The following viscosity terms are not (yet) balanced
      //       by a viscous heating term

      // Relationship between heat conduction and viscosity for neutral
      // gas Chapman, Cowling "The Mathematical Theory of Non-Uniform
      // Gases", CUP 1952 Ferziger, Kaper "Mathematical Theory of
      // Transport Processes in Gases", 1972
      // eta_n = (2. / 5) * kappa_n;

      Field3D viscosity_source = AA * Div_a_Grad_perp_flows(
                                eta_n, Vn,              // Perpendicular viscosity
                                mf_visc_perp_xlow,
                                mf_visc_perp_ylow)    
                              
                              + AA * Div_par_K_Grad_par_mod(               // Parallel viscosity 
                                eta_n, Vn,
                                mf_visc_par_ylow,
                                false) // No viscosity through target boundary
                          ;

      ddt(NVn) += viscosity_source;
      ddt(Pn)  += -(2. /3) * Vn * viscosity_source;
    }

    if (localstate.isSet("momentum_source")) {
      Snv = get<Field3D>(localstate["momentum_source"]);
      ddt(NVn) += Snv;
    } else {
      Snv = 0;
    }

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
      // The purpose is to turn on evolution when nearby cells contain significant density.
      const BoutReal meanNn = (1./6) * (2 * Nn[i] + Nn[i.xp()] + Nn[i.xm()] + Nn[i.yp()] + Nn[i.ym()]);
      const BoutReal factor = exp(- density_floor / meanNn);
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
      set_with_attrs(state[fmt::format("pf{}_adv_perp_xlow", name)], pf_adv_perp_xlow,
                   {{"time_dimension", "t"},
                    {"units", "s^-1"},
                    {"conversion", rho_s0 * SQ(rho_s0) * Nnorm * Omega_ci},
                    {"standard_name", "particle flow"},
                    {"long_name", name + " radial component of perpendicular advection flow."},
                    {"species", name},
                    {"source", "neutral_mixed"}});
    }
    if (pf_adv_perp_ylow.isAllocated()) {
      set_with_attrs(state[fmt::format("pf{}_adv_perp_ylow", name)], pf_adv_perp_ylow,
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
      set_with_attrs(state[fmt::format("mf{}_adv_perp_xlow", name)], mf_adv_perp_xlow,
                   {{"time_dimension", "t"},
                    {"units", "N"},
                    {"conversion", rho_s0 * SQ(rho_s0) * SI::Mp * Nnorm * Cs0 * Omega_ci},
                    {"standard_name", "momentum flow"},
                    {"long_name", name + " radial component of perpendicular momentum advection flow."},
                    {"species", name},
                    {"source", "evolve_momentum"}});
    }
    if (mf_adv_perp_ylow.isAllocated()) {
      set_with_attrs(state[fmt::format("mf{}_adv_perp_ylow", name)], mf_adv_perp_ylow,
                   {{"time_dimension", "t"},
                    {"units", "N"},
                    {"conversion", rho_s0 * SQ(rho_s0) * SI::Mp * Nnorm * Cs0 * Omega_ci},
                    {"standard_name", "momentum flow"},
                    {"long_name", name + " poloidal component of perpendicular momentum advection flow."},
                    {"species", name},
                    {"source", "evolve_momentum"}});
    }
    // This one is awaiting flow implementation into Div_par_fvv

    // if (mf_adv_par_ylow.isAllocated()) {
    //   set_with_attrs(state[fmt::format("mf{}_adv_par_ylow", name)], mf_adv_par_ylow,
    //                {{"time_dimension", "t"},
    //                 {"units", "N"},
    //                 {"conversion", rho_s0 * SQ(rho_s0) * SI::Mp * Nnorm * Cs0 * Omega_ci},
    //                 {"standard_name", "momentum flow"},
    //                 {"long_name", name + " parallel momentum advection flow. Note: May be incomplete."},
    //                 {"species", name},
    //                 {"source", "evolve_momentum"}});
    // }


    // Momentum flows due to viscosity
    if (mf_visc_perp_ylow.isAllocated()) {
      set_with_attrs(state[fmt::format("mf{}_visc_perp_ylow", name)], mf_visc_perp_ylow,
                   {{"time_dimension", "t"},
                    {"units", "N"},
                    {"conversion", rho_s0 * SQ(rho_s0) * SI::Mp * Nnorm * Cs0 * Omega_ci},
                    {"standard_name", "momentum flow"},
                    {"long_name", name + " poloidal component of perpendicular viscosity."},
                    {"species", name},
                    {"source", "evolve_momentum"}});
    }
    if (mf_visc_perp_ylow.isAllocated()) {
      set_with_attrs(state[fmt::format("mf{}_visc_perp_ylow", name)], mf_visc_perp_ylow,
                   {{"time_dimension", "t"},
                    {"units", "N"},
                    {"conversion", rho_s0 * SQ(rho_s0) * SI::Mp * Nnorm * Cs0 * Omega_ci},
                    {"standard_name", "momentum flow"},
                    {"long_name", name + " poloidal component of perpendicular viscosity."},
                    {"species", name},
                    {"source", "evolve_momentum"}});
    }
    if (mf_visc_par_ylow.isAllocated()) {
      set_with_attrs(state[fmt::format("mf{}_visc_par_ylow", name)], mf_visc_par_ylow,
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
      set_with_attrs(state[fmt::format("ef{}_adv_perp_xlow", name)], ef_adv_perp_xlow,
                   {{"time_dimension", "t"},
                    {"units", "W"},
                    {"conversion", rho_s0 * SQ(rho_s0) * Pnorm * Omega_ci},
                    {"standard_name", "power"},
                    {"long_name", name + " radial component of perpendicular energy advection."},
                    {"species", name},
                    {"source", "evolve_pressure"}});
    }
    if (ef_adv_perp_ylow.isAllocated()) {
      set_with_attrs(state[fmt::format("ef{}_adv_perp_ylow", name)], ef_adv_perp_ylow,
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
      set_with_attrs(state[fmt::format("ef{}_cond_perp_xlow", name)], ef_cond_perp_xlow,
                   {{"time_dimension", "t"},
                    {"units", "W"},
                    {"conversion", rho_s0 * SQ(rho_s0) * Pnorm * Omega_ci},
                    {"standard_name", "power"},
                    {"long_name", name + " radial component of perpendicular conduction."},
                    {"species", name},
                    {"source", "evolve_pressure"}});
    }
    if (ef_cond_perp_ylow.isAllocated()) {
      set_with_attrs(state[fmt::format("ef{}_cond_perp_ylow", name)], ef_cond_perp_ylow,
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

void NeutralMixed::precon(const Options& state, BoutReal gamma) {
  if (!precondition) {
    return;
  }

  // First matrix
  //   ( I   0)
  //   (-LE  I)

  Field3D DTdtN = Dnn * Tn * ddt(Nn);
  mesh->communicate(DTdtN);
  DTdtN.applyBoundary("dirichlet");

  ddt(Pn) -= (gamma * 5./3) * FV::Div_a_Grad_perp(DTdtN,
                                                  logPnlim);

  // Second matrix: Invert Pshur
  //   (E^-1   0  )
  //   ( 0    P^-1)
  //
  // d Laplace_perp(x) + a x + (1/c1)Grad(c2) dot Grad_perp(x) = b
  inv->setCoefA(1 - gamma * FV::Div_a_Grad_perp(Dnn, logPnlim));
  inv->setCoefC1(-1. / ((gamma * 5./3) * Dnn));
  inv->setCoefC2(logPnlim);
  inv->setCoefD((-gamma * 5./3) * Dnn);

  //inv->setInnerBoundaryFlags(INVERT_DC_GRAD);
  //inv->setOuterBoundaryFlags(INVERT_DC_GRAD);

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
