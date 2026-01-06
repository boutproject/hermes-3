
#include <bout/constants.hxx>
#include <bout/mesh.hxx>
#include <bout/output_bout_types.hxx>
#include <bout/solver.hxx>

#include "../include/div_ops.hxx"
#include "../include/hermes_utils.hxx"
#include "../include/neutral_full_velocity.hxx"

#include <algorithm>

using bout::globals::mesh;

NeutralFullVelocity::NeutralFullVelocity(const std::string& name, Options& alloptions,
                                         Solver* solver)
    : name(name) {
  AUTO_TRACE();

  // This is used in both transform and finally functions
  coord = mesh->getCoordinates();

  // Normalisations
  const Options& units = alloptions["units"];
  const BoutReal meters = units["meters"];
  const BoutReal seconds = units["seconds"];
  const BoutReal Bnorm = units["Tesla"];
  const BoutReal Nnorm = units["inv_meters_cubed"];
  const BoutReal Tnorm = units["eV"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();

  auto& options = alloptions[name];

  AA = options["AA"].doc("Atomic mass number. Proton = 1").as<int>();

  adiabatic_index =
      options["adiabatic_index"]
          .doc("Ratio of specific heats γ = Cp/Cv [5/3 for monatomic ideal gas]")
          .withDefault(5. / 3);

  constant_transport_coef  = 
      options["constant_transport_coef"]
      .doc("Use constant transport coefficients?")
      .withDefault<bool>(false);

  neutral_viscosity =
      options["viscosity"].doc("Kinematic viscosity [m^2/s]").withDefault(1.0)
      / (meters * meters / seconds);

  neutral_conduction =
      options["conduction"].doc("Heat conduction [m^2/s]").withDefault(1.0)
      / (meters * meters / seconds);

  neutral_gamma = options["neutral_gamma"]
                      .doc("Surface heat transmission coefficient")
                      .withDefault(5. / 4);

  density_floor = options["density_floor"]
                      .doc("A minimum density used when dividing by density."
                           "Normalised units.")
                      .withDefault(1e-8);

  temperature_floor = options["temperature_floor"].doc("Low temperature scale for low_T_diffuse_perp")
    .withDefault<BoutReal>(0.1) / get<BoutReal>(alloptions["units"]["eV"]);

  pressure_floor = density_floor * temperature_floor;

  diagnose =
      options["diagnose"].doc("Output additional diagnostics?").withDefault<bool>(false);

  toroidal_flow =
      options["toroidal_flow"].doc("Evolve toroidal flow?").withDefault<bool>(true);

  momentum_advection = options["momentum_advection"]
                           .doc("Include advection of momentum?")
                           .withDefault<bool>(false);

  curved_torus = options["curved_torus"]
                     .doc("Include toroidal curvature in momentum advection?")
                     .withDefault<bool>(true);

  flux_limit =
      options["flux_limit"]
          .doc("Limit diffusive fluxes to fraction of thermal speed. <0 means off.")
          .withDefault(0.2);

  neutral_lmax = options["neutral_lmax"].doc("Maximum length scale due to the present of walls.")
    .withDefault<BoutReal>(0.1) / get<BoutReal>(alloptions["units"]["meters"]); // Normalised length

  diffusion_limit = options["diffusion_limit"]
                        .doc("Upper limit on diffusion coefficient [m^2/s]. <0 means off")
                        .withDefault(-1.0)
                    / (meters * meters / seconds); // Normalise    

  diffusion_collisions_mode = options["diffusion_collisions_mode"]
      .doc("Can be multispecies: all enabled collisions excl. IZ, or afn: CX, IZ and NN collisions")
      .withDefault<std::string>("multispecies");

  zero_timederivs = options["zero_timederivs"]
                          .doc("Set the time derivatives to zero?")
                          .withDefault<bool>(false);

  // Optionally output time derivatives
  output_ddt =
      options["output_ddt"].doc("Save derivatives to output?").withDefault<bool>(false);      

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

  // Note: We evolve v^x, v^y and v^z because these have magnitudes close to 1
  //       whereas v_x, v_y and v_z are >> 1.
  Vn2D.covariant = false; ///< Evolve contravariant components

  // Evolve 2D density, pressure, and velocity
  solver->add(Nn2D, "N" + name);
  solver->add(Pn2D, "P" + name);
  solver->add(Vn2D, "V" + name);

  // Load necessary metrics for non-orth calculation
  Field2D etaxy, cosbeta;
  if (mesh->get(etaxy, "etaxy")) {
    etaxy = 0.0;
  }
  cosbeta = sqrt(1. - SQ(etaxy));

  // Calculate transformation to Cartesian coordinates
  Field2D Zxy, hthe, Bpxy;

  if (mesh->get(Rxy, "Rxy")) {
    throw BoutException("Fluid neutrals model requires Rxy");
  }
  if (mesh->get(Zxy, "Zxy")) {
    throw BoutException("Fluid neutrals model requires Zxy");
  }
  if (mesh->get(hthe, "hthe")) {
    throw BoutException("Fluid neutrals model requires hthe");
  }
  if (mesh->get(Bpxy, "Bpxy")) {
    throw BoutException("Fluid neutrals model requires Bpxy");
  }

  // Normalise
  Rxy /= meters;
  Zxy /= meters;
  hthe /= meters;
  Bpxy /= Bnorm;

  // Sign of poloidal field
  sigma_Bp = (Bpxy(mesh->xstart, mesh->ystart) > 0.0) ? 1.0 : -1.0;
  output.write("\tsigma_Bp = {}\n", sigma_Bp);

  // Axisymmetric neutrals simplifies things considerably...

  Urx.allocate();
  Ury.allocate();
  Uzx.allocate();
  Uzy.allocate();

  Txr.allocate();
  Txz.allocate();
  Tyr.allocate();
  Tyz.allocate();

  for (int i = 0; i < mesh->LocalNx; i++)
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      // Central differencing of coordinates
      BoutReal dRdtheta, dZdtheta;
      if (j == mesh->ystart) {
        dRdtheta = (Rxy(i, j + 1) - Rxy(i, j)) / (coord->dy(i, j));
        dZdtheta = (Zxy(i, j + 1) - Zxy(i, j)) / (coord->dy(i, j));
      } else if (j == mesh->yend) {
        dRdtheta = (Rxy(i, j) - Rxy(i, j - 1)) / (coord->dy(i, j));
        dZdtheta = (Zxy(i, j) - Zxy(i, j - 1)) / (coord->dy(i, j));
      } else {
        dRdtheta = (Rxy(i, j + 1) - Rxy(i, j - 1)) / (2. * coord->dy(i, j));
        dZdtheta = (Zxy(i, j + 1) - Zxy(i, j - 1)) / (2. * coord->dy(i, j));
      }

      // Match to hthe, 1/|Grad y|
      BoutReal h = sqrt(SQ(dRdtheta) + SQ(dZdtheta));
      BoutReal grady = 1.0 / hthe(i, j);
      dRdtheta = dRdtheta / grady / h;
      dZdtheta = dZdtheta / grady / h;

      BoutReal dRdpsi, dZdpsi;
      if (i == 0) {
        // One-sided differences
        dRdpsi = (Rxy(i + 1, j) - Rxy(i, j)) / (coord->dx(i, j));
        dZdpsi = (Zxy(i + 1, j) - Zxy(i, j)) / (coord->dx(i, j));
      } else if (i == (mesh->LocalNx - 1)) {
        // One-sided differences
        dRdpsi = (Rxy(i, j) - Rxy(i - 1, j)) / (coord->dx(i, j));
        dZdpsi = (Zxy(i, j) - Zxy(i - 1, j)) / (coord->dx(i, j));
      } else {
        dRdpsi = (Rxy(i + 1, j) - Rxy(i - 1, j)) / (2. * coord->dx(i, j));
        dZdpsi = (Zxy(i + 1, j) - Zxy(i - 1, j)) / (2. * coord->dx(i, j));
      }

      // Match to Bp, |Grad psi|. NOTE: this only works if
      // X and Y are orthogonal.
      BoutReal dldpsi = sqrt(SQ(dRdpsi) + SQ(dZdpsi)) * cosbeta(i, j); // ~ 1/(R*Bp)
      dRdpsi /= dldpsi * Bpxy(i, j) * Rxy(i, j);
      dZdpsi /= dldpsi * Bpxy(i, j) * Rxy(i, j);

      Urx(i, j) = dRdpsi;
      Ury(i, j) = dRdtheta;
      Uzx(i, j) = dZdpsi;
      Uzy(i, j) = dZdtheta;

      // Poloidal (R,Z) transformation Jacobian
      BoutReal J = dRdpsi * dZdtheta - dZdpsi * dRdtheta;

      Txr(i, j) = dZdtheta / J;
      Txz(i, j) = -dRdtheta / J;
      Tyr(i, j) = -dZdpsi / J;
      Tyz(i, j) = dRdpsi / J;
    }

  Urx.applyBoundary("neumann");
  Ury.applyBoundary("neumann");
  Uzx.applyBoundary("neumann");
  Uzy.applyBoundary("neumann");

  Txr.applyBoundary("neumann");
  Txz.applyBoundary("neumann");
  Tyr.applyBoundary("neumann");
  Tyz.applyBoundary("neumann");

  // Ensure that guard cells are filled and consistent between processors
  mesh->communicate(Urx, Ury, Uzx, Uzy);
  mesh->communicate(Txr, Txz, Tyr, Tyz);
}

/// Modify the given simulation state
void NeutralFullVelocity::transform(Options& state) {
  AUTO_TRACE();
  mesh->communicate(Nn2D, Vn2D, Pn2D);

  // Boundary conditions
  Nn2D.applyBoundary("neumann");
  Pn2D.applyBoundary("neumann");
  Vn2D.applyBoundary("dirichlet");

  // Floored fields, used for rate coefficients
  Nn2D = floor(Nn2D, 0.0);
  Pn2D = floor(Pn2D, 0.0);

  // Non-zero floor can use a differentiable soft floor
  Field2D Nnlim = softFloor(Nn2D, density_floor);

  Tn2D = Pn2D / Nnlim;
  Tn2D.applyBoundary("neumann");

  //////////////////////////////////////////////////////
  // 2D (X-Y) full velocity model
  //
  // Evolves density Nn2D, velocity vector Vn2D and pressure Pn2D
  //

  for (RangeIterator idwn = mesh->iterateBndryLowerY(); !idwn.isDone(); idwn.next()) {
    // Diriclet conditions on Y
    // Since Vn2D is contravariant, this is the poloidal flow
    Vn2D.y(idwn.ind, mesh->ystart - 1) = -Vn2D.y(idwn.ind, mesh->ystart);

    // Neumann boundary condition on X and Z components
    Vn2D.x(idwn.ind, mesh->ystart - 1) = Vn2D.x(idwn.ind, mesh->ystart);
    Vn2D.z(idwn.ind, mesh->ystart - 1) = Vn2D.z(idwn.ind, mesh->ystart);

    // Neumann conditions on density and pressure
    Nn2D(idwn.ind, mesh->ystart - 1) = Nn2D(idwn.ind, mesh->ystart);
    Pn2D(idwn.ind, mesh->ystart - 1) = Pn2D(idwn.ind, mesh->ystart);
  }

  for (RangeIterator idwn = mesh->iterateBndryUpperY(); !idwn.isDone(); idwn.next()) {
    // Diriclet conditions on Y
    Vn2D.y(idwn.ind, mesh->yend + 1) = -Vn2D.y(idwn.ind, mesh->yend);

    // Neumann boundary condition on X and Z components
    Vn2D.x(idwn.ind, mesh->yend + 1) = Vn2D.x(idwn.ind, mesh->yend);
    Vn2D.z(idwn.ind, mesh->yend + 1) = Vn2D.z(idwn.ind, mesh->yend);

    // Neumann conditions on density and pressure
    Nn2D(idwn.ind, mesh->yend + 1) = Nn2D(idwn.ind, mesh->yend);
    Pn2D(idwn.ind, mesh->yend + 1) = Pn2D(idwn.ind, mesh->yend);
  }

  Vn2D_contravariant = Vn2D;
  Vn2D.toCovariant(); // Convenient to work with covariant components

  // Exchange of parallel momentum. This could be done
  // in a couple of ways, but here we use the fact that
  // Vn2D is covariant and b = e_y / (JB) to write:
  //
  // V_{||n} = b dot V_n = Vn2D.y / (JB)
  Vnpar = Vn2D.y / (coord->J * coord->Bxy);

  // Set values in the state
  auto& localstate = state["species"][name];
  set(localstate["density"], Nn2D);
  set(localstate["AA"], AA); // Atomic mass
  set(localstate["pressure"], Pn2D);
  set(localstate["momentum"], Vnpar * Nn2D * AA);
  set(localstate["velocity"], Vnpar); // Parallel velocity
  set(localstate["temperature"], Tn2D);
}

/// Use the final simulation state to update internal state
/// (e.g. time derivatives)
void NeutralFullVelocity::finally(const Options& state) {
  AUTO_TRACE();
  auto& localstate = state["species"][name];

  ///////////////////////////////////////////////////////
  // Calculate transport coefficeints from collision frequency

  Field2D Nn2D_floor = softFloor(Nn2D, density_floor);
  Field2D Pn2D_floor = softFloor(Pn2D, pressure_floor);
  Field2D Tn2D_floor = softFloor(Tn2D, temperature_floor);

  if (!constant_transport_coef) {

    Field2D Rnn = sqrt(Tn2D_floor / AA) / neutral_lmax; // Neutral-neutral collisions [normalised frequency]

    if (localstate.isSet("collision_frequency")) {

      // Collisionality
      // Braginskii mode: plasma - self collisions and ei, neutrals - CX, IZ
      if (collision_names.empty()) {     // Calculate only once - at the beginning

        if (diffusion_collisions_mode == "afn") {
          for (const auto& collision : localstate["collision_frequencies"].getChildren()) {

            std::string collision_name = collision.second.name();

            if (// Charge exchange
                (collisionSpeciesMatch(    
                  collision_name, name, "+", "cx", "partial")) or
                // Ionisation
                (collisionSpeciesMatch(    
                  collision_name, name, "+", "iz", "partial")) or
                // Neutral-neutral collisions
                (collisionSpeciesMatch(    
                  collision_name, name, name, "coll", "exact"))) {
                    collision_names.push_back(collision_name);
                  }
          }
        // Multispecies mode: all collisions and CX are included
        } else if (diffusion_collisions_mode == "multispecies") {
          for (const auto& collision : localstate["collision_frequencies"].getChildren()) {

            std::string collision_name = collision.second.name();

            if (// Charge exchange
                (collisionSpeciesMatch(    
                  collision_name, name, "", "cx", "partial")) or
                // Any collision (en, in, ee, ii, nn)
                (collisionSpeciesMatch(    
                  collision_name, name, "", "coll", "partial"))) {
                    collision_names.push_back(collision_name);
                  }
          }
          
        } else {
          throw BoutException("\ndiffusion_collisions_mode for {:s} must be either multispecies or braginskii", name);
        }

        if (collision_names.empty()) {
          throw BoutException("\tNo collisions found for {:s} in neutral_full_velocity for selected collisions mode", name);
        }

        // Write chosen collisions to log file
        output_info.write("\t{:s} neutral collisionality mode: '{:s}' using ",
                        name, diffusion_collisions_mode);
        for (const auto& collision : collision_names) {        
          output_info.write("{:s} ", collision);
        }
        output_info.write("\n");
        }

      // Collect the collisionalities based on list of names
      nu = 0;
      for (const auto& collision_name : collision_names) {
        nu += DC(get<Field3D>(localstate["collision_frequencies"][collision_name]));
      }


      // Dnn = Vth^2 / sigma
      Dnn = (Tn2D_floor / AA) / (nu + Rnn);
    } else {
      Dnn = (Tn2D_floor / AA) / Rnn;
    }

    if (diffusion_limit > 0.0) {
      // Impose an upper limit on the diffusion coefficient
      BOUT_FOR(i, Dnn.getRegion("RGN_NOBNDRY")) {
        Dnn[i] = Dnn[i] * diffusion_limit / (Dnn[i] + diffusion_limit);
      }
    }

    // Heat conductivity 
    // Note: This is kappa_n = (5/2) * Pn / (m * nu)
    //       where nu is the collision frequency used in Dnn
    kappa_n = (5. / 2) * Dnn * Nn2D_floor;

    // Viscosity
    // Relationship between heat conduction and viscosity for neutral
    // gas Chapman, Cowling "The Mathematical Theory of Non-Uniform
    // Gases", CUP 1952 Ferziger, Kaper "Mathematical Theory of
    // Transport Processes in Gases", 1972
    // eta_n = (2. / 5) * m_n * kappa_n;
    //
    // eta_n = AA * (2. / 5) * kappa_n;


    if (flux_limit > 0.0) {
      // Thermal velocity of neutrals
      Field2D Vnth = sqrt(Tn2D_floor / AA); 
      
      // Apply flux limit to diffusion,
      // using the local thermal speed and pressure gradient magnitude
      Field2D Dmax = flux_limit * Vnth / (abs(Grad(Pn2D_floor)/Pn2D_floor) + 1. / neutral_lmax);
      BOUT_FOR(i, Dnn.getRegion("RGN_NOBNDRY")) {
        Dnn[i] = Dnn[i] * Dmax[i] / (Dnn[i] + Dmax[i]);
      }

      Field2D kappa_n_max = flux_limit * (3.0 / 2.0 * Vnth * Nn2D_floor) / (abs(Grad(Tn2D_floor))/Tn2D_floor + 1. / neutral_lmax);
      BOUT_FOR(i, kappa_n.getRegion("RGN_NOBNDRY")) {
        kappa_n[i] = kappa_n[i] * kappa_n_max[i] / (kappa_n[i] + kappa_n_max[i]);
      }

      // NOTE: How can we make a good estimate of abs(Grad(Vn)) since we have two velocities?
      // Field2D viscosity_factor = 1.0 / (1.0 + eta_n * abs(Grad(Vn)) / (flux_limit * Pn2D_floor)); 
      // BOUT_FOR(i, eta_n.getRegion("RGN_NOBNDRY")) {
      //   eta_n[i] = eta_n[i] * viscosity_factor[i];
      // }

    }

    eta_n = AA * (2. / 5) * kappa_n;


    mesh->communicate(Dnn);
    Dnn.applyBoundary();

    mesh->communicate(kappa_n);
    kappa_n.applyBoundary("neumann");

    mesh->communicate(eta_n);
    eta_n.applyBoundary("neumann");


    for (RangeIterator idwn = mesh->iterateBndryLowerY(); !idwn.isDone(); idwn.next()) {
      // Neumann conditions for transport coefficients
      kappa_n(idwn.ind, mesh->ystart - 1) = kappa_n(idwn.ind, mesh->ystart);
      eta_n(idwn.ind, mesh->ystart - 1) = eta_n(idwn.ind, mesh->ystart);
    }

    for (RangeIterator idwn = mesh->iterateBndryUpperY(); !idwn.isDone(); idwn.next()) {
      // Neumann conditions for transport coefficients
      kappa_n(idwn.ind, mesh->yend + 1) = kappa_n(idwn.ind, mesh->yend);
      eta_n(idwn.ind, mesh->yend + 1) = eta_n(idwn.ind, mesh->yend);
    }

  }

  /////////////////////////////////////////////////////
  // Neutral density
  ddt(Nn2D) = -Div(Vn2D_contravariant, Nn2D);

  // Particles
  Sn = density_source; // Save for possible output
  if (localstate.isSet("density_source")) {
    Sn += DC(get<Field3D>(localstate["density_source"]));
  }
  ddt(Nn2D) += Sn; // Always add density_source
  
  /////////////////////////////////////////////////////
  // Neutral Velocity  
  // Note: Vn2D.y is proportional to the parallel flow
  //       Vn2D.z is proportional to the toroidal angular momentum
  // Poloidal pressure gradients therefore change the parallel flow
  // without changing the toroidal flow.
  Vector2D GradPn2D = Grad(Pn2D);
  ddt(Vn2D) = GradPn2D / (-AA * Nn2D_floor);
  ASSERT2(ddt(Vn2D).covariant);

  //////////////////////////////////////////////////////
  // Momentum advection

  // Convert to cylindrical coordinates for velocity
  // advection term. This is to avoid Christoffel symbol
  // terms in curvilinear geometry
  Field2D vr = Txr * Vn2D.x + Tyr * Vn2D.y; // Grad R component
  Field2D vz = Txz * Vn2D.x + Tyz * Vn2D.y; // Grad Z component

  // Advect as scalars (no Christoffel symbols needed)
  // X-Y advection
  if (momentum_advection) {
    ddt(vr) = -V_dot_Grad(Vn2D_contravariant, vr);
    ddt(vz) = -V_dot_Grad(Vn2D_contravariant, vz);
  } else {
    ddt(vr) = 0.0;
    ddt(vz) = 0.0;
  }

  // Viscosity
  if (constant_transport_coef) {
    ddt(vr) += Laplace_FV(neutral_viscosity, vr);
    ddt(vz) += Laplace_FV(neutral_viscosity, vz);
  } else {
    ddt(vr) += Laplace_FV(eta_n, vr) / (AA * Nn2D_floor);
    ddt(vz) += Laplace_FV(eta_n, vz) / (AA * Nn2D_floor);
  }


  // Convert back to field-aligned coordinates
  ddt(Vn2D).x += Urx * ddt(vr) + Uzx * ddt(vz);
  ddt(Vn2D).y += Ury * ddt(vr) + Uzy * ddt(vz);

  if (toroidal_flow) {
    Field2D vphi = Vn2D.z / (sigma_Bp * Rxy); // Toroidal component

    if (momentum_advection) {
      ddt(vphi) = -V_dot_Grad(Vn2D_contravariant, vphi);

      if (curved_torus) {
        // Toroidal advection transforms toroidal flow into radial
        ddt(vphi) -= vphi * vr / Rxy; // Conservation of angular momentum
        ddt(vr) += vphi * vphi / Rxy; // Centrifugal force
      }
    } else {
      ddt(vphi) = 0.0;
    }

    // Viscosity
    if (constant_transport_coef) {
      ddt(vphi) += Laplace_FV(neutral_viscosity, vphi);
    } else {
      ddt(vphi) += Laplace_FV(eta_n, vphi) / (AA * Nn2D_floor);
    }

    // Convert back to field-aligned
    ddt(Vn2D).z += ddt(vphi) * sigma_Bp * Rxy;
  }

  /////////////////////////////////////////////////////
  // Atomic processes

  // Momentum. Note need to turn back into covariant form
  // The parallel and toroidal components partly cancel when converting
  // velocity to contravariant form for advection in X-Y.
  if (localstate.isSet("momentum_source")) {
    Snv = DC(get<Field3D>(localstate["momentum_source"]));
    Field2D Fpar_mN = Snv / (AA * Nn2D_floor);
    ddt(Vn2D).y += Fpar_mN * (coord->J * coord->Bxy); // Parallel flow

    if (toroidal_flow) {
      ddt(Vn2D).z += Fpar_mN * coord->g_23 / (coord->J * coord->Bxy); // Toroidal flow
    }

    // NOTE: Should we add the contribution of Sn here? 
    // Sn is introduced in the momentum equation
    // because we solve for Vn instead of AA*Nn*Vn
    // It is propably something like that: 
    // ddt(Vn2D).y += Vn2D.y * Sn / Nn2D_floor * (coord->J * coord->Bxy); // Parallel flow
    // if (toroidal_flow) {
    //   ddt(Vn2D).z += Vn2D.z * Sn / Nn2D_floor * coord->g_23 / (coord->J * coord->Bxy); // Toroidal flow
    // }   

  } else {
    Snv = 0;
  }

  if (localstate.isSet("collision_frequency")) {
    // Damp flow perpendicular to B due to collisions
    Field2D collision_freq = DC(get<Field3D>(localstate["collision_frequency"]));
    // Radial flow
    ddt(Vn2D).x -= Vn2D.x * collision_freq;
    // Binormal flow
    ddt(Vn2D).z -= (Vn2D.z - (coord->g_23 / coord->g_22) * Vn2D.y) * collision_freq;
  }  


  //////////////////////////////////////////////////////
  // Neutral pressure
  ddt(Pn2D) = -adiabatic_index * Div(Vn2D, Pn2D)
              + (adiabatic_index - 1.) * (Vn2D_contravariant * GradPn2D);
  
  if (constant_transport_coef) {
    ddt(Pn2D) += Laplace_FV(Nn2D_floor * neutral_conduction, Tn2D);
  } else {
    ddt(Pn2D) += (adiabatic_index - 1.) * Laplace_FV(kappa_n, Tn2D);    
  }

  // Energy source
  Sp = pressure_source;
  if (localstate.isSet("energy_source")) {
    Sp += (adiabatic_index - 1) * DC(get<Field3D>(localstate["energy_source"]));
  }
  ddt(Pn2D) += Sp;
            
  ///////////////////////////////////////////////////////////////////
  // Boundary condition on fluxes

  for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {

    // Loss of thermal energy to the target.
    // This depends on the reflection coefficient
    // and is controlled by the option neutral_gamma
    //         q = neutral_gamma * n * T * cs

    // Density at the target
    const BoutReal Nnout =
        std::max(0.5 * (Nn2D(r.ind, mesh->ystart) + Nn2D(r.ind, mesh->ystart - 1)), 0.0);
    // Temperature at the target
    const BoutReal Tnout =
        std::max(0.5 * (Tn2D(r.ind, mesh->ystart) + Tn2D(r.ind, mesh->ystart - 1)), 0.0);

    // gamma * n * T * cs
    BoutReal q = neutral_gamma * Nnout * Tnout * sqrt(Tnout);
    // Multiply by cell area to get power
    BoutReal heatflux =
        q * (coord->J(r.ind, mesh->ystart) + coord->J(r.ind, mesh->ystart - 1))
        / (sqrt(coord->g_22(r.ind, mesh->ystart))
           + sqrt(coord->g_22(r.ind, mesh->ystart - 1)));

    // Divide by volume of cell, and multiply by 2/3 (for a monatomic ideal gas) to get
    // pressure
    ddt(Pn2D)(r.ind, mesh->ystart) -=
        (adiabatic_index - 1.) * heatflux
        / (coord->dy(r.ind, mesh->ystart) * coord->J(r.ind, mesh->ystart));
  }

  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {

    // Loss of thermal energy to the target.
    // This depends on the reflection coefficient
    // and is controlled by the option neutral_gamma
    //         q = neutral_gamma * n * T * cs

    // Density at the target
    const BoutReal Nnout =
        std::max(0.5 * (Nn2D(r.ind, mesh->yend) + Nn2D(r.ind, mesh->yend + 1)), 0.0);
    // Temperature at the target
    const BoutReal Tnout =
        std::max(0.5 * (Tn2D(r.ind, mesh->yend) + Tn2D(r.ind, mesh->yend + 1)), 0.0);

    // gamma * n * T * cs
    BoutReal q = neutral_gamma * Nnout * Tnout * sqrt(Tnout);
    // Multiply by cell area to get power
    BoutReal heatflux = q
                        * (coord->J(r.ind, mesh->yend) + coord->J(r.ind, mesh->yend + 1))
                        / (sqrt(coord->g_22(r.ind, mesh->yend))
                           + sqrt(coord->g_22(r.ind, mesh->yend + 1)));

    // Divide by volume of cell, and multiply by 2/3 (for a monatomic ideal gas) to get
    // pressure
    ddt(Pn2D)(r.ind, mesh->yend) -=
        (adiabatic_index - 1.) * heatflux
        / (coord->dy(r.ind, mesh->yend) * coord->J(r.ind, mesh->yend));
  }



#if CHECKLEVEL >= 2
  for (auto& i : Nn2D.getRegion("RGN_NOBNDRY")) {
    if (!std::isfinite(ddt(Nn2D)[i])) {
      throw BoutException("ddt(N{}) non-finite at {}\n", name, i);
    }
    if (!std::isfinite(ddt(Pn2D)[i])) {
      throw BoutException("ddt(P{}) non-finite at {}\n", name, i);
    }
    if (!std::isfinite(ddt(Vn2D).x[i])) {
      throw BoutException("ddt(V{}.x) non-finite at {}\n", name, i);
    }
    if (!std::isfinite(ddt(Vn2D).y[i])) {
      throw BoutException("ddt(V{}.y) non-finite at {}\n", name, i);
    }
    if (!std::isfinite(ddt(Vn2D).z[i])) {
      throw BoutException("ddt(V{}.z) non-finite at {}\n", name, i);
    }
  }
#endif

  // Convert back to contravariant components v^x, v^y, v^z
  ddt(Vn2D).toContravariant();
  Vn2D.toContravariant();

  // Ste time derivatives to zero
  if (zero_timederivs) {
    Field2D zero {0.0};
    ddt(Nn2D) = zero;
    ddt(Pn2D) = zero;
    ddt(Vn2D) = 0.0;
    return;
  }

}

/// Add extra fields for output, or set attributes e.g docstrings
void NeutralFullVelocity::outputVars(Options& state) {
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto Cs0 = get<BoutReal>(state["Cs0"]);
  const BoutReal Pnorm = SI::qe * Tnorm * Nnorm;

  state[std::string("N") + name].setAttributes({{"time_dimension", "t"},
                                                {"units", "m^-3"},
                                                {"conversion", Nnorm},
                                                {"standard_name", "density"},
                                                {"long_name", name + " number density"},
                                                {"species", name},
                                                {"source", "neutral_full_velocity"}});

  state[std::string("P") + name].setAttributes({{"time_dimension", "t"},
                                                {"units", "Pa"},
                                                {"conversion", Pnorm},
                                                {"standard_name", "pressure"},
                                                {"long_name", name + " pressure"},
                                                {"species", name},
                                                {"source", "neutral_full_velocity"}});

  set_with_attrs(state["Urx"], Urx, {{"source", "neutral_full_velocity"}});
  set_with_attrs(state["Ury"], Ury, {{"source", "neutral_full_velocity"}});
  set_with_attrs(state["Uzx"], Uzx, {{"source", "neutral_full_velocity"}});
  set_with_attrs(state["Uzy"], Uzy, {{"source", "neutral_full_velocity"}});
  set_with_attrs(state["Txr"], Txr, {{"source", "neutral_full_velocity"}});
  set_with_attrs(state["Txz"], Txz, {{"source", "neutral_full_velocity"}});
  set_with_attrs(state["Tyr"], Tyr, {{"source", "neutral_full_velocity"}});
  set_with_attrs(state["Tyz"], Tyz, {{"source", "neutral_full_velocity"}});

  if (output_ddt) {
    set_with_attrs(
        state[std::string("ddt(N") + name + std::string(")")], ddt(Nn2D),
        {{"time_dimension", "t"},
         {"units", "m^-3 s^-1"},
         {"conversion", Nnorm * Omega_ci},
         {"long_name", std::string("Rate of change of ") + name + " number density"},
         {"source", "neutral_full_velocity"}});
    set_with_attrs(state[std::string("ddt(P") + name + std::string(")")], ddt(Pn2D),
                   {{"time_dimension", "t"},
                    {"units", "Pa s^-1"},
                    {"conversion", Pnorm * Omega_ci},
                    {"source", "neutral_full_velocity"}});
    // set_with_attrs(state[std::string("ddt(NV") + name + std::string(")")], ddt(Vn2D),
    //                {{"time_dimension", "t"},
    //                 {"units", "kg m^-2 s^-2"},
    //                 {"conversion", SI::Mp * Nnorm * Cs0 * Omega_ci},
    //                 {"source", "neutral_full_velocity"}});
  }  
  if (diagnose) {
    set_with_attrs(state[std::string("T") + name], Tn2D,
                   {{"time_dimension", "t"},
                    {"units", "eV"},
                    {"conversion", Tnorm},
                    {"standard_name", "temperature"},
                    {"long_name", name + " temperature"},
                    {"species", name},
                    {"source", "neutral_full_velocity"}});

    set_with_attrs(state[std::string("V") + name], Vnpar,
                   {{"time_dimension", "t"},
                    {"units", "m / s"},
                    {"conversion", Cs0},
                    {"standard_name", "velocity"},
                    {"long_name", name + " parallel velocity"},
                    {"species", name},
                    {"source", "neutral_full_velocity"}});

    if (Dnn.isAllocated()) {
      set_with_attrs(state[std::string("Dnn") + name], Dnn,
                    {{"time_dimension", "t"},
                      {"units", "m^2/s"},
                      {"conversion", Cs0 * Cs0 / Omega_ci},
                      {"standard_name", "diffusion coefficient"},
                      {"long_name", name + " diffusion coefficient"},
                      {"source", "neutral_full_velocity"}});
    }
    set_with_attrs(state[std::string("SN") + name], Sn,
                   {{"time_dimension", "t"},
                    {"units", "m^-3 s^-1"},
                    {"conversion", Nnorm * Omega_ci},
                    {"standard_name", "density source"},
                    {"long_name", name + " number density source"},
                    {"source", "neutral_full_velocity"}});
    set_with_attrs(state[std::string("SP") + name], Sp,
                   {{"time_dimension", "t"},
                    {"units", "Pa s^-1"},
                    {"conversion", SI::qe * Tnorm * Nnorm * Omega_ci},
                    {"standard_name", "pressure source"},
                    {"long_name", name + " pressure source"},
                    {"source", "neutral_full_velocity"}});
    set_with_attrs(state[std::string("SNV") + name], Snv,
                   {{"time_dimension", "t"},
                    {"units", "kg m^-2 s^-2"},
                    {"conversion", SI::Mp * Nnorm * Cs0 * Omega_ci},
                    {"standard_name", "momentum source"},
                    {"long_name", name + " momentum source"},
                    {"source", "neutral_full_velocity"}});
    set_with_attrs(state[std::string("S") + name + std::string("_src")], density_source,
                   {{"time_dimension", "t"},
                    {"units", "m^-3 s^-1"},
                    {"conversion", Nnorm * Omega_ci},
                    {"standard_name", "density source"},
                    {"long_name", name + " number density source"},
                    {"species", name},
                    {"source", "neutral_full_velocity"}});
    set_with_attrs(state[std::string("P") + name + std::string("_src")], pressure_source,
                   {{"time_dimension", "t"},
                    {"units", "Pa s^-1"},
                    {"conversion", Pnorm * Omega_ci},
                    {"standard_name", "pressure source"},
                    {"long_name", name + " pressure source"},
                    {"species", name},
                    {"source", "neutral_full_velocity"}});
  }
}
