/// Ion viscosity model

#include <bout/constants.hxx>
#include <bout/fv_ops.hxx>
#include <bout/difops.hxx>
#include <bout/output_bout_types.hxx>
#include <bout/mesh.hxx>
#include "../include/hermes_utils.hxx"

#include "../include/ion_viscosity.hxx"
#include "../include/div_ops.hxx"

using bout::globals::mesh;

IonViscosity::IonViscosity(std::string name, Options& alloptions, Solver*) {
  auto& options = alloptions[name];

  eta_limit_alpha = options["eta_limit_alpha"]
    .doc("Viscosity flux limiter coefficient. <0 = turned off")
    .withDefault(-1.0);

  diagnose = options["diagnose"]
      .doc("Output additional diagnostics?")
      .withDefault<bool>(false);

  perpendicular = options["perpendicular"]
    .doc("Include perpendicular flow? (Requires phi)")
    .withDefault<bool>(false);

  viscosity_collisions_mode = options["viscosity_collisions_mode"]
      .doc("Can be multispecies: all collisions, or braginskii: self collisions")
      .withDefault<std::string>("multispecies");

  bounce_frequency = options["bounce_frequency"]
    .doc("Include neoclassical modification of the collision time?")
    .withDefault<bool>(false);

  bounce_frequency_q95 = options["bounce_frequency_q95"]
    .doc("Include input for q95 when using bounce frequency modification to viscosity")
    .withDefault(4.0);

  bounce_frequency_epsilon = options["bounce_frequency_epsilon"]
    .doc("Include input for inverse aspect ratio epsilon when using bounce frequency modification to viscosity")
    .withDefault(0.3);

  bounce_frequency_R = options["bounce_frequency_R"]
    .doc("Include input for major radius R when using bounce frequency modification to viscosity")
    .withDefault(2.0);

  if (perpendicular) {
    // Read curvature vector
    try {
      Curlb_B.covariant = false; // Contravariant
      mesh->get(Curlb_B, "bxcv");
    } catch (BoutException& e) {
      // May be 2D, reading as 3D
      Vector2D curv2d;
      curv2d.covariant = false;
      mesh->get(curv2d, "bxcv");
      Curlb_B = curv2d;
    }

    if (Options::root()["mesh"]["paralleltransform"]["type"].as<std::string>()
        == "shifted") {
      Field2D I;
      mesh->get(I, "sinty");
      Curlb_B.z += I * Curlb_B.x;
    }

    // Normalisations
    const Options& units = alloptions["units"];
    const BoutReal Bnorm = units["Tesla"];
    const BoutReal Lnorm = units["meters"];

    auto coord = mesh->getCoordinates();



    Curlb_B.x /= Bnorm;
    Curlb_B.y *= SQ(Lnorm);
    Curlb_B.z *= SQ(Lnorm);
    
    Curlb_B *= 2. / coord->Bxy;
  }
  if (bounce_frequency) {
    const Options& units = alloptions["units"];
    const BoutReal Lnorm = units["meters"];
    bounce_frequency_R /= Lnorm;
  }
  

}

void IonViscosity::transform(Options &state) {
  AUTO_TRACE();

  Options& allspecies = state["species"];

  auto coord = mesh->getCoordinates();
  const Field2D Bxy = coord->Bxy;
  const Field2D sqrtB = sqrt(Bxy);
  const Field2D Grad_par_logB = Grad_par(log(Bxy));

  // Loop through all species
  for (auto& kv : allspecies.getChildren()) {
    if (kv.first == "e") {
      continue; // Skip electrons -> only ions
    }
    const auto& species_name = kv.first;

    Options& species = allspecies[species_name];

    if (!(isSetFinal(species["pressure"], "ion_viscosity") and
          isSetFinal(species["velocity"], "ion_viscosity") and
          isSetFinal(species["charge"], "ion_viscosity"))) {
      // Species doesn't have a pressure, velocity and charge => Skip
      continue;
    }

    if (std::fabs(get<BoutReal>(species["charge"])) < 1e-3) {
      // No charge
      continue;
    }

    if (collision_names.empty()) {     /// Calculate only once - at the beginning

      if (viscosity_collisions_mode == "braginskii") {
        for (const auto& collision : species["collision_frequencies"].getChildren()) {

          std::string collision_name = collision.second.name();

          if (/// Self-collisions
              (collisionSpeciesMatch(    
                collision_name, species.name(), species.name(), "coll", "exact")) or
              /// Ion-electron collisions
              (collisionSpeciesMatch(    
                collision_name, species.name(), "+", "coll", "partial"))) {
                  
                  collision_names.push_back(collision_name);
                }
        }
      // Multispecies mode: all collisions and CX are included
      } else if (viscosity_collisions_mode == "multispecies") {
        for (const auto& collision : species["collision_frequencies"].getChildren()) {

          std::string collision_name = collision.second.name();

          if (/// Charge exchange
              (collisionSpeciesMatch(    
                collision_name, species.name(), "", "cx", "partial")) or
              /// Any collision (en, in, ee, ii, nn)
              (collisionSpeciesMatch(    
                collision_name, species.name(), "", "coll", "partial"))) {
                  
                  collision_names.push_back(collision_name);
                }
        }
      } else {
        throw BoutException("\tviscosity_collisions_mode for {:s} must be either multispecies or braginskii", species.name());
      }

      if (collision_names.empty()) {
        throw BoutException("\tNo collisions found for {:s} in ion_viscosity for selected collisions mode", species.name());
      }

      /// Write chosen collisions to log file
      output_info.write("\t{:s} viscosity collisionality mode: '{:s}' using ",
                      species.name(), viscosity_collisions_mode);
      for (const auto& collision : collision_names) {        
        output_info.write("{:s} ", collision);
      }

      output_info.write("\n");

      }

    /// Collect the collisionalities based on list of names
    nu = 0;
    for (const auto& collision_name : collision_names) {
      nu += GET_VALUE(Field3D, species["collision_frequencies"][collision_name]);
    }

    const Field3D tau = 1. / nu;
    const Field3D P = get<Field3D>(species["pressure"]);
    const Field3D V = get<Field3D>(species["velocity"]);

    // Parallel ion viscosity (4/3 * 0.96 coefficient)
    Field3D eta = 1.28 * P * tau;
    
    Field2D bounce_factor = 1.0; // if bounce_frequency = false, this factor does nothing to anything
    Field2D nu_star = 1.0;

    if (bounce_frequency) {
      // Need to collect the DC density and temperature, ion collision time and thermal velocity to calculate the bounce frequency, with the equation taken as in Rozhansky 2009. 
      
      const Field2D N_av = DC(get<Field3D>(species["density"]));
      const Field2D T_av = DC(get<Field3D>(species["temperature"]));
      const Field2D tau_av = DC(tau);

      // calculating ion thermal velocity
      BoutReal mass = get<BoutReal>(species["AA"]);
      const Field2D v_thermal = sqrt(2.0 * T_av / mass);

      nu_star *=  (bounce_frequency_R * bounce_frequency_q95) / ( tau_av * pow(bounce_frequency_epsilon, 1.5) * v_thermal);

      bounce_factor *= (1 / (1 + (1./nu_star))) * (1 / (1 + (1. / pow(bounce_frequency_epsilon, 1.5)) * (1./nu_star)));
      eta *= bounce_factor;
    } 
    
    if (eta_limit_alpha > 0.) {
      // SOLPS-style flux limiter
      // Values of alpha ~ 0.5 typically

      const Field3D q_cl = eta * Grad_par(V);   // Collisional value
      const Field3D q_fl = eta_limit_alpha * P; // Flux limit

      eta = eta / (1. + abs(q_cl / q_fl));

      eta.getMesh()->communicate(eta);
      eta.applyBoundary("neumann");
    }
    
    // This term is the parallel flow part of
    // -(2/3) B^(3/2) Grad_par(Pi_ci / B^(3/2))
    const Field3D div_Pi_cipar = sqrtB * FV::Div_par_K_Grad_par(eta / Bxy, sqrtB * V);

    add(species["momentum_source"], div_Pi_cipar);
    subtract(species["energy_source"], V * div_Pi_cipar); // Internal energy

    if (!perpendicular) {
      if (diagnose) {
        const Field2D P_av = DC(P);
        const Field2D tau_av = DC(tau);
        const Field2D V_av = DC(V);

        // Parallel ion stress tensor component, calculated here because before it was only div_Pi_cipar
        Field2D Pi_cipar = -0.96 * P_av * tau_av * bounce_factor *
                          (2. * Grad_par(V_av) + V_av * Grad_par_logB);
        Field2D Pi_ciperp = 0 * Pi_cipar; // Perpendicular components and divergence of current J equal to 0 for only parallel viscosity case
        Field2D DivJ = 0 * Pi_cipar;
        // Find the diagnostics struct for this species
        auto search = diagnostics.find(species_name);
        if (search == diagnostics.end()) {
          // First time, create diagnostic
          diagnostics.emplace(species_name, Diagnostics {Pi_ciperp, Pi_cipar, DivJ, bounce_factor});
        } else {
          // Update diagnostic values
          auto& d = search->second;
          d.Pi_ciperp = Pi_ciperp;
          d.Pi_cipar = Pi_cipar;
          d.DivJ = DivJ;
          d.bounce_factor = bounce_factor;
        }
      }
      continue; // Skip perpendicular flow parts below
    }

    // The following calculation is performed on the axisymmetric components

    // Need electrostatic potential for ExB flow
    const Field2D phi_av = DC(get<Field3D>(state["fields"]["phi"]));

    // Density and temperature needed for diamagnetic flows
    const Field2D N_av = DC(get<Field3D>(species["density"]));
    const Field2D T_av = DC(get<Field3D>(species["temperature"]));
    const Field2D P_av = DC(P);
    const Field2D tau_av = DC(tau);
    const Field2D V_av = DC(V);

    // Parallel ion stress tensor component
    Field2D Pi_cipar = -0.96 * P_av * tau_av * bounce_factor *
                          (2. * Grad_par(V_av) + V_av * Grad_par_logB);
    // Could also be written as:
    // Pi_cipar = -0.96*Pi*tau*2.*Grad_par(sqrt(Bxy)*Vi)/sqrt(Bxy);

    // Perpendicular ion stress tensor
    // 0.96 P tau kappa * (V_E + V_di + 1.61 b x Grad(T)/B )
    // Note: Heat flux terms are neglected for now
    Field2D Pi_ciperp = -0.5 * 0.96 * P_av * tau_av * bounce_factor * 
      (Curlb_B * Grad(phi_av + 1.61 * T_av) - Curlb_B * Grad(P_av) / N_av);

    // Limit size of stress tensor components
    // If the off-diagonal components of the pressure tensor are large compared
    // to the scalar pressure, then the model is probably breaking down.
    // This can happen in very low density regions
    BOUT_FOR(i, Pi_ciperp.getRegion("RGN_NOBNDRY")) {
      if (fabs(Pi_ciperp[i]) > P_av[i]) {
        Pi_ciperp[i] = SIGN(Pi_ciperp[i]) * P_av[i];
      }
      if (fabs(Pi_cipar[i]) > P_av[i]) {
        Pi_cipar[i] = SIGN(Pi_cipar[i]) * P_av[i];
      }
    }

    // Apply perpendicular boundary conditions
    Pi_ciperp.applyBoundary("neumann");
    Pi_cipar.applyBoundary("neumann");

    // Apply parallel boundary conditions
    int jy = 1;
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      Pi_ciperp(r.ind, jy) = Pi_ciperp(r.ind, jy + 1);
      Pi_cipar(r.ind, jy) = Pi_cipar(r.ind, jy + 1);
    }
    jy = mesh->yend + 1;
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      Pi_ciperp(r.ind, jy) = Pi_ciperp(r.ind, jy - 1);
      Pi_cipar(r.ind, jy) = Pi_cipar(r.ind, jy - 1);
    }

    mesh->communicate(Pi_ciperp, Pi_cipar);

    const Field3D div_Pi_ciperp = - (2. / 3) * Grad_par(Pi_ciperp) + Pi_ciperp * Grad_par_logB;
    //const Field3D div_Pi_ciperp = - (2. / 3) * B32 * Grad_par(Pi_ciperp / B32);

    add(species["momentum_source"], div_Pi_ciperp);
    subtract(species["energy_source"], V_av * div_Pi_ciperp);

    // Total scalar ion viscous pressure
    Field2D Pi_ci = Pi_cipar + Pi_ciperp;

#if CHECKLEVEL >= 1
    for (auto& i : N_av.getRegion("RGN_NOBNDRY")) {
      if (!std::isfinite(Pi_cipar[i])) {
        throw BoutException("{} Pi_cipar non-finite at {}.\n", species_name, i);
      }
      if (!std::isfinite(Pi_ciperp[i])) {
        throw BoutException("{} Pi_ciperp non-finite at {}.\n", species_name, i);
      }
    }
#endif

    Pi_ci.applyBoundary("neumann");

    // Divergence of current in vorticity equation
    Field3D DivJ = Div(0.5 * Pi_ci * Curlb_B) -
      Div_n_bxGrad_f_B_XPPM(1. / 3, Pi_ci, false, true);
    add(state["fields"]["DivJextra"], DivJ);

    // Transfer of energy between ion internal energy and ExB flow
    subtract(species["energy_source"],
             Field3D(0.5 * Pi_ci * Curlb_B * Grad(phi_av + P_av)
                     - (1 / 3) * bracket(Pi_ci, phi_av + P_av, BRACKET_STD)));

    if (diagnose) {
      // Find the diagnostics struct for this species
      auto search = diagnostics.find(species_name);
      if (search == diagnostics.end()) {
        // First time, create diagnostic
        diagnostics.emplace(species_name, Diagnostics {Pi_ciperp, Pi_cipar, DivJ, bounce_factor, nu_star});
      } else {
        // Update diagnostic values
        auto& d = search->second;
        d.Pi_ciperp = Pi_ciperp;
        d.Pi_cipar = Pi_cipar;
        d.DivJ = DivJ;
        d.bounce_factor = bounce_factor;
        d.nu_star = nu_star;
      }
    }
  }
}

void IonViscosity::outputVars(Options &state) {
  AUTO_TRACE();

  if (diagnose) {
    // Normalisations
    auto Nnorm = get<BoutReal>(state["Nnorm"]);
    auto Tnorm = get<BoutReal>(state["Tnorm"]);
    auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
    BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

    for (const auto& it : diagnostics) {
      const std::string& species_name = it.first;
      const auto& d = it.second;

      set_with_attrs(state[std::string("P") + species_name + std::string("_ciperp")], d.Pi_ciperp,
                     {{"time_dimension", "t"},
                      {"units", "Pa"},
                      {"conversion", Pnorm},
                      {"standard_name", "Viscous pressure"},
                      {"long_name", species_name + " perpendicular collisional viscous pressure"},
                      {"species", species_name},
                      {"source", "ion_viscosity"}});

      set_with_attrs(state[std::string("P") + species_name + std::string("_cipar")], d.Pi_cipar,
                     {{"time_dimension", "t"},
                      {"units", "Pa"},
                      {"conversion", Pnorm},
                      {"standard_name", "Viscous pressure"},
                      {"long_name", species_name + " parallel collisional viscous pressure"},
                      {"species", species_name},
                      {"source", "ion_viscosity"}});

      set_with_attrs(state[std::string("DivJvis_") + species_name], d.DivJ,
                     {{"time_dimension", "t"},
                      {"units", "A m^-3"},
                      {"conversion", SI::qe * Nnorm * Omega_ci},
                      {"long_name", std::string("Divergence of viscous current due to species") + species_name},
                      {"species", species_name},
                      {"source", "ion_viscosity"}});

      set_with_attrs(state[std::string("bounce_factor_") + species_name], d.bounce_factor,
                     {{"time_dimension", "t"},
                      {"units", "none"},
                      {"conversion", 1},
                      {"long_name", std::string("Bounce factor for viscosity calculation") + species_name},
                      {"species", species_name},
                      {"source", "ion_viscosity"}});

      set_with_attrs(state[std::string("nu_star_") + species_name], d.nu_star,
                    {{"time_dimension", "t"},
                      {"units", "none"},
                      {"conversion", 1},
                      {"long_name", std::string("nu star") + species_name},
                      {"species", species_name},
                      {"source", "ion_viscosity"}});
    }
  }
}
