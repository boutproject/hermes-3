
#include <bout/constants.hxx>
#include <bout/derivs.hxx>
#include <bout/difops.hxx>
#include <bout/fv_ops.hxx>
#include <bout/output_bout_types.hxx>

#include "../include/div_ops.hxx"
#include "../include/hermes_build_config.hxx"
#include "../include/neutral_full_velocity_curv.hxx"

using bout::globals::mesh;

using ParLimiter = FV::Upwind;


NeutralFullVelocityCurv::NeutralFullVelocityCurv(const std::string& name, Options& alloptions, Solver* solver)
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
  yboundary.init(options);

  AA = options["AA"].doc("Particle atomic mass. Proton = 1").withDefault(1.0);
  
  // Evolving variables e.g name is "h" or "h+"
  solver->add(Nn, std::string("N") + name);


  evolve_momentum = options["evolve_momentum"]
                        .doc("Evolve parallel neutral momentum?")
                        .withDefault<bool>(true);

  evolve_momentum_xz = options["evolve_momentum_xz"]
                        .doc("Evolve poloidal neutral momentum?")
                        .withDefault<bool>(true);

  momentum_advection = options["momentum_advection"]
                        .doc("Use non-linear advection of momentum?")
                        .withDefault<bool>(false);

  momentum_loss = options["momentum_loss"]
                        .doc("Use perpendicular momentu loss by CX with plasma, assuming plasma has no perpendicular velocity?")
                        .withDefault<bool>(false);
  
  evolve_pressure = options["evolve_pressure"]
                        .doc("Evolve neutral pressure?")
                        .withDefault<bool>(false);

  inherited_T = options["inherited_T"]
                        .doc("Inherit temperature from h+?")
                        .withDefault<bool>(false);
  
  isMMS = options["isMMS"]
                        .doc("Is this MMS? If yes, stop sources and sinks")
                        .withDefault<bool>(false);

  if (evolve_momentum) {
    solver->add(NVn, std::string("NV") + name);
  } else {
    output_warn.write(
        "WARNING: Not evolving neutral parallel momentum. NVn and Vn set to zero\n");
    initial_Vn = options["initial_Vn"]
                        .doc("Initial neutral velocity in parallel direction?")
                        .withDefault(Field3D{0.0}) / (meters/seconds);
    Vn = initial_Vn;
    NVn = AA * Nn * Vn;

  }



  if (evolve_momentum_xz) {
    solver->add(NVn_x, std::string("NVx") + name);
    solver->add(NVn_z, std::string("NVz") + name);
  } else {
    output_warn.write(
        "WARNING: Not evolving neutral parallelpoloidal momentum. Set to zero\n");
    initial_Vn_x = options["initial_Vn_x"]
                        .doc("Initial neutral velocity in x?")
                        .withDefault(Field3D{0.0}) / (meters/seconds);

    initial_Vn_z = options["initial_Vn_z"]
                        .doc("Initial neutral velocity in z?")
                        .withDefault(Field3D{0.0}) / (meters/seconds);
    

    Vn_x = initial_Vn_x;
    Vn_z = initial_Vn_z;
    NVn_x = AA * Nn * Vn_x;
    NVn_z = AA * Nn * Vn_z;
    
  }

  



  if (evolve_pressure) {
    solver->add(Pn, std::string("P") + name);
  } else {
    output_warn.write(
        "WARNING: Not evolving neutral pressure!");

    initial_Tn = options["initial_Tn"]
                        .doc("Initial neutral temperature when pressure is not evolved?")
			.withDefault(Field3D{1.0}) / Tnorm;

    Tn = initial_Tn;
    Pn = Tn * Nn;
  }


  
  
  density_floor = options["density_floor"]
                 .doc("A minimum density used when dividing NVn by Nn. "
                      "Normalised units.")
                 .withDefault(1e-8);

  dissipative = options["dissipative"]
                 .doc("Use strong dissipation in parallel divergence?")
                 .withDefault(true);

  disable_dndt = options["disable_dndt"]
                 .doc("Disable dn/dt?")
                 .withDefault(false);
  
  use_finite_difference = options["use_finite_difference"]
                   .doc("Use finite difference for perpendicular diffusion?")
                   .withDefault<bool>(false);


  parallel_dirichlet = options["parallel_dirichlet"]
                   .doc("Use parallel dirichlet boundary conditions for the plasma?")
                   .withDefault<bool>(true);
  

  
  n_lowsource = options["n_lowsource"].withDefault(-1.0) / Nnorm;
  T_lowsource = options["T_lowsource"].withDefault(-1.0) / Tnorm;
  lowsource_scale = options["lowsource_scale"].withDefault(1e-5) * Omega_ci;
  
  neutral_lmax = options["neutral_lmax"].doc("Largest distance to the target, limits diffusion").withDefault<BoutReal>(0.1) / meters;
  
  temperature_floor = options["temperature_floor"].doc("Low temperature scale for low_T_diffuse_perp")
    .withDefault<BoutReal>(0.1) / get<BoutReal>(alloptions["units"]["eV"]);
  
  pressure_floor = density_floor * temperature_floor;

  precondition = options["precondition"]
                     .doc("Enable preconditioning in neutral model?")
                     .withDefault<bool>(false);

  lax_flux = options["lax_flux"]
                     .doc("Enable stabilising lax flux?")
                     .withDefault<bool>(true);


  neutral_viscosity = options["neutral_viscosity"]
                          .doc("Include neutral gas viscosity?")
                          .withDefault<bool>(false);

  neutral_conduction = options["neutral_conduction"]
                          .doc("Include neutral gas heat conduction?")
                          .withDefault<bool>(false);

  flux_limit = options["flux_limit"].doc("Limit the flux by taking harmonic mean of velocity and sound speed?").withDefault<BoutReal>(-1.0);
  
  BoutReal diffusion_norm = meters * meters / seconds;

  include_D = options.isSet("anomalous_D");
  anomalous_D = options["anomalous_D"].withDefault(Field3D{0.0}) / diffusion_norm;

  include_nu = options.isSet("anomalous_nu");
  anomalous_nu = options["anomalous_nu"].withDefault(Field3D{0.0}) / diffusion_norm;

  anomalous_D.applyBoundary("neumann");
  anomalous_nu.applyBoundary("neumann");
  mesh->communicate(anomalous_D, anomalous_nu);
  anomalous_D.applyParallelBoundary("parallel_neumann_o1");
  anomalous_nu.applyParallelBoundary("parallel_neumann_o1");
  
  if (precondition) {
    inv = Laplacian::create(&options["precon_laplace"]);
    inv->setCoefA(1.0);
  }

  // Optionally output time derivatives
  output_ddt =
      options["output_ddt"].doc("Save derivatives to output?").withDefault<bool>(false);

  diagnose =
      options["diagnose"].doc("Save additional diagnostics?").withDefault<bool>(false);


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


  if (Nn.isFci()) {
    dagp = FCI::getDagp_fv(alloptions, mesh);
  }
}

void NeutralFullVelocityCurv::transform(Options& state) {
  AUTO_TRACE();

  Nn.applyBoundary();
  mesh->communicate(Nn);
  Nn.applyParallelBoundary();

  if (!evolve_momentum) {
    NVn = AA * initial_Vn * Nn;
  }

  NVn.applyBoundary();
  mesh->communicate(NVn);
  NVn.applyParallelBoundary();

  if (!evolve_momentum_xz) {
    NVn_x = AA * initial_Vn_x * Nn;
    NVn_z = AA * initial_Vn_z * Nn;
  }

  NVn_x.applyBoundary();
  NVn_z.applyBoundary();
  mesh->communicate(NVn_x, NVn_z);
  NVn_x.applyParallelBoundary();
  NVn_z.applyParallelBoundary();
  
  if (!evolve_pressure) {
    if (inherited_T) {
      Options& allspecies = state["species"];
      Options& donor_species = allspecies["h+"];
      const auto donor_T = GET_NOBOUNDARY(Field3D, donor_species["temperature"]);
      Pn = donor_T * Nn;
    } else {
      Pn = initial_Tn * Nn;
    }
  }

  Pn.applyBoundary();
  mesh->communicate(Pn);
  Pn.applyParallelBoundary();
  
  Nn = floor(Nn, 1e-3 * density_floor);
  Pn = floor(Pn, 1e-3 * pressure_floor);

  Nnlim = floor(Nn, density_floor);
  Pnlim = floor(Pn, pressure_floor);  
  
  Vn = NVn / (AA * Nnlim);
  Vn_x = NVn_x / (AA * Nnlim);
  Vn_z = NVn_z / (AA * Nnlim);


  Vnlim = 1.0 * Vn;
  Vn_xlim = 1.0 * Vn_x;
  Vn_zlim = 1.0 * Vn_z;


  
  Tn = Pn / Nnlim;

  /////////////////////////////////////////////////////
  // Parallel boundary conditions
  if (!isMMS && parallel_dirichlet) {
    TRACE("Neutral boundary conditions");
    yboundary.iter_pnts([&](auto& pnt) {
      // Free boundary (constant gradient) density
      pnt.dirichlet_o2(Nn, pnt.extrapolate_sheath_o2(Nn));

      // Zero gradient temperature, heat flux added later
      pnt.neumann_o2(Tn,0.0);

      // Zero-gradient pressure
      pnt.neumann_o1(Pn,0.0);
      pnt.neumann_o1(Pnlim,0.0);
    
      // No flow into wall
      pnt.dirichlet_o2(Vn,0.0); 
      pnt.dirichlet_o2(NVn,0.0);
    
    }); // end yboundary.iter_pnts()
  } 


  // Set values in the state
  auto& localstate = state["species"][name];
  set(localstate["density"], Nn);
  set(localstate["AA"], AA); // Atomic mass
  set(localstate["pressure"], Pn);
  set(localstate["momentum"], NVn);
  set(localstate["velocity"], Vn);
  set(localstate["temperature"], Tn);
  
  set(localstate["momentum_x"], NVn_x);
  set(localstate["momentum_z"], NVn_z);
  set(localstate["velocity_x"], Vn_x);
  set(localstate["velocity_z"], Vn_z);
}

void NeutralFullVelocityCurv::finally(const Options& state) {
  AUTO_TRACE();
  auto& localstate = state["species"][name];

  logPnlim = log(Pnlim);

  Field3D Tnlim = floor(Tn, temperature_floor);
  
  sound_speed = 0;
  if (lax_flux) {
    sound_speed = sqrt(Tnlim * (5. / 3) / AA);
  }

  if (isMMS) {
    sound_speed = 0.0;
  }

  if (flux_limit > 0.0) {
    BOUT_FOR(i, Nn.getRegion("RGN_NOY")) {
      const auto iyp = i.yp();
      const auto iym = i.ym();

      BoutReal vmax = flux_limit * sound_speed[i];
      
      if (Vn.yup()[iyp] > 0.0) {
	Vnlim.yup()[iyp] = Vnlim.yup()[iyp] * vmax / (Vnlim.yup()[iyp] + vmax);
      } else {
	Vnlim.yup()[iyp] = -fabs(Vnlim.yup()[iyp]) * vmax / (fabs(Vnlim.yup()[iyp]) + vmax);
      }

      if (Vn.ydown()[iym] > 0.0) {
	Vnlim.ydown()[iym] = Vnlim.ydown()[iym] * vmax / (Vnlim.ydown()[iym] + vmax);
      } else {
        Vnlim.ydown()[iym] = -fabs(Vnlim.ydown()[iym]) * vmax / (fabs(Vnlim.ydown()[iym]) + vmax);
      }

      if (Vn_x[i] > 0.0) {
	Vn_xlim[i] = Vn_x[i] * vmax / (Vn_x[i] + vmax);
      } else {
	Vn_xlim[i] = -fabs(Vn_x[i]) * vmax / (fabs(Vn_x[i]) + vmax);
      }

      if (Vn_z[i] > 0.0) {
	Vn_zlim[i] = Vn_z[i] * vmax / (Vn_z[i] + vmax);
      }	else {
	Vn_zlim[i] = -fabs(Vn_z[i]) * vmax / (fabs(Vn_z[i]) + vmax);
      }
      
      
      
    }
      
  }

  
  Field3D Rnn = sqrt(Tnlim / AA) / neutral_lmax;

  if (localstate.isSet("collision_frequency")) {
    Dnn = (Tnlim / AA) / (get<Field3D>(localstate["collision_frequency"]) + Rnn);
  } else {
    Dnn = (Tnlim / AA) / Rnn;
  }
  
  if (isMMS) {
    Dnn.applyBoundary("free_o2");
  } else {
    Dnn.applyBoundary("neumann");
  }
  
  mesh->communicate(Dnn);
  Dnn.applyParallelBoundary("parallel_neumann_o1");  
  
  kappa_n = (5. / 2) * Dnn*Nn;

  eta_n = AA * (2. / 5) * kappa_n;
  
  
  /////////////////////////////////////////////////////
  // Neutral density
  TRACE("Neutral density");
  Field3D dummy;

  
  // Parallel advection
  ddt(Nn) = -FV::Div_par_mod<hermes::Limiter>(Nn, Vnlim, sound_speed, dummy, dissipative, false);

  
  // Perpendicular advection
  ddt(Nn) += -Div_perp(Nn, Vn_xlim, Vn_zlim, sound_speed, dissipative);  


  // Sources
  Sn = density_source; // Save for possible output
  if (localstate.isSet("density_source")) {
    Sn += get<Field3D>(localstate["density_source"]);
  }
  if (!isMMS) {
    ddt(Nn) += Sn; // Always add density_source
  }


  // Lowsource
  if (n_lowsource > 0.0) {
    ddt(Nn) += low_sourceterm(Nn, n_lowsource, lowsource_scale);
  } 


  // Anomalous diffusion
  if (include_D) {
    ddt(Nn) += Div_a_Grad_perp(anomalous_D, Nn, use_finite_difference);
  }


  // Disabling time derivate for testing purposes
  if (disable_dndt) {
    ddt(Nn) = 0.0;
  }


  /////////////////////////////////////////////////////                                                                                                                                                           
  // Neutral pressure

  
  TRACE("Neutral pressure");

  if (evolve_pressure) {
    Field3D dummy_Pn;
    // Parallel advection                                                                                                                                                                                          
    if (!isMMS) {
      ddt(Pn) = -(5.0 / 3.0) * FV::Div_par_mod<hermes::Limiter>(Pn, Vnlim, sound_speed, dummy_Pn, dissipative);      // Parallel advection                                                                           
    } else {
      ddt(Pn) = -(5.0 / 3.0) * Div_par(Pn * Vnlim);
    }
    ddt(Pn) += (2. / 3) * Vnlim * Grad_par(Pn);


    // Perpendicular advection

    ddt(Pn) -= (5.0 / 3.0) * Div_perp(Pn, Vn_xlim, Vn_zlim, sound_speed, dissipative);
    
    ddt(Pn) += (2.0 / 3.0) * (Vn_xlim * Grad_x(Pn) + Vn_z * Grad_z(Pn));


    // Thermal conduction
    if (neutral_conduction) {
      ddt(Pn) += (2.0/3.0) * Div_par_K_Grad_par_mod(kappa_n, Tn, dummy_Pn, false);
      ddt(Pn) += (2.0/3.0) * Div_a_Grad_perp(kappa_n, Tn, use_finite_difference);
    }

    Sp = pressure_source;
    if (localstate.isSet("energy_source")) {
      Sp += (2. / 3) * get<Field3D>(localstate["energy_source"]);
    }
    
    if (!isMMS) {
      ddt(Pn) += Sp;
    }

    if (T_lowsource > 0.0) {
      ddt(Pn) += low_sourceterm(Tn, T_lowsource, lowsource_scale);
    }
    
    
  }

  
  /////////////////////////////////////////////////////                                                                                                                                                           
  // Neutral parallel momentum 

  if (evolve_momentum) {
    Field3D dummy_NVn;

    
    // Parallel gradient
    ddt(NVn) = -Grad_par(Pn);

    
    // Parallel advection
    if (momentum_advection) {                                                                                                                        
      ddt(NVn) += -AA * FV::Div_par_fvv<hermes::Limiter>(Nnlim, Vnlim, sound_speed);
    }

    
    // Perpendicular advection
    if (include_D) {
      ddt(NVn) += Div_par_K_Grad_par_mod(AA * anomalous_D * Vn, Nn, dummy_NVn, false);
    }
    
    if (include_nu) {
      ddt(NVn) += Div_par_K_Grad_par_mod(AA * anomalous_nu * Nn, Vn, dummy_NVn, false);
    }

    if (neutral_viscosity) {
      Field3D viscosity_source = Div_par_K_Grad_par_mod(eta_n , Vn , mf_visc_par_ylow , false);
      viscosity_source += Div_a_Grad_perp(eta_n , Vn, use_finite_difference);
      ddt(NVn) += viscosity_source;
      if (evolve_pressure) {
	ddt(Pn)  += -(2. /3) * Vn * viscosity_source;
      }
    }
        
  }


  
  /////////////////////////////////////////////////////                                                                                                                                                           
  // Neutral perpendicular momentum 

  if (evolve_momentum_xz) {
    Field3D dummy_NVnxz;


    // Perpendicular gradient
    ddt(NVn_x) = -Grad_x(Pn);
    ddt(NVn_z) = -Grad_z(Pn);
    

    // Perpendicular momentum advection
    if (momentum_advection) {
      ddt(NVn_x) += -AA * Div_perp_flows(Nnlim, Vn_x, Vn_xlim, Vn_zlim, sound_speed);
      ddt(NVn_z) += -AA * Div_perp_flows(Nnlim, Vn_z, Vn_xlim, Vn_zlim, sound_speed);
    }


    // Anomalous diffusion
    if (include_D) {
      ddt(NVn_x) += Div_a_Grad_perp(AA * anomalous_D * Vn_x, Nn, use_finite_difference);
      ddt(NVn_z) += Div_a_Grad_perp(AA * anomalous_D * Vn_z, Nn, use_finite_difference);
    }
    

    // Anomalous viscosity
    if (include_nu) {
      ddt(NVn_x) += Div_a_Grad_perp(AA * anomalous_nu * Nn, Vn_x, use_finite_difference);
      ddt(NVn_z) += Div_a_Grad_perp(AA * anomalous_nu * Nn, Vn_z, use_finite_difference);
    }

    if (momentum_loss && localstate.isSet("collision_frequency") ) {
      const Options& allspecies = state["species"];
      const Options& donor_species = allspecies["h+"];
      const auto N_ion = GET_NOBOUNDARY(Field3D, donor_species["density"]);
      Field3D mom_lossfactor = -AA * N_ion * get<Field3D>(localstate["collision_frequency"]);
      ddt(NVn_x) += mom_lossfactor * Vn_x;
      ddt(NVn_z) += mom_lossfactor * Vn_z;	    
    }
			   
    
    // Collisional viscosity
    if (neutral_viscosity) {
      Field3D viscosity_source_x = Div_par_K_Grad_par_mod(eta_n , Vn_x , dummy_NVnxz, false);
      viscosity_source_x += Div_a_Grad_perp(eta_n , Vn_x, use_finite_difference);
      ddt(NVn_x) += viscosity_source_x;
      if (evolve_pressure) {
        ddt(Pn)  += -(2. /3) * Vn_x * viscosity_source_x;
      }

      Field3D viscosity_source_z = Div_par_K_Grad_par_mod(eta_n , Vn_z , dummy_NVnxz, false);
      viscosity_source_z += Div_a_Grad_perp(eta_n , Vn_z, use_finite_difference);
      ddt(NVn_z) += viscosity_source_z;
      if (evolve_pressure) {
        ddt(Pn)  += -(2. /3) * Vn_z * viscosity_source_z;
      }
    }
  } 

  
  

}

void NeutralFullVelocityCurv::outputVars(Options& state) {
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

  if (diagnose) {
    set_with_attrs(state[std::string("Vx") + name], Vn_x,
		   {{"source", "neutral_full_velocity_curv"}});
    set_with_attrs(state[std::string("Vz") + name], Vn_z,
                   {{"source", "neutral_full_velocity_curv"}});

  
  
  
  }
  
  
}

void NeutralFullVelocityCurv::precon(const Options& state, BoutReal gamma) {
  if (!precondition) {
    return;
  }
  const auto& species = state["species"][name];
  const Field3D N = get<Field3D>(species["density"]);

  // Set the coefficient in Div_par( B * Grad_par )
  Field3D coef = - gamma * Dnn;

  inv->setCoefD(coef);
  Field3D dT = ddt(Pn);
  dT.applyBoundary("neumann");
  mesh->communicate(dT);
  Field3D dummy = 0.0;
  ddt(Pn) = inv->solve(dT, ddt(Pn));

}
