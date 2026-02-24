#include <bout/fv_ops.hxx>
#include <bout/solver.hxx>

using bout::globals::mesh;
#include <bout/constants.hxx>
#include "../include/div_ops.hxx"
#include "../include/relax_potential.hxx"
#include "../include/hermes_build_config.hxx"

RelaxPotential::RelaxPotential(std::string name, Options& alloptions, Solver* solver) {
  AUTO_TRACE();

  auto* coord = mesh->getCoordinates();

  // Normalisations
  const Options& units = alloptions["units"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();
  const BoutReal Bnorm = units["Tesla"];
  const BoutReal Lnorm = units["meters"];
  auto& options = alloptions[name];

  exb_advection = options["exb_advection"]
                      .doc("Include nonlinear ExB advection?")
                      .withDefault<bool>(true);

  scale_ExB = options["scale_ExB"]
                   .doc("Scale ExB flow?")
                   .withDefault<BoutReal>(1.0);
  
  diamagnetic =
      options["diamagnetic"].doc("Include diamagnetic current?").withDefault<bool>(true);

  diamagnetic_polarisation =
      options["diamagnetic_polarisation"]
          .doc("Include diamagnetic drift in polarisation current?")
          .withDefault<bool>(true);

  diamagnetic_bracketform = options["diamagnetic_bracketform"]
                        .doc("Include diamagnetic form that uses arakawa brackets? FCI version")
                        .withDefault<bool>(mesh->isFci());
  
  boussinesq = options["boussinesq"]
                   .doc("Use the Boussinesq approximation?")
                   .withDefault<bool>(true);

  floating_boundary = options["floating_boundary"]
                   .doc("Calculate boundary condition for phi from electron temperature?")
                   .withDefault<bool>(false);

  disable_ddt = options["disable_ddt"]
                   .doc("Disable timevolution to iterate")
                   .withDefault<bool>(false);

  viscosity_core = options["viscosity_core"]
    .doc("Kinematic viscosity [m^2/s]")
    .withDefault<Field3D>(0.0)
    / (Lnorm * Lnorm * Omega_ci);

  viscosity_core.applyBoundary("neumann");
  mesh->communicate(viscosity_core);
  viscosity_core.applyParallelBoundary("parallel_neumann_o1");
  
  viscosity = options["viscosity"]
    .doc("Kinematic viscosity [m^2/s]")
    .withDefault<Field3D>(0.0)
    / (Lnorm * Lnorm * Omega_ci);
  
  viscosity.applyBoundary("neumann");
  mesh->communicate(viscosity);
  viscosity.applyParallelBoundary("parallel_neumann_o1");

  viscosity_par = options["viscosity_par"]
    .doc("Parallel kinematic viscosity [m^2/s]")
    .withDefault<Field3D>(0.0)
    / (Lnorm * Lnorm * Omega_ci);

  viscosity_par.applyBoundary("neumann");
  mesh->communicate(viscosity_par);
  viscosity_par.applyParallelBoundary("parallel_neumann_o1");
  
  phi_dissipation = options["phi_dissipation"]
                        .doc("Parallel dissipation of potential [Recommended]")
                        .withDefault<bool>(true);

  average_atomic_mass = options["average_atomic_mass"]
                            .doc("Weighted average atomic mass, for polarisaion current "
                                 "(Boussinesq approximation)")
                            .withDefault<BoutReal>(2.0); // Deuterium

  poloidal_flows =
      options["poloidal_flows"].doc("Include poloidal ExB flow").withDefault<bool>(true);

  lambda_1 = options["lambda_1"].doc("λ_1 > 1").withDefault<BoutReal>(100.0);
  lambda_2 = options["lambda_2"].doc("λ_2 > λ_1").withDefault<BoutReal>(10000.0);

  vort_timedissipation = options["vort_timedissipation"].doc("Vort dissipation").withDefault<BoutReal>(-1.0) * Omega_ci;
  
  solver->add(Vort, "Vort"); // Vorticity evolving
  solver->add(phi1, "phi1"); // Evolving scaled potential ϕ_1 = λ_2 ϕ

  if (diamagnetic and !diamagnetic_bracketform) {
    // Read curvature vector
    try {
      Curlb_B.covariant = false; // Contravariant
      mesh->get(Curlb_B, "bxcv");

    } catch (BoutException& e) {
      // May be 2D, reading as 3D
      Vector2D curv2d;
      curv2d.covariant = false;
      if (mesh->get(curv2d, "bxcv")) {
        throw BoutException("Curvature vector not found in input");
      }
      Curlb_B = curv2d;
    }

    if (Options::root()["mesh"]["paralleltransform"]["type"].as<std::string>()
        == "shifted") {
      Field2D I;
      mesh->get(I, "sinty");
      Curlb_B.z += I * Curlb_B.x;
    }

    Options& units = alloptions["units"];
    BoutReal Bnorm = units["Tesla"];
    BoutReal Lnorm = units["meters"];

    if (mesh->isFci()) {
      // All coordinates (x,y,z) are dimensionless
      // -> e_x has dimensions of length
      Curlb_B.x *= SQ(Lnorm);
    } else {
      // Field-aligned (Clebsch) coordinates
      Curlb_B.x /= Bnorm;
    }

    Curlb_B.y *= SQ(Lnorm);
    Curlb_B.z *= SQ(Lnorm);

    Curlb_B *= 2. / coord->Bxy;
  }

  logB = log(coord->Bxy);
  logB.applyBoundary("neumann");
  mesh->communicate(logB);
  logB.applyParallelBoundary("parallel_neumann_o1");
  
  Bsq = SQ(coord->Bxy);
  if (Vort.isFci()) {
    dagp = FCI::getDagp_fv(alloptions, mesh);

    const auto coord = mesh->getCoordinates();
    // Note: This is 1 for a Clebsch coordinate system
    //       Remove parallel slices before operations
    bracket_factor = sqrt(coord->g_22.withoutParallelSlices()) / (coord->J.withoutParallelSlices() * coord->Bxy);
  } else {
    bracket_factor = 1.0;
  }

  ones = 1.0;
  ones.applyBoundary("neumann");
  mesh->communicate(ones);
  ones.applyParallelBoundary("parallel_neumann_o1");
  
}

void RelaxPotential::transform(Options& state) {
  AUTO_TRACE();

  // Scale potential
  
  phi1.applyBoundary("neumann");
  if (floating_boundary) {
    Field3D Te = GET_NOBOUNDARY(Field3D, state["species"]["e"]["temperature"]);
    BoutReal Me_Mp = get<BoutReal>(state["species"]["e"]["AA"]);
    BoutReal sheathmult = log(0.5 * sqrt(1. / (Me_Mp * PI)));

    if (mesh->lastX()) {
      int n = mesh->LocalNx;
      for (int j = mesh->ystart; j <= mesh->yend; j++) {
	for (int k = 0; k < mesh->LocalNz; k++) {
	  BoutReal phi1value = lambda_2 * sheathmult * 0.5 * (Te(n-2, j, k) + Te(n-3, j, k));
	  phi1(n-2, j, k) = 2.0 * phi1value - phi1(n-3, k, k);
	  phi1(n-1, j, k) = phi1(n-2, j, k);
	}
      }
    }    
  }
  
  
  Vort.applyBoundary();

  mesh->communicate(Vort, phi1);

  if (phi.isFci()){
    phi1.applyParallelBoundary("parallel_neumann_o1");
    Vort.applyParallelBoundary("parallel_neumann_o1");
  }
  auto& fields = state["fields"];

  phi = phi1 / lambda_2;
  
  ddt(Vort) = 0.0;

  if (diamagnetic) {
    // Diamagnetic current. This is calculated here so that the energy sources/sinks
    // can be calculated for the evolving species.

    if (!diamagnetic_bracketform){
    
      Vector3D Jdia;
      Jdia.x = 0.0;
      Jdia.y = 0.0;
      Jdia.z = 0.0;
      Jdia.covariant = Curlb_B.covariant;

      Options& allspecies = state["species"];
      
      // Pre-calculate this rather than calculate for each species
      Vector3D Grad_phi = Grad(phi);
      
      for (auto& kv : allspecies.getChildren()) {
	Options& species = allspecies[kv.first]; // Note: need non-const
	
	if (!(IS_SET_NOBOUNDARY(species["pressure"]) and IS_SET(species["charge"])
	      and (get<BoutReal>(species["charge"]) != 0.0))) {
	  continue; // No pressure or charge -> no diamagnetic current
	}
	// Note that the species must have a charge, but charge is not used,
	// because it cancels out in the expression for current

	auto P = GET_NOBOUNDARY(Field3D, species["pressure"]);

	Vector3D Jdia_species = P * Curlb_B; // Diamagnetic current for this species

	// This term energetically balances diamagnetic term
	// in the vorticity equation
	subtract(species["energy_source"], Jdia_species * Grad_phi);

	Jdia += Jdia_species; // Collect total diamagnetic current
      }
      
      // Note: This term is central differencing so that it balances
      // the corresponding compression term in the species pressure equations
      if (phi.isFci()) {
	mesh->communicate(Jdia);
	Jdia.applyBoundary("neumann");
	Jdia.y.applyParallelBoundary("parallel_neumann_o2");
      }
      Field3D DivJdia = Div(Jdia);
      ddt(Vort) += DivJdia;
      
      if (diamagnetic_polarisation) {
	// Calculate energy exchange term nonlinear in pressure
	// ddt(Pi) += Pi * Div((Pe + Pi) * Curlb_B);
	for (auto& kv : allspecies.getChildren()) {
	  Options& species = allspecies[kv.first]; // Note: need non-const
	  
	  if (!(IS_SET_NOBOUNDARY(species["pressure"]) and IS_SET(species["charge"])
		and IS_SET(species["AA"]))) {
	    continue; // No pressure, charge or mass -> no polarisation current due to
	    // rate of change of diamagnetic flow
	  }
	  auto P = GET_NOBOUNDARY(Field3D, species["pressure"]);

	  add(species["energy_source"], (3. / 2) * P * DivJdia);
	}
      }
      
      set(fields["DivJdia"], DivJdia);
    
    } else {                    // use diamagnetic bracketform
      Options& allspecies = state["species"];
      
      for (auto& kv : allspecies.getChildren()) {
	Options& species = allspecies[kv.first]; // Note: need non-const
	
	if (!(IS_SET_NOBOUNDARY(species["pressure"]) and IS_SET(species["charge"]))) {
	  continue; // No pressure or charge -> no diamagnetic current                                                                                
	}
	
	if (fabs(get<BoutReal>(species["charge"])) < 1e-5) {
	  // No charge                                                                                                                                
	  continue;
	}

	auto P = GET_NOBOUNDARY(Field3D, species["pressure"]);
	Field3D DivJdia_species = 2.0 * bracket(logB, P, BRACKET_ARAKAWA) * bracket_factor;
	ddt(Vort) += DivJdia_species;
	// Balance of this term in the species energy equation
	
	if (diamagnetic_polarisation){
	  add(species["energy_source"], P * DivJdia_species);
	}
	subtract(species["energy_source"], P * 2.0 * bracket(logB, phi) * bracket_factor);
	
      }
    }
  }
  set(fields["vorticity"], Vort);
  set(fields["phi"], phi);
}

void RelaxPotential::finally(const Options& state) {
  AUTO_TRACE();

  const Options& allspecies = state["species"];

  phi = get<Field3D>(state["fields"]["phi"]);
  Vort = get<Field3D>(state["fields"]["vorticity"]);

  if (exb_advection) {
    ddt(Vort) -= scale_ExB * Div_n_bxGrad_f_B_XPPM(Vort, phi, bndry_flux, poloidal_flows) * bracket_factor;
  }

  if (state.isSection("fields") and state["fields"].isSet("DivJextra")) {
    auto DivJextra = get<Field3D>(state["fields"]["DivJextra"]);

    // Parallel current is handled here, to allow different 2D or 3D closures
    // to be used
    ddt(Vort) += DivJextra;
  }

  // Parallel current due to species parallel flow
  for (auto& kv : allspecies.getChildren()) {
    const Options& species = kv.second;

    if (!species.isSet("charge") or !species.isSet("momentum")) {
      continue; // Not charged, or no parallel flow
    }
    const BoutReal Z = get<BoutReal>(species["charge"]);
    if (fabs(Z) < 1e-5) {
      continue; // Not charged
    }

    const Field3D N = get<Field3D>(species["density"]);
    const Field3D NV = get<Field3D>(species["momentum"]);
    const BoutReal A = get<BoutReal>(species["AA"]);

    // Note: Using NV rather than N*V so that the cell boundary flux is correct
    ddt(Vort) += Div_par((Z / A) * NV);
  }

  if (phi_dissipation) {
    // Adds dissipation term like in other equations, but depending on gradient of
    // potential
    Field3D sound_speed = get<Field3D>(state["sound_speed"]);

    Field3D zero {0.0};
    zero.splitParallelSlices();
    zero.yup() = 0.0;
    zero.ydown() = 0.0;

    Field3D dummy;
    ddt(Vort) -= FV::Div_par_mod<hermes::Limiter>(-phi, zero, sound_speed, dummy);
  }

  // Viscosity
  Field3D flow_xlow,flow_zlow;
  ddt(Vort) += (*dagp)(viscosity, Vort, flow_xlow, flow_zlow, false);
  ddt(Vort) += viscosity_core * (*dagp)(ones, Vort, flow_xlow, flow_zlow, false);
  
  Field3D dummy;
  ddt(Vort) += Div_par_K_Grad_par_mod(viscosity_par, Vort, dummy, false);
  
  // Solve diffusion equation for potential

  if (vort_timedissipation > 0.0) {
    ddt(Vort) -= vort_timedissipation * Vort;
  } 
  
  if (boussinesq) {

    Field3D flow_xlow_phi,flow_zlow_phi;
    ddt(phi1) = lambda_1 *( (*dagp)(average_atomic_mass / Bsq, phi, flow_xlow_phi, flow_zlow_phi, false) - Vort);

    if (diamagnetic_polarisation) {
      for (auto& kv : allspecies.getChildren()) {
        // Note: includes electrons (should it?)

        const Options& species = kv.second;
        if (!species.isSet("charge")) {
          continue; // Not charged
        }
        const BoutReal Z = get<BoutReal>(species["charge"]);
        if (fabs(Z) < 1e-5) {
          continue; // Not charged
        }
        if (!species.isSet("pressure")) {
          continue; // No pressure
        }
        const BoutReal A = get<BoutReal>(species["AA"]);
        const Field3D P = get<Field3D>(species["pressure"]);
	Field3D flow_xlow_dia,flow_zlow_dia;
	ddt(phi1) += lambda_1 * (*dagp)(A / Bsq, P, flow_xlow_dia, flow_zlow_dia, false);
      }
    }
  } else {
    // Non-Boussinesq. Calculate mass density by summing over species
    throw BoutException("Non_boussinesq not implemented");
    // Calculate vorticity from potential phi
    Field3D phi_vort = 0.0;
    for (auto& kv : allspecies.getChildren()) {
      const Options& species = kv.second;

      if (!species.isSet("charge")) {
        continue; // Not charged
      }
      const BoutReal Zi = get<BoutReal>(species["charge"]);
      if (fabs(Zi) < 1e-5) {
        continue; // Not charged
      }

      const BoutReal Ai = get<BoutReal>(species["AA"]);
      const Field3D Ni = get<Field3D>(species["density"]);
      phi_vort += Div_a_Grad_perp((Ai / Bsq) * Ni, phi);

      if (diamagnetic_polarisation and species.isSet("pressure")) {
        // Calculate the diamagnetic flow contribution
        const Field3D Pi = get<Field3D>(species["pressure"]);
        phi_vort += Div_a_Grad_perp(Ai / Bsq, Pi);
      }
    }

    ddt(phi1) = lambda_1 * (phi_vort - Vort);
  }

  if (disable_ddt) {
    ddt(Vort) = 0.0;
    ddt(phi1) = 0.0;
  }
  
}

void RelaxPotential::outputVars(Options& state) {
  AUTO_TRACE();
  // Normalisations
  auto Tnorm = get<BoutReal>(state["Tnorm"]);

  set_with_attrs(state["phi"], phi,
                 {{"time_dimension", "t"},
                  {"units", "V"},
                  {"conversion", Tnorm},
                  {"standard_name", "potential"},
                  {"long_name", "plasma potential"},
                  {"source", "relax_potential"}});
}
