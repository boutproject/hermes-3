
#include "../include/vorticity.hxx"
#include "../include/div_ops.hxx"

#include <invert_laplace.hxx>
#include <bout/invert/laplacexy.hxx>
#include <bout/constants.hxx>
#include <difops.hxx>

using bout::globals::mesh;

Vorticity::Vorticity(std::string name, Options &alloptions, Solver *solver) {
  AUTO_TRACE();
  
  solver->add(Vort, "Vort");

  SAVE_REPEAT(phi);
  
  auto& options = alloptions[name];

  exb_advection = options["exb_advection"]
                      .doc("Include ExB advection (nonlinear term)?")
                      .withDefault<bool>(true);

  diamagnetic = options["diamagnetic"]
                    .doc("Include diamagnetic current?")
                    .withDefault<bool>(true);

  sheath_boundary =
      options["sheath_boundary"]
          .doc("Set potential to j=0 sheath at boundaries? (default = 0)")
          .withDefault<bool>(false);

  diamagnetic_polarisation =
      options["diamagnetic_polarisation"]
          .doc("Include diamagnetic drift in polarisation current?")
          .withDefault<bool>(true);

  average_atomic_mass =
      options["average_atomic_mass"]
          .doc("Weighted average atomic mass, for polarisaion current "
               "(Boussinesq approximation)")
          .withDefault<BoutReal>(2.0); // Deuterium

  bndry_flux = options["bndry_flux"]
                   .doc("Allow flows through radial boundaries")
                   .withDefault<bool>(true);

  poloidal_flows = options["poloidal_flows"]
                       .doc("Include poloidal ExB flow")
                       .withDefault<bool>(true);

  split_n0 = options["split_n0"]
                 .doc("Split phi into n=0 and n!=0 components")
                 .withDefault<bool>(false);

  auto coord = mesh->getCoordinates();
  
  if (split_n0) {
    // Create an XY solver for n=0 component
    laplacexy = new LaplaceXY(mesh);
    // Set coefficients for Boussinesq solve
    laplacexy->setCoefs(1. / SQ(coord->Bxy), 0.0);
    phi2D = 0.0; // Starting guess
  }
  phiSolver = Laplacian::create(&options["laplacian"]);
  // Set coefficients for Boussinesq solve
  phiSolver->setCoefC(1. / SQ(coord->Bxy));
  phi = 0.0;
  
  // Read curvature vector
  try {
    Curlb_B.covariant = false; // Contravariant
    mesh->get(Curlb_B, "bxcv");
  
  } catch (BoutException &e) {
    try {
      // May be 2D, reading as 3D
      Vector2D curv2d;
      curv2d.covariant = false;
      mesh->get(curv2d, "bxcv");
      Curlb_B = curv2d;
    } catch (BoutException &e) {
      if (diamagnetic) {
        // Need curvature
        throw;
      } else {
        output_warn.write("No curvature vector in input grid");
        Curlb_B = 0.0;
      }
    }
  }

  if (Options::root()["mesh"]["paralleltransform"].as<std::string>() == "shifted") {
    Field2D I;
    mesh->get(I, "sinty");
    Curlb_B.z += I * Curlb_B.x;
  }
  
  Options& units = alloptions["units"];
  BoutReal Bnorm = units["Tesla"];
  BoutReal Lnorm = units["meters"];
  
  Curlb_B.x /= Bnorm;
  Curlb_B.y *= SQ(Lnorm);
  Curlb_B.z *= SQ(Lnorm);

  Curlb_B *= 2. / coord->Bxy;
  
  Bsq = SQ(coord->Bxy);
}

void Vorticity::transform(Options &state) {
  AUTO_TRACE();
  
  auto& fields = state["fields"];

  // Sheath multiplier Te -> phi (2.84522 for Deuterium)
  BoutReal sheathmult = 0.0;
  if (sheath_boundary) {
    sheathmult = log(0.5 * sqrt(SI::Mp / SI::Me / PI));
  }
  
  Field3D Te; // Electron temperature, use for outer boundary conditions
  if (state["species"]["e"].isSet("temperature")) {
    // Electron temperature set
    Te = get<Field3D>(state["species"]["e"]["temperature"]);
  } else {
    Te = 0.0;
  }
  
  // Calculate potential 
  if (split_n0) {
    ////////////////////////////////////////////
    // Split into axisymmetric and non-axisymmetric components
    Field2D Vort2D = DC(Vort); // n=0 component
    
    // Set the boundary to 2.8*Te
    phi2D.setBoundaryTo(sheathmult * DC(Te));
    
    phi2D = laplacexy->solve(Vort2D, phi2D);
    
    // Solve non-axisymmetric part using X-Z solver
    phi = phi2D + phiSolver->solve((Vort - Vort2D) * Bsq,
                                   sheathmult * (Te - DC(Te)));
    
  } else {
    phi = phiSolver->solve(Vort * Bsq, (sheathmult * Te));
  }

  if (diamagnetic_polarisation) {
    // Diamagnetic term in vorticity. Note this is weighted by the mass/charge ratio
    // This includes all species, including electrons
    Options& allspecies = state["species"];
    for (auto& kv : allspecies.getChildren()) {
      Options& species = allspecies[kv.first]; // Note: need non-const
      
      if (!(species.isSet("pressure") and species.isSet("charge") and species.isSet("AA"))) {
        continue; // No pressure, charge or mass -> no polarisation current
      }
      
      auto P = get<Field3D>(species["pressure"]);
      auto AA = get<BoutReal>(species["AA"]);
      auto charge = get<BoutReal>(species["charge"]);
      
      phi -= P * (AA / charge);
    }
  }

  // Account for mix of ions
  phi  /= average_atomic_mass;
  
  ddt(Vort) = 0.0;
  
  if (diamagnetic) {
    // Diamagnetic current. This is calculated here so that the energy sources/sinks
    // can be calculated for the evolving species.

    Vector3D Jdia; Jdia.x = 0.0; Jdia.y = 0.0; Jdia.z = 0.0;
    Jdia.covariant = Curlb_B.covariant;

    Options& allspecies = state["species"];
    
    for (auto& kv : allspecies.getChildren()) {
      Options& species = allspecies[kv.first]; // Note: need non-const

      if (!(species.isSet("pressure") and species.isSet("charge"))) {
        continue; // No pressure or charge -> no diamagnetic current
      }
      // Note that the species must have a charge, but charge is not used,
      // because it cancels out in the expression for current
      
      auto P = get<Field3D>(species["pressure"]);

      Vector3D Jdia_species = P * Curlb_B; // Diamagnetic current for this species
      
      // This term energetically balances diamagnetic term
      // in the vorticity equation
      subtract(species["energy_source"],
               Jdia_species * Grad(phi));

      Jdia += Jdia_species; // Collect total diamagnetic current
    }
    
    // Note: This term is central differencing so that it balances
    // the corresponding compression term in the species pressure equations
    Field3D DivJdia = Div(Jdia);
    ddt(Vort) += DivJdia;

    if (diamagnetic_polarisation) {
      // Calculate energy exchange term nonlinear in pressure
      // ddt(Pi) += Pi * Div((Pe + Pi) * Curlb_B);
      for (auto& kv : allspecies.getChildren()) {
        Options& species = allspecies[kv.first]; // Note: need non-const

        if (!(species.isSet("pressure") and species.isSet("charge") and species.isSet("AA"))) {
          continue; // No pressure, charge or mass -> no polarisation current due to diamagnetic flow
        }
        auto P = get<Field3D>(species["pressure"]);
        auto AA = get<BoutReal>(species["AA"]);
        auto charge = get<BoutReal>(species["charge"]);
        
        add(species["energy_source"],
            (3./2) * P * (AA / average_atomic_mass / charge) * DivJdia);
      }
    }

    set(fields["DivJdia"], DivJdia);
  }

  set(fields["vorticity"], Vort);
  set(fields["phi"], phi);
}
  
void Vorticity::finally(const Options &state) {
  AUTO_TRACE();
  
  if (exb_advection) {
    ddt(Vort) -= Div_n_bxGrad_f_B_XPPM(Vort, phi, bndry_flux, poloidal_flows);

    /*
      ddt(Vort) -=
        Div_n_bxGrad_f_B_XPPM(0.5 * Vort, phi, bndry_flux, poloidal_flows);

    // V_ExB dot Grad(Pi)
    Field3D vEdotGradPi = bracket(phi, Pi, BRACKET_ARAKAWA);
    vEdotGradPi.applyBoundary("free_o2");
    // delp2(phi) term
    Field3D DelpPhi_2B2 = 0.5 * Delp2(phi) / Bsq;
    DelpPhi_2B2.applyBoundary("free_o2");
    
    mesh->communicate(vEdotGradPi, DelpPhi_2B2);
    
    ddt(Vort) -= FV::Div_a_Laplace_perp(0.5 / Bsq, vEdotGradPi);
    */
  }
  
  if (state.isSection("fields") and state["fields"].isSet("DivJextra")) {
    auto DivJextra = get<Field3D>(state["fields"]["DivJextra"]);

    // Parallel current is handled here, to allow different 2D or 3D closures
    // to be used
    ddt(Vort) += DivJextra;
    
    // This term is central differencing so that it balances the parallel gradient
    // of the potential in Ohm's law
    //ddt(Vort) += Div_par(Jpar);
  }
  
}


