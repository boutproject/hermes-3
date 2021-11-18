#include <bout/fv_ops.hxx>
#include <bout/solver.hxx>

using bout::globals::mesh;

#include "../include/div_ops.hxx"
#include "../include/relax_potential.hxx"

RelaxPotential::RelaxPotential(std::string name, Options& alloptions, Solver* solver) {
  AUTO_TRACE();

  auto* coord = mesh->getCoordinates();

  auto& options = alloptions[name];

  exb_advection = options["exb_advection"]
                      .doc("Include nonlinear ExB advection?")
                      .withDefault<bool>(true);

  diamagnetic =
      options["diamagnetic"].doc("Include diamagnetic current?").withDefault<bool>(true);

  diamagnetic_polarisation =
      options["diamagnetic_polarisation"]
          .doc("Include diamagnetic drift in polarisation current?")
          .withDefault<bool>(true);

  boussinesq = options["boussinesq"]
                   .doc("Use the Boussinesq approximation?")
                   .withDefault<bool>(true);

  average_atomic_mass = options["average_atomic_mass"]
                            .doc("Weighted average atomic mass, for polarisaion current "
                                 "(Boussinesq approximation)")
                            .withDefault<BoutReal>(2.0); // Deuterium

  poloidal_flows =
      options["poloidal_flows"].doc("Include poloidal ExB flow").withDefault<bool>(true);

  lambda_1 = options["lambda"].doc("λ_1").withDefault(10.0);
  lambda_2 = options["lambda_2"].doc("λ_2").withDefault(10.0);

  solver->add(Vort, "Vort"); // Vorticity evolving
  solver->add(phi, "phi1");  // Evolving scaled potential ϕ_1 = λ_2 ϕ
  SAVE_REPEAT(phi);

  if (diamagnetic) {
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
  }

  Bsq = SQ(coord->Bxy);
}

void RelaxPotential::transform(Options& state) {
  AUTO_TRACE();

  // Scale potential
  phi = phi1 / lambda_2;

  auto& fields = state["fields"];

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

void RelaxPotential::finally(const Options& state) {
  AUTO_TRACE();

  if (exb_advection) {
    ddt(Vort) -= Div_n_bxGrad_f_B_XPPM(Vort, phi, bndry_flux, poloidal_flows);
  }

  if (state.isSection("fields") and state["fields"].isSet("DivJextra")) {
    auto DivJextra = get<Field3D>(state["fields"]["DivJextra"]);

    // Parallel current is handled here, to allow different 2D or 3D closures
    // to be used
    ddt(Vort) += DivJextra;
  }

  // Solve diffusion equation for potential

  if (boussinesq) {
    ddt(phi1) =
        lambda_1 * (FV::Div_a_Laplace_perp(average_atomic_mass / Bsq, phi) - Vort);
  } else {
    // Non-Boussinesq. Calculate mass density by summing over species

    ddt(phi1) =
        lambda_1 * (FV::Div_a_Laplace_perp(average_atomic_mass / Bsq, phi) - Vort);
  }
}
