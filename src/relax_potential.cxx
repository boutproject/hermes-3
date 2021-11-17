#include <bout/fv_ops.hxx>
#include <bout/solver.hxx>

using bout::globals::mesh;

#include "../include/relax_potential.hxx"
#include "../include/div_ops.hxx"

RelaxPotential::RelaxPotential(std::string name, Options& alloptions, Solver* solver) {
  AUTO_TRACE();

  solver->add(Vort, "Vort"); // Vorticity evolving
  solver->add(phi, "phi1");  // Evolving scaled potential ϕ_1 = λ_2 ϕ
  SAVE_REPEAT(phi);

  // Read curvature vector
  try {
    Curlb_B.covariant = false; // Contravariant
    mesh->get(Curlb_B, "bxcv");

  } catch (BoutException& e) {
    try {
      // May be 2D, reading as 3D
      Vector2D curv2d;
      curv2d.covariant = false;
      mesh->get(curv2d, "bxcv");
      Curlb_B = curv2d;
    } catch (BoutException& e) {
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

  auto* coord = mesh->getCoordinates();

  Curlb_B *= 2. / coord->Bxy;

  Bsq = SQ(coord->Bxy);
}

void RelaxPotential::transform(Options& state) {
  AUTO_TRACE();

  auto& fields = state["fields"];

  // Scale potential
  phi = phi1 / lambda_2;

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
