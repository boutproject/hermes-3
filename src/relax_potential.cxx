
#include "../include/relax_potential.hxx"
#include "../include/component.hxx"
#include "../include/div_ops.hxx"
#include "../include/guarded_options.hxx"
#include "../include/hermes_utils.hxx"
#include "../include/permissions.hxx"

#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/constants.hxx>
#include <bout/derivs.hxx>
#include <bout/difops.hxx>
#include <bout/field.hxx>
#include <bout/field2d.hxx>
#include <bout/field3d.hxx>
#include <bout/fv_ops.hxx>
#include <bout/globals.hxx>
#include <bout/options.hxx>
#include <bout/output.hxx>
#include <bout/solver.hxx>
#include <bout/sys/range.hxx>
#include <bout/utils.hxx>
#include <bout/vecops.hxx>
#include <bout/vector3d.hxx>

#include <cmath>
#include <string>

using bout::globals::mesh;

namespace {
/// Limited free gradient of log of a quantity
/// This ensures that the guard cell values remain positive
/// while also ensuring that the quantity never increases
///
///  fm  fc | fp
///         ^ boundary
///
/// exp( 2*log(fc) - log(fm) )
///
BoutReal limitFree(BoutReal fm, BoutReal fc) {
  if (fm < fc) {
    return fc; // Neumann rather than increasing into boundary
  }
  if (fm < 1e-10) {
    return fc; // Low / no density condition
  }
  BoutReal fp = SQ(fc) / fm;
#if CHECKLEVEL >= 2
  if (!std::isfinite(fp)) {
    throw BoutException("SheathBoundary limitFree: {}, {} -> {}", fm, fc, fp);
  }
#endif

  return fp;
}

void applyParallelNeumannBoundary(Field3D& f) {
  if (f.hasParallelSlices()) {
    Field3D& f_ydown = f.ydown();
    Field3D& f_yup = f.yup();
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        f_ydown(r.ind, mesh->ystart - 1, jz) = f(r.ind, mesh->ystart, jz);
      }
    }
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        f_yup(r.ind, mesh->yend + 1, jz) = f(r.ind, mesh->yend, jz);
      }
    }
    return;
  }

  Field3D f_fa = toFieldAligned(f);
  for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
    for (int jz = 0; jz < mesh->LocalNz; jz++) {
      f_fa(r.ind, mesh->ystart - 1, jz) = f_fa(r.ind, mesh->ystart, jz);
    }
  }
  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
    for (int jz = 0; jz < mesh->LocalNz; jz++) {
      f_fa(r.ind, mesh->yend + 1, jz) = f_fa(r.ind, mesh->yend, jz);
    }
  }
  f = fromFieldAligned(f_fa);
}
} // namespace

RelaxPotential::RelaxPotential(std::string name, Options& alloptions, Solver* solver)
    : Component({
          readWrite("fields:vorticity"),
          readWrite("fields:phi"),
          readIfSet("species:{all_species}:charge"),
          readOnly("species:{charged}:AA"),
      }) {

  auto& options = alloptions[name];

  evolve_vorticity = options["evolve_vorticity"].doc("").withDefault<bool>(true);

  if (evolve_vorticity) {
    solver->add(Vort, "Vort"); // Vorticity evolving
  }
  solver->add(phi1, "phi1"); // Evolving scaled potential ϕ_1 = λ_2 ϕ

  // Normalisations
  const Options& units = alloptions["units"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();
  const BoutReal Bnorm = units["Tesla"];
  const BoutReal Lnorm = units["meters"];
  const BoutReal Tnorm = units["eV"];
  const BoutReal Nnorm = units["inv_meters_cubed"];

  exb_advection = options["exb_advection"]
                      .doc("Include nonlinear ExB advection?")
                      .withDefault<bool>(true);

  exb_advection_simplified = options["exb_advection_simplified"]
                                 .doc("Simplify nonlinear ExB advection form?")
                                 .withDefault<bool>(true);

  diamagnetic =
      options["diamagnetic"].doc("Include diamagnetic current?").withDefault<bool>(true);

  sheath_boundary =
      options["sheath_boundary"]
          .doc("Set potential to j=0 sheath at radial boundaries? (default = 0)")
          .withDefault<bool>(false);

  diamagnetic_polarisation =
      options["diamagnetic_polarisation"]
          .doc("Include diamagnetic drift in polarisation current?")
          .withDefault<bool>(true);

  collisional_friction =
      options["collisional_friction"]
          .doc("Damp vorticity based on mass-weighted collision frequency")
          .withDefault<bool>(false);

  boussinesq = options["boussinesq"]
                   .doc("Use the Boussinesq approximation?")
                   .withDefault<bool>(true);

  average_atomic_mass = options["average_atomic_mass"]
                            .doc("Weighted average atomic mass, for polarisation current "
                                 "(Boussinesq approximation)")
                            .withDefault<BoutReal>(2.0); // Deuterium

  bndry_flux = options["bndry_flux"]
                   .doc("Allow flows through radial boundaries")
                   .withDefault<bool>(true);

  poloidal_flows =
      options["poloidal_flows"].doc("Include poloidal ExB flow").withDefault<bool>(true);

  viscosity = 0.0;
  viscosity = options["viscosity"]
                  .doc("Perpendicular Kinematic viscosity [m^2/s]")
                  .withDefault(viscosity)
              / (Lnorm * Lnorm * Omega_ci);
  viscosity.applyBoundary("dirichlet");

  viscosity_par = 0.0;
  viscosity_par = options["viscosity_par"]
                      .doc("Parallel Kinematic viscosity [m^2/s]")
                      .withDefault(viscosity_par)
                  / (Lnorm * Lnorm * Omega_ci);
  viscosity_par.applyBoundary("dirichlet");
  viscosity_par.applyParallelBoundary("parallel_dirichlet_o2");

  hyper_z = options["hyper_z"].doc("Hyper-viscosity in Z. < 0 -> off").withDefault(-1.0);

  // Numerical dissipation terms
  // These are required to suppress parallel zig-zags in
  // cell centred formulations. Essentially adds (hopefully small)
  // parallel currents

  vort_dissipation = options["vort_dissipation"]
                         .doc("Parallel dissipation of vorticity")
                         .withDefault<bool>(false);

  phi_dissipation = options["phi_dissipation"]
                        .doc("Parallel dissipation of potential [Recommended]")
                        .withDefault<bool>(true);

  phi_boundary_relax = options["phi_boundary_relax"]
                           .doc("Relax x boundaries of phi towards Neumann?")
                           .withDefault<bool>(false);

  phi_sheath_dissipation = options["phi_sheath_dissipation"]
                               .doc("Add dissipation when phi < 0.0 at the sheath")
                               .withDefault<bool>(false);

  damp_core_vorticity = options["damp_core_vorticity"]
                            .doc("Damp vorticity at the core boundary?")
                            .withDefault<bool>(false);

  lambda_1 = options["lambda_1"].doc("λ_1 > 1").withDefault(3e8)
             / (Tnorm * Omega_ci / SI::qe / Nnorm);
  lambda_2 = options["lambda_2"].doc("λ_2 > λ_1").withDefault(1.0);

  // Add phi to restart files so that the value in the boundaries
  // is restored on restart. This is done even when phi is not evolving,
  // so that phi can be saved and re-loaded

  // Set initial value. Will be overwritten if restarting
  phi1 = 0.0;
  Vort = 0.0;

  const auto* coord = mesh->getCoordinates();

  if (phi_boundary_relax) {
    // Set the last update time to -1, so it will reset
    // the first time RHS function is called
    phi_boundary_last_update = -1.;

    phi_core_averagey = options["phi_core_averagey"]
                            .doc("Average phi core boundary in Y?")
                            .withDefault<bool>(false)
                        and mesh->periodicY(mesh->xstart);

    phi_boundary_timescale = options["phi_boundary_timescale"]
                                 .doc("Timescale for phi boundary relaxation [seconds]")
                                 .withDefault(1e-4)
                             / get<BoutReal>(alloptions["units"]["seconds"]);
  } else {
    setPermissions(readIfSet("species:{charged}:temperature", Regions::Interior));
  }

  if (diamagnetic) {
    // FIXME: These will only be read if BOTH charge and pressure are set
    setPermissions(readIfSet("species:{charged}:pressure", Regions::Interior));
    setPermissions(readIfSet("species:{all_species}:charge"));
    // FIXME: The weay transform_impl is currently written,
    // energy_source is set for neutral species with an explicit
    // charge declared as 0 if diamagnetic_polarisation == true. I
    // suspect that's a mistake though.
    setPermissions(readWrite("species:{all_species}:energy_source"));
    setPermissions(readWrite("fields:DivJdia"));

    // Read curvature vector
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
      }
      output_warn.write("No curvature vector in input grid");
      Curlb_B = 0.0;
    }
  }

  if (collisional_friction) {
    setPermissions(readOnly("species:{charged}:density", Regions::Interior));
    setPermissions(readIfSet("species:{positive_ions}:collision_frequency"));
    setPermissions(readWrite("fields:DivJcol"));
  }

  if (Options::root()["mesh"]["paralleltransform"]["type"].as<std::string>()
      == "shifted") {
    Field2D I;
    mesh->get(I, "sinty");
    Curlb_B.z += I * Curlb_B.x;
  }

  Curlb_B.x /= Bnorm;
  Curlb_B.y *= SQ(Lnorm);
  Curlb_B.z *= SQ(Lnorm);

  Curlb_B *= 2. / coord->Bxy;

  Bsq = SQ(coord->Bxy);

  diagnose =
      options["diagnose"].doc("Output additional diagnostics?").withDefault<bool>(false);
}

Field3D RelaxPotential::calculatePihat(GuardedOptions allspecies) {
  Pi_hat = 0.0;
  if (!diamagnetic_polarisation) {
    Pi_hat.applyBoundary("neumann");
    return Pi_hat;
  }

  for (auto& kv : allspecies.getChildren()) {
    const GuardedOptions species = allspecies[kv.first];

    if (!(IS_SET_NOBOUNDARY(species["pressure"]) and species.isSet("charge")
          and species.isSet("AA"))) {
      continue;
    }

    const auto charge = get<BoutReal>(species["charge"]);
    if (fabs(charge) < 1e-5) {
      continue;
    }

    const auto P = GET_NOBOUNDARY(Field3D, species["pressure"]);
    const auto AA = get<BoutReal>(species["AA"]);
    Pi_hat += P * (AA / average_atomic_mass / charge);
  }

  Pi_hat.applyBoundary("neumann");
  return Pi_hat;
}

void RelaxPotential::applyPhiBoundary(Field3D& phi, GuardedOptions state) {
  if (phi_boundary_relax) {
    const BoutReal time = get<BoutReal>(state["time"]);

    if (phi_boundary_last_update < 0.0) {
      phi_boundary_last_update = time;
      return;
    }
    if (time <= phi_boundary_last_update) {
      return;
    }

    const BoutReal weight =
        exp(-(time - phi_boundary_last_update) / phi_boundary_timescale);
    phi_boundary_last_update = time;

    if (mesh->firstX()) {
      BoutReal phivalue = 0.0;
      if (phi_core_averagey) {
        BoutReal philocal = 0.0;
        for (int j = mesh->ystart; j <= mesh->yend; j++) {
          for (int k = 0; k < mesh->LocalNz; k++) {
            philocal += phi(mesh->xstart, j, k);
          }
        }
        MPI_Comm comm_inner = mesh->getYcomm(0);
        int np;
        MPI_Comm_size(comm_inner, &np);
        MPI_Allreduce(&philocal, &phivalue, 1, MPI_DOUBLE, MPI_SUM, comm_inner);
        phivalue /= (np * mesh->LocalNz * mesh->LocalNy);
      }

      for (int j = mesh->ystart; j <= mesh->yend; j++) {
        if (!phi_core_averagey) {
          phivalue = 0.0;
          for (int k = 0; k < mesh->LocalNz; k++) {
            phivalue += phi(mesh->xstart, j, k);
          }
          phivalue /= mesh->LocalNz;
        }

        const BoutReal oldvalue =
            0.5 * (phi(mesh->xstart - 1, j, 0) + phi(mesh->xstart, j, 0));
        const BoutReal newvalue = (weight * oldvalue) + ((1. - weight) * phivalue);

        for (int k = 0; k < mesh->LocalNz; k++) {
          phi(mesh->xstart - 1, j, k) = 2. * newvalue - phi(mesh->xstart, j, k);
        }
        if (mesh->xstart > 1) {
          for (int k = 0; k < mesh->LocalNz; k++) {
            phi(mesh->xstart - 2, j, k) = phi(mesh->xstart - 1, j, k);
          }
        }
      }
    }

    if (mesh->lastX()) {
      for (int j = mesh->ystart; j <= mesh->yend; j++) {
        BoutReal phivalue = 0.0;
        for (int k = 0; k < mesh->LocalNz; k++) {
          phivalue += phi(mesh->xend, j, k);
        }
        phivalue /= mesh->LocalNz;

        const BoutReal oldvalue =
            0.5 * (phi(mesh->xend + 1, j, 0) + phi(mesh->xend, j, 0));
        const BoutReal newvalue = (weight * oldvalue) + ((1. - weight) * phivalue);

        for (int k = 0; k < mesh->LocalNz; k++) {
          phi(mesh->xend + 1, j, k) = 2. * newvalue - phi(mesh->xend, j, k);
        }
        if (mesh->LocalNx - mesh->xend > 2) {
          for (int k = 0; k < mesh->LocalNz; k++) {
            phi(mesh->xend + 2, j, k) = phi(mesh->xend + 1, j, k);
          }
        }
      }
    }
    return;
  }

  BoutReal sheathmult = 0.0;
  if (sheath_boundary) {
    const BoutReal Me_Mp = get<BoutReal>(state["species"]["e"]["AA"]);
    sheathmult = log(0.5 * sqrt(1. / (Me_Mp * PI)));
  }

  Field3D Te;
  if (state["species"]["e"].isSet("temperature")) {
    Te = GET_NOBOUNDARY(Field3D, state["species"]["e"]["temperature"]);
  } else {
    Te = 0.0;
  }

  if (mesh->firstX()) {
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      BoutReal teavg = 0.0;
      for (int k = 0; k < mesh->LocalNz; k++) {
        teavg += Te(mesh->xstart, j, k);
      }
      teavg /= mesh->LocalNz;
      const BoutReal phivalue = sheathmult * teavg;
      for (int k = 0; k < mesh->LocalNz; k++) {
        phi(mesh->xstart - 1, j, k) = 2. * phivalue - phi(mesh->xstart, j, k);
      }
      if (mesh->xstart > 1) {
        for (int k = 0; k < mesh->LocalNz; k++) {
          phi(mesh->xstart - 2, j, k) = phi(mesh->xstart - 1, j, k);
        }
      }
    }
  }

  if (mesh->lastX()) {
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      BoutReal teavg = 0.0;
      for (int k = 0; k < mesh->LocalNz; k++) {
        teavg += Te(mesh->xend, j, k);
      }
      teavg /= mesh->LocalNz;
      const BoutReal phivalue = sheathmult * teavg;
      for (int k = 0; k < mesh->LocalNz; k++) {
        phi(mesh->xend + 1, j, k) = 2. * phivalue - phi(mesh->xend, j, k);
      }
      if (mesh->xend < mesh->LocalNx - 2) {
        for (int k = 0; k < mesh->LocalNz; k++) {
          phi(mesh->xend + 2, j, k) = phi(mesh->xend + 1, j, k);
        }
      }
    }
  }
}

Field3D RelaxPotential::calculateDivJdia(Field3D& phi, GuardedOptions allspecies) {
  Vector3D Jdia;
  Jdia.x = 0.0;
  Jdia.y = 0.0;
  Jdia.z = 0.0;
  Jdia.covariant = Curlb_B.covariant;

  if (phi.hasParallelSlices()) {
    Field3D& phi_ydown = phi.ydown();
    Field3D& phi_yup = phi.yup();
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        phi_ydown(r.ind, mesh->ystart - 1, jz) =
            2 * phi(r.ind, mesh->ystart, jz) - phi_yup(r.ind, mesh->ystart + 1, jz);
      }
    }
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        phi_yup(r.ind, mesh->yend + 1, jz) =
            2 * phi(r.ind, mesh->yend, jz) - phi_ydown(r.ind, mesh->yend - 1, jz);
      }
    }
  } else {
    Field3D phi_fa = toFieldAligned(phi);
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        phi_fa(r.ind, mesh->ystart - 1, jz) =
            2 * phi_fa(r.ind, mesh->ystart, jz) - phi_fa(r.ind, mesh->ystart + 1, jz);
      }
    }
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        phi_fa(r.ind, mesh->yend + 1, jz) =
            2 * phi_fa(r.ind, mesh->yend, jz) - phi_fa(r.ind, mesh->yend - 1, jz);
      }
    }
    phi = fromFieldAligned(phi_fa);
  }

  const Vector3D Grad_phi = Grad(phi);

  for (auto& kv : allspecies.getChildren()) {
    const GuardedOptions species = allspecies[kv.first];

    if (!(IS_SET_NOBOUNDARY(species["pressure"]) and IS_SET(species["charge"]))) {
      continue;
    }
    if (fabs(get<BoutReal>(species["charge"])) < 1e-5) {
      continue;
    }

    auto P = GET_NOBOUNDARY(Field3D, species["pressure"]);

    if (P.hasParallelSlices()) {
      Field3D& P_ydown = P.ydown();
      Field3D& P_yup = P.yup();
      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          P_ydown(r.ind, mesh->ystart - 1, jz) =
              2 * P(r.ind, mesh->ystart, jz) - P_yup(r.ind, mesh->ystart + 1, jz);
        }
      }
      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          P_yup(r.ind, mesh->yend + 1, jz) =
              2 * P(r.ind, mesh->yend, jz) - P_ydown(r.ind, mesh->yend - 1, jz);
        }
      }
    } else {
      Field3D P_fa = toFieldAligned(P);
      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          auto i = indexAt(P_fa, r.ind, mesh->ystart, jz);
          P_fa[i.ym()] = limitFree(P_fa[i.yp()], P_fa[i]);
        }
      }
      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          auto i = indexAt(P_fa, r.ind, mesh->yend, jz);
          P_fa[i.yp()] = limitFree(P_fa[i.ym()], P_fa[i]);
        }
      }
      P = fromFieldAligned(P_fa);
    }

    const Vector3D Jdia_species = P * Curlb_B;
    subtract(species["energy_source"], Jdia_species * Grad_phi);
    Jdia += Jdia_species;
  }

  return Div(Jdia);
}

Field3D RelaxPotential::calculateDivJcol(const Field3D& phi, const Field3D& pi_hat,
                                         GuardedOptions allspecies) {
  Field3D sum_A_nu_n = zeroFrom(phi);
  Field3D sum_A_n = zeroFrom(phi);

  for (const auto& kv : allspecies.getChildren()) {
    const auto species = kv.second;

    if (!(species.isSet("charge") and species.isSet("AA"))) {
      continue;
    }
    if (fabs(get<BoutReal>(species["charge"])) < 1e-5) {
      continue;
    }

    const BoutReal A = get<BoutReal>(species["AA"]);
    const Field3D N = GET_NOBOUNDARY(Field3D, species["density"]);
    const Field3D AN = A * N;
    sum_A_n += AN;
    if (IS_SET(species["collision_frequency"])) {
      sum_A_nu_n += AN * GET_VALUE(Field3D, species["collision_frequency"]);
    }
  }

  Field3D weighted_collision_frequency = sum_A_nu_n / sum_A_n;
  weighted_collision_frequency.applyBoundary("neumann");

  return -FV::Div_a_Grad_perp(weighted_collision_frequency * average_atomic_mass / Bsq,
                              phi + pi_hat);
}

Field3D RelaxPotential::calculateExBAdvectionSource(const Field3D& vort,
                                                    const Field3D& phi,
                                                    const Field3D& pi_hat) {
  Field3D result = zeroFrom(vort);

  if (!exb_advection) {
    return result;
  }

  if (exb_advection_simplified) {
    result -= Div_n_bxGrad_f_B_XPPM(vort, phi, bndry_flux, poloidal_flows);
    return result;
  }

  result -= Div_n_bxGrad_f_B_XPPM(0.5 * vort, phi, bndry_flux, poloidal_flows);

  Field3D vEdotGradPi = bracket(phi, pi_hat, BRACKET_ARAKAWA);
  vEdotGradPi.applyBoundary("free_o2");

  Field3D DelpPhi_2B2 = 0.5 * average_atomic_mass * Delp2(phi) / Bsq;
  DelpPhi_2B2.applyBoundary("free_o2");

  mesh->communicate(vEdotGradPi, DelpPhi_2B2);

  result -= FV::Div_a_Grad_perp(0.5 * average_atomic_mass / Bsq, vEdotGradPi);
  result -= Div_n_bxGrad_f_B_XPPM(DelpPhi_2B2, phi + pi_hat, bndry_flux, poloidal_flows);

  return result;
}

Field3D RelaxPotential::calculateParallelCurrentSource(const Options& state) {
  const auto* coord = mesh->getCoordinates();

  Field3D result{0.0};
  bool initialised = false;
  auto ensure_initialised = [&](const Field3D& like) {
    if (!initialised) {
      result = zeroFrom(like);
      initialised = true;
    }
  };

  if (state["fields"].isSet("DivJextra")) {
    const auto DivJextra = get<Field3D>(state["fields"]["DivJextra"]);
    ensure_initialised(DivJextra);
    result += DivJextra;
  }

  if (!state.isSet("species")) {
    return result;
  }

  for (const auto& kv : state["species"].getChildren()) {
    const Options& species = kv.second;

    if (!species.isSet("charge") or !species.isSet("momentum")) {
      continue;
    }
    const BoutReal Z = get<BoutReal>(species["charge"]);
    if (fabs(Z) < 1e-5) {
      continue;
    }

    const Field3D NV = get<Field3D>(species["momentum"]);
    ensure_initialised(NV);
    const BoutReal A = get<BoutReal>(species["AA"]);
    const Field3D jpar = (Z / A) * NV;
    result += Div_par(jpar);

    if (state["fields"].isSet("Apar_flutter")) {
      const Field3D Apar_flutter = get<Field3D>(state["fields"]["Apar_flutter"]);
      result += coord->Bxy * bracket(jpar / coord->Bxy, Apar_flutter, BRACKET_ARAKAWA);
    }
  }

  return result;
}

Field3D RelaxPotential::calculateDissipationSource(const Options& state,
                                                   const Field3D& vort,
                                                   const Field3D& phi) {
  Field3D result = zeroFrom(vort);

  result += FV::Div_a_Grad_perp(viscosity, vort);
  result += FV::Div_par_K_Grad_par(viscosity_par, vort);

  if (vort_dissipation) {
    const Field3D sound_speed = get<Field3D>(state["sound_speed"]);
    result -= FV::Div_par(vort, 0.0, sound_speed);
  }

  if (phi_dissipation) {
    const Field3D sound_speed = get<Field3D>(state["sound_speed"]);
    result -= FV::Div_par(-phi, 0.0, sound_speed);
  }

  if (hyper_z > 0) {
    auto* coord = vort.getCoordinates();
    result -= hyper_z * SQ(SQ(coord->dz)) * D4DZ4(vort);
  }

  if (phi_sheath_dissipation) {
    auto phi_fa = toFieldAligned(phi);
    Field3D dissipation{zeroFrom(phi_fa)};
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        auto i = indexAt(phi_fa, r.ind, mesh->ystart, jz);
        const BoutReal phisheath = 0.5 * (phi_fa[i] + phi_fa[i.ym()]);
        dissipation[i] = -floor(-phisheath, 0.0);
      }
    }

    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        auto i = indexAt(phi_fa, r.ind, mesh->yend, jz);
        const BoutReal phisheath = 0.5 * (phi_fa[i] + phi_fa[i.yp()]);
        dissipation[i] = -floor(-phisheath, 0.0);
      }
    }
    result += fromFieldAligned(dissipation);
  }

  if (damp_core_vorticity and mesh->firstX() and mesh->periodicY(mesh->xstart)) {
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      BoutReal vort_avg = 0.0;
      for (int k = 0; k < mesh->LocalNz; k++) {
        vort_avg += vort(mesh->xstart, j, k);
      }
      vort_avg /= mesh->LocalNz;
      for (int k = 0; k < mesh->LocalNz; k++) {
        result(mesh->xstart, j, k) -= 0.01 * vort_avg;
      }
    }
  }

  return result;
}

Field3D RelaxPotential::calculatePhi1Source(const Field3D& vort, const Field3D& vort_rhs,
                                            const Field3D& vort_from_phi) const {
  if (evolve_vorticity) {
    return lambda_1 * (vort_from_phi - vort);
  }
  return -lambda_1 * vort_rhs;
}

void RelaxPotential::transform_impl(GuardedOptions& state) {

  auto allspecies = state["species"];

  phi.name = "phi";
  auto fields = state["fields"];

  // Set the boundary of phi1.
  phi1.applyBoundary("neumann");

  // Scale potential
  phi = phi1 / lambda_2;

  // Calculate vorticity from phi.
  this->Vort_from_phi = vorticity(phi, allspecies);
  if (!evolve_vorticity) {
    // Not evolving vorticity so calculate from phi
    Vort = this->Vort_from_phi;
  }

  // Set the boundary of Vort. Needed only if dissipation terms are included.
  Vort.applyBoundary("neumann");
  applyParallelNeumannBoundary(Vort);

  Pi_hat = calculatePihat(allspecies);
  applyPhiBoundary(phi, state);

  // Ensure that potential is set in the communication guard cells
  mesh->communicate(phi, phi1, Vort); // Should we communicate phi1?

  // Vorticity equation

  ddt(Vort) = 0.0;

  if (diamagnetic) {
    DivJdia = calculateDivJdia(phi, allspecies);
    ddt(Vort) += DivJdia;
    set(fields["DivJdia"], DivJdia);
  }

  if (collisional_friction) {
    DivJcol = calculateDivJcol(phi, Pi_hat, allspecies);
    ddt(Vort) += DivJcol;
    set(fields["DivJcol"], DivJcol);
  }

  set(fields["vorticity"], Vort);
  set(fields["phi"], phi);
}

void RelaxPotential::finally(const Options& state) {
  phi = get<Field3D>(state["fields"]["phi"]);

  ddt(Vort) += calculateExBAdvectionSource(Vort, phi, Pi_hat);
  ddt(Vort) += calculateParallelCurrentSource(state);
  ddt(Vort) += calculateDissipationSource(state, Vort, phi);
  ddt(phi1) = calculatePhi1Source(Vort, ddt(Vort), Vort_from_phi);
}

void RelaxPotential::outputVars(Options& state) {
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);

  state["Vort"].setAttributes({{"time_dimension", "t"},
                               {"units", "C m^-3"},
                               {"conversion", SI::qe * Nnorm},
                               {"long_name", "vorticity"},
                               {"source", "relax_potential"}});

  set_with_attrs(state["phi"], phi,
                 {{"time_dimension", "t"},
                  {"units", "V"},
                  {"conversion", Tnorm},
                  {"standard_name", "potential"},
                  {"long_name", "plasma potential"},
                  {"source", "relax_potential"}});

  if (diagnose) {
    set_with_attrs(state["ddt(Vort)"], ddt(Vort),
                   {{"time_dimension", "t"},
                    {"units", "A m^-3"},
                    {"conversion", SI::qe * Nnorm * Omega_ci},
                    {"long_name", "Rate of change of vorticity"},
                    {"source", "relax_potential"}});
    set_with_attrs(state["ddt(phi)"], ddt(phi1),
                   {{"time_dimension", "t"},
                    {"units", "V/s"},
                    {"conversion", Tnorm * Omega_ci},
                    {"long_name", "Rate of change of electrostatic potential"},
                    {"source", "relax_potential"}});
    if (diamagnetic) {
      set_with_attrs(state["DivJdia"], DivJdia,
                     {{"time_dimension", "t"},
                      {"units", "A m^-3"},
                      {"conversion", SI::qe * Nnorm * Omega_ci},
                      {"long_name", "Divergence of diamagnetic current"},
                      {"source", "relax_potential"}});
    }
    if (collisional_friction) {
      set_with_attrs(state["DivJcol"], DivJcol,
                     {{"time_dimension", "t"},
                      {"units", "A m^-3"},
                      {"conversion", SI::qe * Nnorm * Omega_ci},
                      {"long_name", "Divergence of collisional current"},
                      {"source", "relax_potential"}});
    }
  }
}

Field3D RelaxPotential::vorticity(const Field3D& phi, GuardedOptions& allspecies) {
  if (boussinesq) {
    Field3D phi_vort = FV::Div_a_Grad_perp(average_atomic_mass / Bsq, phi);

    if (diamagnetic_polarisation) {
      for (const auto& kv : allspecies.getChildren()) {
        // Note: includes electrons (should it?)

        const GuardedOptions& species = kv.second;
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
        const Field3D P = GET_NOBOUNDARY(Field3D, species["pressure"]);
        phi_vort += FV::Div_a_Grad_perp(A / Bsq, P);
      }
    }
    return phi_vort;
  }
  // Non-Boussinesq. Calculate mass density by summing over species

  // Calculate vorticity from potential phi
  Field3D phi_vort = 0.0;
  for (const auto& kv : allspecies.getChildren()) {
    const GuardedOptions& species = kv.second;

    if (!species.isSet("charge")) {
      continue; // Not charged
    }
    const BoutReal Zi = get<BoutReal>(species["charge"]);
    if (fabs(Zi) < 1e-5) {
      continue; // Not charged
    }

    const BoutReal Ai = get<BoutReal>(species["AA"]);
    const Field3D Ni = get<Field3D>(species["density"]);

    phi_vort += FV::Div_a_Grad_perp((Ai / Bsq) * Ni, phi);

    if (diamagnetic_polarisation and species.isSet("pressure")) {
      // Calculate the diamagnetic flow contribution
      const Field3D Pi = get<Field3D>(species["pressure"]);
      phi_vort += FV::Div_a_Grad_perp(Ai / Bsq / Zi, Pi);
    }
  }
  return phi_vort;
}
