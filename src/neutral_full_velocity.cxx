
#include <bout/constants.hxx>
#include <bout/output_bout_types.hxx>

#include "../include/div_ops.hxx"
#include "../include/hermes_utils.hxx"
#include "../include/neutral_full_velocity.hxx"

#include "bout/mesh.hxx"
#include "bout/solver.hxx"

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
  const BoutReal Tnorm = units["eV"];

  auto& options = alloptions[name];

  AA = options["AA"].doc("Atomic mass number. Proton = 1").as<int>();

  adiabatic_index =
      options["adiabatic_index"]
          .doc("Ratio of specific heats γ = Cp/Cv [5/3 for monatomic ideal gas]")
          .withDefault(5. / 3);

  dynamic_coefficients = options["dynamic_coefficients"]
    .doc("Recalculate conduction coeffient?")
    .withDefault<bool>(false);

  if (!dynamic_coefficients) {
    // Fixed conduction coefficient

    neutral_conduction =
      options["conduction"].doc("Heat conduction [m^2/s]").withDefault(1.0)
      / (meters * meters / seconds);
    neutral_conduction.applyBoundary("dirichlet");
  }

  neutral_viscosity =
      options["viscosity"].doc("Kinematic viscosity [m^2/s]").withDefault(1.0)
      / (meters * meters / seconds);

  neutral_gamma = options["neutral_gamma"]
                      .doc("Surface heat transmission coefficient")
                      .withDefault(5. / 4);

  density_floor = options["density_floor"]
                      .doc("A minimum density used when dividing by density."
                           "Normalised units.")
                      .withDefault(1e-5);

  temperature_floor = options["temperature_floor"]
                          .doc("Low temperature scale")
                          .withDefault<BoutReal>(0.1)
                      / get<BoutReal>(alloptions["units"]["eV"]);

  precondition = options["precondition"]
                     .doc("Enable preconditioning in neutral model?")
                     .withDefault<bool>(false);

  if (precondition) {
    laplacian = std::unique_ptr<Laplacian>(Laplacian::create(&options["precon_laplace"]));

    laplacian->setInnerBoundaryFlags(INVERT_DC_GRAD | INVERT_AC_GRAD);
    laplacian->setOuterBoundaryFlags(INVERT_DC_GRAD | INVERT_AC_GRAD);
  }

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

  // Non-zero floor can use a differentiable soft floor
  Field2D Nnlim = softFloor(Nn2D, density_floor);

  // Equation of state at low density
  //
  // e = Cv T Nlim / N
  // p = N T
  //
  // The internal energy evolution of (N * e) is therefore
  // evolving Pn2D_solver = Nlim * T
  // rather than pressure Pn2D = N * T

  Tn2D = Pn2D / Nnlim;
  Tn2D.applyBoundary("neumann");
  Pn2D_solver = Pn2D; // Save solver variable to restore later
  Pn2D = Tn2D * Nn2D; // Equation of state

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

  Field2D Nn2D_floor = softFloor(Nn2D, density_floor);
  Field2D Tn2D_floor = softFloor(Tn2D, temperature_floor);
  Field2D Pn2D_floor = Nn2D_floor * Tn2D_floor;

  if (dynamic_coefficients) {
    if (!localstate.isSet("collision_frequency")) {
      throw BoutException("dynamic_coefficients requires collision_frequency");
    }
    Field2D collision_freq = DC(get<Field3D>(localstate["collision_frequency"]));

    // Temperature diffusion
    // Note: This is chi_n = (5/2) * Tn / (m * nu)
    //       where nu is the collision frequency
    neutral_conduction = (5. / 2) * Tn2D / (AA * collision_freq);

    // Set to zero at the boundaries
    neutral_conduction.applyBoundary("dirichlet");

    // Kinematic viscosity
    // Relationship between heat conduction and viscosity for neutral
    // gas Chapman, Cowling "The Mathematical Theory of Non-Uniform
    // Gases", CUP 1952 Ferziger, Kaper "Mathematical Theory of
    // Transport Processes in Gases", 1972
    //neutral_viscosity = (2. / 5) * neutral_conduction;
  }

  // Density
  ddt(Nn2D) = -Div(Vn2D_contravariant, Nn2D);

  // Velocity
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
  // vr.applyBoundary("dirichlet");
  // vz.applyBoundary("dirichlet");

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
  Field2D mass_density = AA * Nn2D_floor;
  Field2D dynamic_viscosity = AA * Nn2D * neutral_viscosity;
  dynamic_viscosity.applyBoundary("neumann");
  Field2D momentum_flux_limit = Pn2D_floor;

  ddt(vr) += Laplace_FV(dynamic_viscosity, vr, momentum_flux_limit) / mass_density;
  ddt(vz) += Laplace_FV(dynamic_viscosity, vz, momentum_flux_limit) / mass_density;

  // Convert back to field-aligned coordinates
  ddt(Vn2D).x += Urx * ddt(vr) + Uzx * ddt(vz);
  ddt(Vn2D).y += Ury * ddt(vr) + Uzy * ddt(vz);

  if (toroidal_flow) {
    Field2D vphi = Vn2D.z / (sigma_Bp * Rxy); // Toroidal component
    vphi.applyBoundary("neumann");

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
    ddt(vphi) += Laplace_FV(dynamic_viscosity, vphi, momentum_flux_limit) / mass_density;

    // Convert back to field-aligned
    ddt(Vn2D).z += ddt(vphi) * sigma_Bp * Rxy;
  }

  //////////////////////////////////////////////////////
  // Pressure
  // Note that the equation of state is modified by density floor
  //

  Field2D heat_flux_limit = Pn2D_floor * sqrt(Tn2D_floor * (2 / (PI * AA)));

  ddt(Pn2D) =
      -Div(Vn2D, Nn2D_floor * Tn2D + (adiabatic_index - 1.) * Pn2D)
      + (adiabatic_index - 1.) * (Vn2D_contravariant * GradPn2D)
    + (adiabatic_index - 1.) * Laplace_FV(Nn2D * neutral_conduction, Tn2D, heat_flux_limit);

  ///////////////////////////////////////////////////////////////////
  // Boundary condition on fluxes

  for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {

    // Loss of thermal energy to the target.
    // This depends on the reflection coefficient
    // and is controlled by the option neutral_gamma
    //         q = neutral_gamma * n * T * cs

    // Density at the target
    BoutReal Nnout = 0.5 * (Nn2D(r.ind, mesh->ystart) + Nn2D(r.ind, mesh->ystart - 1));
    if (Nnout < 0.0)
      Nnout = 0.0;
    // Temperature at the target
    BoutReal Tnout = 0.5 * (Tn2D(r.ind, mesh->ystart) + Tn2D(r.ind, mesh->ystart - 1));
    if (Tnout < 0.0)
      Tnout = 0.0;

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
    BoutReal Nnout = 0.5 * (Nn2D(r.ind, mesh->yend) + Nn2D(r.ind, mesh->yend + 1));
    if (Nnout < 0.0)
      Nnout = 0.0;
    // Temperature at the target
    BoutReal Tnout = 0.5 * (Tn2D(r.ind, mesh->yend) + Tn2D(r.ind, mesh->yend + 1));
    if (Tnout < 0.0)
      Tnout = 0.0;

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

  /////////////////////////////////////////////////////
  // Atomic processes

  // Particles
  if (localstate.isSet("density_source")) {
    ddt(Nn2D) += DC(get<Field3D>(localstate["density_source"]));
  }

  // Momentum. Note need to turn back into covariant form
  // The parallel and toroidal components partly cancel when converting
  // velocity to contravariant form for advection in X-Y.
  if (localstate.isSet("momentum_source")) {
    Field2D Fpar_mN = DC(get<Field3D>(localstate["momentum_source"])) / (AA * Nn2D_floor);
    ddt(Vn2D).y += Fpar_mN * (coord->J * coord->Bxy); // Parallel flow

    if (toroidal_flow) {
      ddt(Vn2D).z += Fpar_mN * coord->g_23 / (coord->J * coord->Bxy); // Toroidal flow
    }
  }

  if (localstate.isSet("collision_frequency")) {
    // Damp flow perpendicular to B due to collisions
    Field2D collision_freq = DC(get<Field3D>(localstate["collision_frequency"]));
    // Radial flow
    ddt(Vn2D).x -= Vn2D.x * collision_freq;
    // Binormal flow
    ddt(Vn2D).z -= (Vn2D.z - (coord->g_23 / coord->g_22) * Vn2D.y) * collision_freq;
  }

  // Energy
  if (localstate.isSet("energy_source")) {
    ddt(Pn2D) += (adiabatic_index - 1) * DC(get<Field3D>(localstate["energy_source"]));
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

  // Restore solver Pn2D so that restart files have the correcte value
  Pn2D = Pn2D_solver;
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

    if (dynamic_coefficients) {
      set_with_attrs(state[std::string("chi_") + name], neutral_conduction,
                     {{"time_dimension", "t"},
                      {"units", "m^2/s"},
                      {"conversion", Cs0 * Cs0 / Omega_ci},
                      {"standard_name", "conduction coefficient"},
                      {"long_name", name + " heat conduction coefficient"},
                      {"source", "neutral_full_velocity"}});
    } else {
      set_with_attrs(state[std::string("chi_") + name], neutral_conduction,
                     {{"units", "m^2/s"},
                      {"conversion", Cs0 * Cs0 / Omega_ci},
                      {"standard_name", "conduction coefficient"},
                      {"long_name", name + " heat conduction coefficient"},
                      {"source", "neutral_full_velocity"}});
    }
  }
}

void NeutralFullVelocity::precon(const Options& state, BoutReal gamma) {
  if (!precondition) {
    return;
  }

  Field2D Nn2D_floor = softFloor(Nn2D, density_floor);
  Field2D kappa = Nn2D * neutral_conduction;
  Tn2D.applyBoundary("neumann");

  // Solve cross-field heat diffusion
  //
  // dP/dt = (g - 1) Div(kappa Grad(P/N))
  //
  // The preconditioner inverts an operator:
  //
  // (I - dt J) = x - dt (g - 1)Div(kappa Grad(x / N))
  //
  // To remove the inner N we rearrange this to
  //
  // (N x' - dt (g - 1) Div(kappa Grad(x')))
  //
  // where x' = x / N
  //
  // The coefficients in the solver are:
  //
  //   d Laplace_perp(x) + a x + (1/c1)Grad(c2) dot Grad_perp(x) = b
  //
  Field2D logTn = log(Tn2D);
  laplacian->setCoefA(Nn2D_floor - gamma * (adiabatic_index - 1.) * FV::Div_a_Grad_perp(kappa, logTn));
  laplacian->setCoefC1(-1. / (gamma * (adiabatic_index - 1.)));
  laplacian->setCoefC2(kappa * logTn);
  laplacian->setCoefD((-gamma * (adiabatic_index - 1.)) * kappa);

  // This factor of N converts from x' to x
  ddt(Pn2D) = laplacian->solve(ddt(Pn2D)) * Nn2D_floor;
}
