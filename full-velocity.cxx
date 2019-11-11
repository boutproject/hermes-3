
#include "full-velocity.hxx"

#include <boutexception.hxx>
#include <options.hxx>

#include "div_ops.hxx"

using bout::globals::mesh;

FullVelocity::FullVelocity(Solver *solver, Mesh *mesh, Options &options)
    : NeutralModel(options) {
  /*! 2D (X-Y) full velocity model
   *
   * Vn2D is covariant (the default), so has components
   * V_x, V_y, V_z
   */

  coord = mesh->getCoordinates();
  
  options.get("gamma_ratio", gamma_ratio, 5. / 3);
  options.get("viscosity", neutral_viscosity, 1e-2);
  options.get("bulk", neutral_bulk, 1e-2);
  options.get("conduction", neutral_conduction, 1e-2);

  OPTION(options, outflow_ydown, false); // Allow outflowing neutrals?

  OPTION(options, neutral_gamma, 5. / 4);

  // Evolve 2D density, pressure, and velocity
  solver->add(Nn2D, "Nn");
  solver->add(Pn2D, "Pn");
  solver->add(Vn2D, "Vn");

  DivV2D.setBoundary("Pn"); // Same boundary condition as Pn
  SAVE_REPEAT(DivV2D);

  // Load necessary metrics for non-orth calculation
  Field2D etaxy, cosbeta;
  if (mesh->get(etaxy, "etaxy")) {
    etaxy = 0.0;
  }
  cosbeta = sqrt(1. - SQ(etaxy));

  // Calculate transformation to Cartesian coordinates
  Field2D Rxy, Zxy, hthe, Bpxy;

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
  Rxy /= Lnorm;
  Zxy /= Lnorm;
  hthe /= Lnorm;
  Bpxy /= Bnorm;

  // Axisymmetric neutrals simplifies things considerably...

  Urx.allocate();
  Ury.allocate();
  Uzx.allocate();
  Uzy.allocate();

  Txr.allocate();
  Txz.allocate();
  Tyr.allocate();
  Tyz.allocate();

  Coordinates *coord = mesh->getCoordinates();
  
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
      BoutReal dldpsi =
          sqrt(SQ(dRdpsi) + SQ(dZdpsi)) * cosbeta(i, j); // ~ 1/(R*Bp)
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

  SAVE_ONCE(Urx, Ury, Uzx, Uzy);
  SAVE_ONCE(Txr, Txz, Tyr, Tyz);

  // Atomic processes
  S = 0;
  F = 0;
  Qi = 0;
  Rp = 0;
}

void FullVelocity::update(const Field3D &Ne, const Field3D &Te,
                          const Field3D &Ti, const Field3D &Vi) {

  mesh->communicate(Nn2D, Vn2D, Pn2D);
  
  // Navier-Stokes for axisymmetric neutral gas profiles
  // Nn2D, Pn2D and Tn2D are unfloored
  Tn2D = Pn2D / Nn2D;

  // Nn and Tn are 3D floored fields, used for rate coefficients
  Field3D Nn = floor(Nn2D, 1e-8);
  Field3D Tn = Pn2D / Nn;
  Tn = floor(Tn, 0.01 / Tnorm);

  //////////////////////////////////////////////////////
  // 2D (X-Y) full velocity model
  //
  // Evolves density Nn2D, velocity vector Vn2D and pressure Pn2D
  //

  if (outflow_ydown) {
    // Outflowing boundaries at ydown. If flow direction is
    // into domain then zero value is set. If flow is out of domain
    // then Neumann conditions are set

    for (RangeIterator idwn = mesh->iterateBndryLowerY(); !idwn.isDone();
         idwn.next()) {

      if (Vn2D.y(idwn.ind, mesh->ystart) < 0.0) {
        // Flowing out of domain
        Vn2D.y(idwn.ind, mesh->ystart - 1) = Vn2D.y(idwn.ind, mesh->ystart);
      } else {
        // Flowing into domain
        Vn2D.y(idwn.ind, mesh->ystart - 1) = -Vn2D.y(idwn.ind, mesh->ystart);
      }
      // Neumann boundary condition on X and Z components
      Vn2D.x(idwn.ind, mesh->ystart - 1) = Vn2D.x(idwn.ind, mesh->ystart);
      Vn2D.z(idwn.ind, mesh->ystart - 1) = Vn2D.z(idwn.ind, mesh->ystart);

      // Neumann conditions on density and pressure
      Nn2D(idwn.ind, mesh->ystart - 1) = Nn2D(idwn.ind, mesh->ystart);
      Pn2D(idwn.ind, mesh->ystart - 1) = Pn2D(idwn.ind, mesh->ystart);
    }
  }

  // Density
  ddt(Nn2D) = -Div(Vn2D, Nn2D);

  Field2D Nn2D_floor = floor(Nn2D, 1e-2);
  // Velocity
  ddt(Vn2D) = -Grad(Pn2D) / Nn2D_floor;

  //////////////////////////////////////////////////////
  // Momentum advection

  // Convert to cylindrical coordinates for velocity
  // advection term. This is to avoid Christoffel symbol
  // terms in curvilinear geometry
  Field2D vr = Txr * Vn2D.x + Tyr * Vn2D.y; // Grad R component
  Field2D vz = Txz * Vn2D.x + Tyz * Vn2D.y; // Grad Z component

  // Advect as scalars (no Christoffel symbols needed)
  ddt(vr) = -V_dot_Grad(Vn2D, vr);
  ddt(vz) = -V_dot_Grad(Vn2D, vz);

  // Convert back to field-aligned coordinates
  ddt(Vn2D).x += Urx * ddt(vr) + Uzx * ddt(vz);
  ddt(Vn2D).y += Ury * ddt(vr) + Uzy * ddt(vz);

  //////////////////////////////////////////////////////
  // Viscosity
  // This includes dynamic ( neutral_viscosity)
  // and bulk/volume viscosity ( neutral_bulk )

  ddt(vr) = Laplace_FV(neutral_viscosity, vr);
  ddt(vz) = Laplace_FV(neutral_viscosity, vz);

  ddt(Vn2D).x += Urx * ddt(vr) + Uzx * ddt(vz);
  ddt(Vn2D).y += Ury * ddt(vr) + Uzy * ddt(vz);

  DivV2D = Div(Vn2D);
  DivV2D.applyBoundary(0.0);
  mesh->communicate(DivV2D);

  // ddt(Vn2D) += Grad( (neutral_viscosity/3. + neutral_bulk) * DivV2D ) /
  // Nn2D_floor;

  //////////////////////////////////////////////////////
  // Pressure
  ddt(Pn2D) = -Div(Vn2D, Pn2D) -
              (gamma_ratio - 1.) * Pn2D * DivV2D * floor(Nn2D, 0) / Nn2D_floor +
              Laplace_FV(neutral_conduction, Pn2D / Nn2D);

  ///////////////////////////////////////////////////////////////////
  // Boundary condition on fluxes

  for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {

    // Loss of thermal energy to the target.
    // This depends on the reflection coefficient
    // and is controlled by the option neutral_gamma
    //         q = neutral_gamma * n * T * cs

    // Density at the target
    BoutReal Nnout =
        0.5 * (Nn2D(r.ind, mesh->ystart) + Nn2D(r.ind, mesh->ystart - 1));
    if (Nnout < 0.0)
      Nnout = 0.0;
    // Temperature at the target
    BoutReal Tnout =
        0.5 * (Tn2D(r.ind, mesh->ystart) + Tn2D(r.ind, mesh->ystart - 1));
    if (Tnout < 0.0)
      Tnout = 0.0;

    // gamma * n * T * cs
    BoutReal q = neutral_gamma * Nnout * Tnout * sqrt(Tnout);
    // Multiply by cell area to get power
    BoutReal heatflux =
        q * (coord->J(r.ind, mesh->ystart) + coord->J(r.ind, mesh->ystart - 1)) /
        (sqrt(coord->g_22(r.ind, mesh->ystart)) +
         sqrt(coord->g_22(r.ind, mesh->ystart - 11)));

    // Divide by volume of cell, and multiply by 2/3 to get pressure
    ddt(Pn2D)(r.ind, mesh->ystart) -=
        (2. / 3) * heatflux /
        (coord->dy(r.ind, mesh->ystart) * coord->J(r.ind, mesh->ystart));
  }

  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {

    // Loss of thermal energy to the target.
    // This depends on the reflection coefficient
    // and is controlled by the option neutral_gamma
    //         q = neutral_gamma * n * T * cs

    // Density at the target
    BoutReal Nnout =
        0.5 * (Nn2D(r.ind, mesh->yend) + Nn2D(r.ind, mesh->yend + 1));
    if (Nnout < 0.0)
      Nnout = 0.0;
    // Temperature at the target
    BoutReal Tnout =
        0.5 * (Tn2D(r.ind, mesh->yend) + Tn2D(r.ind, mesh->yend + 1));
    if (Tnout < 0.0)
      Tnout = 0.0;

    // gamma * n * T * cs
    BoutReal q = neutral_gamma * Nnout * Tnout * sqrt(Tnout);
    // Multiply by cell area to get power
    BoutReal heatflux =
        q * (coord->J(r.ind, mesh->yend) + coord->J(r.ind, mesh->yend + 1)) /
        (sqrt(coord->g_22(r.ind, mesh->yend)) +
         sqrt(coord->g_22(r.ind, mesh->yend + 1)));

    // Divide by volume of cell, and multiply by 2/3 to get pressure
    ddt(Pn2D)(r.ind, mesh->yend) -=
        (2. / 3) * heatflux /
        (coord->dy(r.ind, mesh->yend) * coord->J(r.ind, mesh->yend));
  }

  // Exchange of parallel momentum. This could be done
  // in a couple of ways, but here we use the fact that
  // Vn2D is covariant and b = e_y / (JB) to write:
  //
  // V_{||n} = b dot V_n = Vn2D.y / (JB)
  Field2D Vnpar = Vn2D.y / (coord->J * coord->Bxy);

  /////////////////////////////////////////////////////
  // Atomic processes

  Field3D Riz, Rrc, Rcx;
  neutral_rates(Ne, Te, Ti, Vi, Nn, Tn, Vnpar, S, F, Qi, Rp, Riz, Rrc, Rcx);

  Fperp = Rrc + Rcx; // Friction for vorticity

  // Loss of momentum in the X and Z directions
  ddt(Vn2D).x -= (DC(Rcx) + DC(Riz)) * Vn2D.x / Nn2D_floor;
  ddt(Vn2D).z -= (DC(Rcx) + DC(Riz)) * Vn2D.z / Nn2D_floor;

  // Particles
  ddt(Nn2D) += DC(S); // Average over toroidal angle z

  // Energy
  ddt(Pn2D) += (2. / 3) * DC(Qi);

  // Momentum. Note need to turn back into covariant form
  ddt(Vn2D).y += DC(F) * (coord->J * coord->Bxy) / Nn2D_floor;

  // Density evolution
  for (auto &i : Nn2D.getRegion("RGN_ALL")) {
    if ((Nn2D[i] < 1e-8) && (ddt(Nn2D)[i] < 0.0)) {
      ddt(Nn2D)[i] = 0.0;
    }
  }
}

void FullVelocity::addDensity(int x, int y, int UNUSED(z), BoutReal dndt) {
  ddt(Nn2D)(x, y) += dndt / mesh->LocalNz; // Average over Z, so Z index unused
}

void FullVelocity::addPressure(int x, int y, int UNUSED(z), BoutReal dpdt) {
  ddt(Pn2D)(x, y) += dpdt / mesh->LocalNz;
}

void FullVelocity::addMomentum(int x, int y, int UNUSED(z), BoutReal dnvdt) {
  // Momentum added in the direction of the magnetic field
  // Vn2D is covariant and b = e_y / (JB) to write:
  //
  // V_{||n} = b dot V_n = Vn2D.y / (JB)

  ddt(Vn2D.y)(x, y) += dnvdt * (coord->J(x, y) * coord->Bxy (x, y)) / Nn2D(x, y); 
}
