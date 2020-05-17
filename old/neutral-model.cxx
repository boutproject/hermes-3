/*!
 *
 */

#include "neutral-model.hxx"

#include "diffusion2d.hxx"
#include "full-velocity.hxx"
#include "mixed.hxx"
#include "none.hxx"
#include "recycling.hxx"

using bout::globals::mesh;

NeutralModel *NeutralModel::create(Solver *solver, Mesh *mesh,
                                   Options &options) {
  // Decide which neutral model to use
  std::string type  = options["type"].withDefault<std::string>("none");

  if (type == "none") {
    // Neutral model which does nothing
    return NULL; // new NeutralNone(solver, mesh, options);
  } else if (type == "diffusion2d") {
    // Diffusion in X-Z only
    return new Diffusion2D(solver, mesh, options);
  } else if (type == "recycling") {
    // Recycling at target, assumes exponential neutral profile
    return new NeutralRecycling(solver, mesh, options);
  } else if (type == "fullvelocity") {
    // 3D Navier-Stokes
    return new FullVelocity(solver, mesh, options);
  } else if (type == "mixed") {
    // Diffusive in X-Z, fluid in Y. Similar to UEDGE
    return new NeutralMixed(solver, mesh, options);
  }
  throw BoutException("Unrecognised neutral model '%s'", type.c_str());
}

/*!
 * Atomic processes
 *
 * This code integrates atomic cross-sections over each grid cell
 *
 * NOTE: Currently this only integrates in Y, but should integrate in 3D
 */
void NeutralModel::neutral_rates(
    const Field3D &Ne, const Field3D &Te, const Field3D &Ti,
    const Field3D &Vi, // Plasma quantities
    const Field3D &Nn, const Field3D &Tn, const Field3D &Vnpar, // Neutral gas
    Field3D &S, Field3D &F, Field3D &Qi, Field3D &R, // Transfer rates
    Field3D &Riz, Field3D &Rrc, Field3D &Rcx) {      // Rates

  // Allocate output fields
  S = 0.0;
  F = 0.0;
  Qi = 0.0;
  R = 0.0;

  Riz = 0.0;
  Rrc = 0.0;
  Rcx = 0.0;

  Coordinates *coord = mesh->getCoordinates();

  for (int i = mesh->xstart; i <= mesh->xend; i++)
    for (int j = mesh->ystart; j <= mesh->yend; j++)
      for (int k = 0; k < mesh->LocalNz; k++) {
        // Integrate rates over each cell in Y
        // NOTE: This should integrate over (x,y,z)

        // Calculate cell centre (C), left (L) and right (R) values
        BoutReal Te_C = Te(i, j, k),
                 Te_L = 0.5 * (Te(i, j - 1, k) + Te(i, j, k)),
                 Te_R = 0.5 * (Te(i, j, k) + Te(i, j + 1, k));
        BoutReal Ti_C = Ti(i, j, k),
                 Ti_L = 0.5 * (Ti(i, j - 1, k) + Ti(i, j, k)),
                 Ti_R = 0.5 * (Ti(i, j, k) + Ti(i, j + 1, k));
        BoutReal Ne_C = Ne(i, j, k),
                 Ne_L = 0.5 * (Ne(i, j - 1, k) + Ne(i, j, k)),
                 Ne_R = 0.5 * (Ne(i, j, k) + Ne(i, j + 1, k));
        BoutReal Vi_C = Vi(i, j, k),
                 Vi_L = 0.5 * (Vi(i, j - 1, k) + Vi(i, j, k)),
                 Vi_R = 0.5 * (Vi(i, j, k) + Vi(i, j + 1, k));
        BoutReal Tn_C = Tn(i, j, k),
                 Tn_L = 0.5 * (Tn(i, j - 1, k) + Tn(i, j, k)),
                 Tn_R = 0.5 * (Tn(i, j, k) + Tn(i, j + 1, k));
        BoutReal Nn_C = Nn(i, j, k),
                 Nn_L = 0.5 * (Nn(i, j - 1, k) + Nn(i, j, k)),
                 Nn_R = 0.5 * (Nn(i, j, k) + Nn(i, j + 1, k));
        BoutReal Vn_C = Vnpar(i, j, k),
                 Vn_L = 0.5 * (Vnpar(i, j - 1, k) + Vnpar(i, j, k)),
                 Vn_R = 0.5 * (Vnpar(i, j, k) + Vnpar(i, j + 1, k));

        if (Ne_C < 0.)
          Ne_C = 0.0;
        if (Ne_L < 0.)
          Ne_L = 0.0;
        if (Ne_R < 0.)
          Ne_R = 0.0;
        if (Nn_C < 0.)
          Nn_C = 0.0;
        if (Nn_L < 0.)
          Nn_L = 0.0;
        if (Nn_R < 0.)
          Nn_R = 0.0;

        // Jacobian (Cross-sectional area)
        BoutReal J_C = coord->J(i, j),
                 J_L = 0.5 * (coord->J(i, j - 1) + coord->J(i, j)),
                 J_R = 0.5 * (coord->J(i, j) + coord->J(i, j + 1));

        ///////////////////////////////////////////
        // Charge exchange

        BoutReal R_cx_L = Ne_L * Nn_L * hydrogen.chargeExchange(Te_L * Tnorm) *
                          (Nnorm / Fnorm);
        BoutReal R_cx_C = Ne_C * Nn_C * hydrogen.chargeExchange(Te_C * Tnorm) *
                          (Nnorm / Fnorm);
        BoutReal R_cx_R = Ne_R * Nn_R * hydrogen.chargeExchange(Te_R * Tnorm) *
                          (Nnorm / Fnorm);

        // Power transfer from plasma to neutrals
        // Factor of 3/2 to convert temperature to energy

        Qi(i, j, k) = (3. / 2) * (J_L * (Ti_L - Tn_L) * R_cx_L +
                                  4. * J_C * (Ti_C - Tn_C) * R_cx_C +
                                  J_R * (Ti_R - Tn_R) * R_cx_R) /
                      (6. * J_C);

        // Plasma-neutral friction
        F(i, j, k) =
            (J_L * (Vi_L - Vn_L) * R_cx_L + 4. * J_C * (Vi_C - Vn_C) * R_cx_C +
             J_R * (Vi_R - Vn_R) * R_cx_R) /
            (6. * J_C);

        // Cell-averaged rate
        Rcx(i, j, k) =
            (J_L * R_cx_L + 4. * J_C * R_cx_C + J_R * R_cx_R) / (6. * J_C);

        ///////////////////////////////////////
        // Recombination

        BoutReal R_rc_L = hydrogen.recombination(Ne_L * Nnorm, Te_L * Tnorm) *
                          SQ(Ne_L) * Nnorm / Fnorm;
        BoutReal R_rc_C = hydrogen.recombination(Ne_C * Nnorm, Te_C * Tnorm) *
                          SQ(Ne_C) * Nnorm / Fnorm;
        BoutReal R_rc_R = hydrogen.recombination(Ne_R * Nnorm, Te_R * Tnorm) *
                          SQ(Ne_R) * Nnorm / Fnorm;

        // Radiated power from plasma
        // Factor of 1.09 so that recombination becomes an energy source at
        // 5.25eV
        R(i, j, k) = (J_L * (1.09 * Te_L - 13.6 / Tnorm) * R_rc_L +
                      4. * J_C * (1.09 * Te_C - 13.6 / Tnorm) * R_rc_C +
                      J_R * (1.09 * Te_R - 13.6 / Tnorm) * R_rc_R) /
                     (6. * J_C);

        // Plasma sink / neutral source
        S(i, j, k) =
            (J_L * R_rc_L + 4. * J_C * R_rc_C + J_R * R_rc_R) / (6. * J_C);

        // Transfer of ion momentum to neutrals
        F(i, j, k) += (J_L * Vi_L * R_rc_L + 4. * J_C * Vi_C * R_rc_C +
                       J_R * Vi_R * R_rc_R) /
                      (6. * J_C);

        // Transfer of ion energy to neutrals
        Qi(i, j, k) += (3. / 2) *
                       (J_L * Ti_L * R_rc_L + 4. * J_C * Ti_C * R_rc_C +
                        J_R * Ti_R * R_rc_R) /
                       (6. * J_C);

        // Cell-averaged rate
        Rrc(i, j, k) =
            (J_L * R_rc_L + 4. * J_C * R_rc_C + J_R * R_rc_R) / (6. * J_C);

        ///////////////////////////////////////
        // Ionisation

        BoutReal R_iz_L =
            Ne_L * Nn_L * hydrogen.ionisation(Te_L * Tnorm) * Nnorm / Fnorm;
        BoutReal R_iz_C =
            Ne_C * Nn_C * hydrogen.ionisation(Te_C * Tnorm) * Nnorm / Fnorm;
        BoutReal R_iz_R =
            Ne_R * Nn_R * hydrogen.ionisation(Te_R * Tnorm) * Nnorm / Fnorm;

        // Neutral sink, plasma source
        S(i, j, k) -=
            (J_L * R_iz_L + 4. * J_C * R_iz_C + J_R * R_iz_R) / (6. * J_C);

        // Transfer of neutral momentum to ions
        F(i, j, k) -= (J_L * Vn_L * R_iz_L + 4. * J_C * Vn_C * R_iz_C +
                       J_R * Vn_R * R_iz_R) /
                      (6. * J_C);

        // Transfer of neutral energy to ions
        Qi(i, j, k) -= (3. / 2) *
                       (J_L * Tn_L * R_iz_L + 4. * J_C * Tn_C * R_iz_C +
                        J_R * Tn_R * R_iz_R) /
                       (6. * J_C);

        // Ionisation and electron excitation energy
        R(i, j, k) += (Eionize / Tnorm) *
                      (J_L * R_iz_L + 4. * J_C * R_iz_C + J_R * R_iz_R) /
                      (6. * J_C);

        // Cell-averaged rate
        Riz(i, j, k) =
            (J_L * R_iz_L + 4. * J_C * R_iz_C + J_R * R_iz_R) / (6. * J_C);
      }
}
