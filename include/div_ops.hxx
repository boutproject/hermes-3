/*
  Finite volume discretisations of divergence operators

  ***********

    Copyright B.Dudson, J.Leddy, University of York, September 2016
              email: benjamin.dudson@york.ac.uk

    This file is part of Hermes.

    Hermes is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Hermes is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Hermes.  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef __DIV_OPS_H__
#define __DIV_OPS_H__

#include <bout/field3d.hxx>
#include <bout/vector3d.hxx>
#include <bout/fv_ops.hxx>

/*!
 * Diffusion in index space
 *
 * Similar to using Div_par_diffusion(SQ(mesh->dy)*mesh->g_22, f)
 *
 * @param[in] The field to be differentiated
 * @param[in] bndry_flux  Are fluxes through the boundary calculated?
 */
const Field3D Div_par_diffusion_index(const Field3D& f, bool bndry_flux = true);

const Field3D Div_n_bxGrad_f_B_XPPM(const Field3D& n, const Field3D& f,
                                    bool bndry_flux = true, bool poloidal = false,
                                    bool positive = false);

const Field3D Div_Perp_Lap_FV_Index(const Field3D& a, const Field3D& f, bool xflux);

const Field3D Div_Z_FV_Index(const Field3D& a, const Field3D& f);

// 4th-order flux conserving term, in index space
const Field3D D4DX4_FV_Index(const Field3D& f, bool bndry_flux = false);
const Field3D D4DZ4_Index(const Field3D& f);

// Div ( k * Grad(f) )
const Field2D Laplace_FV(const Field2D& k, const Field2D& f);

/// Perpendicular diffusion including X and Y directions
const Field3D Div_a_Grad_perp_upwind(const Field3D& a, const Field3D& f);

namespace FV {

/// Superbee limiter
///
/// This corresponds to the limiter function
///    φ(r) = max(0, min(2r, 1), min(r,2)
///
/// The value at cell right (i.e. i + 1/2) is:
///
///   n.R = n.c - φ(r) (n.c - (n.p + n.c)/2)
///       = n.c + φ(r) (n.p - n.c)/2
///
/// Four regimes:
///  a) r < 1/2 -> φ(r) = 2r
///     n.R = n.c + gL
///  b) 1/2 < r < 1 -> φ(r) = 1
///     n.R = n.c + gR/2
///  c) 1 < r < 2 -> φ(r) = r
///     n.R = n.c + gL/2
///  d) 2 < r  -> φ(r) = 2
///     n.R = n.c + gR
///
///  where the left and right gradients are:
///   gL = n.c - n.m
///   gR = n.p - n.c
///
struct Superbee {
  void operator()(Stencil1D& n) {
    BoutReal gL = n.c - n.L;
    BoutReal gR = n.R - n.c;

    // r = gL / gR
    // Limiter is φ(r)
    if (gL * gR < 0) {
      // Different signs => Zero gradient
      n.L = n.R = n.c;
    } else {
      BoutReal sign = SIGN(gL);
      gL = fabs(gL);
      gR = fabs(gR);
      BoutReal half_slope = sign * BOUTMAX(BOUTMIN(gL, 0.5*gR),
                                           BOUTMIN(gR, 0.5*gL));
      n.L = n.c - half_slope;
      n.R = n.c + half_slope;
    }
  }
};

template <typename CellEdges = MC>
const Field3D Div_par_fvv(const Field3D& f_in, const Field3D& v_in,
                          const Field3D& wave_speed_in, bool fixflux = true) {

  ASSERT1(areFieldsCompatible(f_in, v_in));
  ASSERT1(areFieldsCompatible(f_in, wave_speed_in));

  Mesh* mesh = f_in.getMesh();

  CellEdges cellboundary;

  /// Ensure that f, v and wave_speed are field aligned
  Field3D f = toFieldAligned(f_in, "RGN_NOX");
  Field3D v = toFieldAligned(v_in, "RGN_NOX");
  Field3D wave_speed = toFieldAligned(wave_speed_in, "RGN_NOX");

  Coordinates* coord = f_in.getCoordinates();

  Field3D result{zeroFrom(f)};

  // Only need one guard cell, so no need to communicate fluxes
  // Instead calculate in guard cells to preserve fluxes
  int ys = mesh->ystart - 1;
  int ye = mesh->yend + 1;

  for (int i = mesh->xstart; i <= mesh->xend; i++) {

    if (!mesh->firstY(i) || mesh->periodicY(i)) {
      // Calculate in guard cell to get fluxes consistent between processors
      ys = mesh->ystart - 1;
    } else {
      // Don't include the boundary cell. Note that this implies special
      // handling of boundaries later
      ys = mesh->ystart;
    }

    if (!mesh->lastY(i) || mesh->periodicY(i)) {
      // Calculate in guard cells
      ye = mesh->yend + 1;
    } else {
      // Not in boundary cells
      ye = mesh->yend;
    }

    for (int j = ys; j <= ye; j++) {
      // Pre-calculate factors which multiply fluxes

      // For right cell boundaries
      BoutReal common_factor = (coord->J(i, j) + coord->J(i, j + 1))
                               / (sqrt(coord->g_22(i, j)) + sqrt(coord->g_22(i, j + 1)));

      BoutReal flux_factor_rc = common_factor / (coord->dy(i, j) * coord->J(i, j));
      BoutReal flux_factor_rp =
          common_factor / (coord->dy(i, j + 1) * coord->J(i, j + 1));

      // For left cell boundaries
      common_factor = (coord->J(i, j) + coord->J(i, j - 1))
                      / (sqrt(coord->g_22(i, j)) + sqrt(coord->g_22(i, j - 1)));

      BoutReal flux_factor_lc = common_factor / (coord->dy(i, j) * coord->J(i, j));
      BoutReal flux_factor_lm =
          common_factor / (coord->dy(i, j - 1) * coord->J(i, j - 1));

      for (int k = 0; k < mesh->LocalNz; k++) {

        ////////////////////////////////////////////
        // Reconstruct f at the cell faces
        // This calculates s.R and s.L for the Right and Left
        // face values on this cell

        // Reconstruct f at the cell faces
        Stencil1D s;
        s.c = f(i, j, k);
        s.m = f(i, j - 1, k);
        s.p = f(i, j + 1, k);

        cellboundary(s); // Calculate s.R and s.L

        // Reconstruct v at the cell faces
        Stencil1D sv;
        sv.c = v(i, j, k);
        sv.m = v(i, j - 1, k);
        sv.p = v(i, j + 1, k);

        cellboundary(sv);

        ////////////////////////////////////////////
        // Right boundary

        // Calculate velocity at right boundary (y+1/2)
        BoutReal vpar = 0.5 * (v(i, j, k) + v(i, j + 1, k));
        BoutReal flux;

        if (mesh->lastY(i) && (j == mesh->yend) && !mesh->periodicY(i)) {
          // Last point in domain

          BoutReal bndryval = 0.5 * (s.c + s.p);
          if (fixflux) {
            // Use mid-point to be consistent with boundary conditions
            flux = bndryval * vpar * vpar;
          } else {
            // Add flux due to difference in boundary values
            flux =
                s.R * vpar * sv.R
                + BOUTMAX(wave_speed(i, j, k), fabs(v(i, j, k)), fabs(v(i, j + 1, k)))
                  * (s.R * sv.R - bndryval * vpar);
          }
        } else {
          // Maximum wave speed in the two cells
          BoutReal amax = BOUTMAX(wave_speed(i, j, k), wave_speed(i, j + 1, k),
                                  fabs(v(i, j, k)), fabs(v(i, j + 1, k)));

          flux = s.R * 0.5 * (sv.R + amax) * sv.R;
        }

        result(i, j, k) += flux * flux_factor_rc;
        result(i, j + 1, k) -= flux * flux_factor_rp;

        ////////////////////////////////////////////
        // Calculate at left boundary

        vpar = 0.5 * (v(i, j, k) + v(i, j - 1, k));

        if (mesh->firstY(i) && (j == mesh->ystart) && !mesh->periodicY(i)) {
          // First point in domain
          BoutReal bndryval = 0.5 * (s.c + s.m);
          if (fixflux) {
            // Use mid-point to be consistent with boundary conditions
            flux = bndryval * vpar * vpar;
          } else {
            // Add flux due to difference in boundary values
            flux =
                s.L * vpar * sv.L
                - BOUTMAX(wave_speed(i, j, k), fabs(v(i, j, k)), fabs(v(i, j - 1, k)))
                  * (s.L * sv.L - bndryval * vpar);
          }
        } else {
          // Maximum wave speed in the two cells
          BoutReal amax = BOUTMAX(wave_speed(i, j, k), wave_speed(i, j - 1, k),
                                  fabs(v(i, j, k)), fabs(v(i, j - 1, k)));

          flux = s.L * 0.5 * (sv.L - amax) * sv.L;
        }

        result(i, j, k) -= flux * flux_factor_lc;
        result(i, j - 1, k) += flux * flux_factor_lm;
      }
    }
  }
  return fromFieldAligned(result, "RGN_NOBNDRY");
}
/// Finite volume parallel divergence
///
/// NOTE: Modified version, applies limiter to velocity and field
///       Performs better (smaller overshoots) than Div_par
///
/// Preserves the sum of f*J*dx*dy*dz over the domain
///
/// @param[in] f_in   The field being advected.
///                   This will be reconstructed at cell faces
///                   using the given CellEdges method
/// @param[in] v_in   The advection velocity.
///                   This will be interpolated to cell boundaries
///                   using linear interpolation
/// @param[in] wave_speed_in  Local maximum speed of all waves in the system at each
//                            point in space
/// @param[in] fixflux     Fix the flux at the boundary to be the value at the
///                        midpoint (for boundary conditions)
///
/// NB: Uses to/from FieldAligned coordinates
template <typename CellEdges = MC>
const Field3D Div_par_mod(const Field3D& f_in, const Field3D& v_in,
                          const Field3D& wave_speed_in, bool fixflux = true) {

  ASSERT1_FIELDS_COMPATIBLE(f_in, v_in);
  ASSERT1_FIELDS_COMPATIBLE(f_in, wave_speed_in);

  Mesh* mesh = f_in.getMesh();

  CellEdges cellboundary;

  ASSERT2(f_in.getDirectionY() == v_in.getDirectionY());
  ASSERT2(f_in.getDirectionY() == wave_speed_in.getDirectionY());
  const bool are_unaligned =
      ((f_in.getDirectionY() == YDirectionType::Standard)
       and (v_in.getDirectionY() == YDirectionType::Standard)
       and (wave_speed_in.getDirectionY() == YDirectionType::Standard));

  Field3D f = are_unaligned ? toFieldAligned(f_in, "RGN_NOX") : f_in;
  Field3D v = are_unaligned ? toFieldAligned(v_in, "RGN_NOX") : v_in;
  Field3D wave_speed =
      are_unaligned ? toFieldAligned(wave_speed_in, "RGN_NOX") : wave_speed_in;

  Coordinates* coord = f_in.getCoordinates();

  Field3D result{zeroFrom(f)};

  // Only need one guard cell, so no need to communicate fluxes
  // Instead calculate in guard cells to preserve fluxes
  int ys = mesh->ystart - 1;
  int ye = mesh->yend + 1;

  for (int i = mesh->xstart; i <= mesh->xend; i++) {

    if (!mesh->firstY(i) || mesh->periodicY(i)) {
      // Calculate in guard cell to get fluxes consistent between processors
      ys = mesh->ystart - 1;
    } else {
      // Don't include the boundary cell. Note that this implies special
      // handling of boundaries later
      ys = mesh->ystart;
    }

    if (!mesh->lastY(i) || mesh->periodicY(i)) {
      // Calculate in guard cells
      ye = mesh->yend + 1;
    } else {
      // Not in boundary cells
      ye = mesh->yend;
    }

    for (int j = ys; j <= ye; j++) {
      // Pre-calculate factors which multiply fluxes
#if not(BOUT_USE_METRIC_3D)
      // For right cell boundaries
      BoutReal common_factor = (coord->J(i, j) + coord->J(i, j + 1))
                               / (sqrt(coord->g_22(i, j)) + sqrt(coord->g_22(i, j + 1)));

      BoutReal flux_factor_rc = common_factor / (coord->dy(i, j) * coord->J(i, j));
      BoutReal flux_factor_rp =
          common_factor / (coord->dy(i, j + 1) * coord->J(i, j + 1));

      // For left cell boundaries
      common_factor = (coord->J(i, j) + coord->J(i, j - 1))
                      / (sqrt(coord->g_22(i, j)) + sqrt(coord->g_22(i, j - 1)));

      BoutReal flux_factor_lc = common_factor / (coord->dy(i, j) * coord->J(i, j));
      BoutReal flux_factor_lm =
          common_factor / (coord->dy(i, j - 1) * coord->J(i, j - 1));
#endif
      for (int k = 0; k < mesh->LocalNz; k++) {
#if BOUT_USE_METRIC_3D
        // For right cell boundaries
        BoutReal common_factor =
            (coord->J(i, j, k) + coord->J(i, j + 1, k))
            / (sqrt(coord->g_22(i, j, k)) + sqrt(coord->g_22(i, j + 1, k)));

        BoutReal flux_factor_rc =
            common_factor / (coord->dy(i, j, k) * coord->J(i, j, k));
        BoutReal flux_factor_rp =
            common_factor / (coord->dy(i, j + 1, k) * coord->J(i, j + 1, k));

        // For left cell boundaries
        common_factor = (coord->J(i, j, k) + coord->J(i, j - 1, k))
                        / (sqrt(coord->g_22(i, j, k)) + sqrt(coord->g_22(i, j - 1, k)));

        BoutReal flux_factor_lc =
            common_factor / (coord->dy(i, j, k) * coord->J(i, j, k));
        BoutReal flux_factor_lm =
            common_factor / (coord->dy(i, j - 1, k) * coord->J(i, j - 1, k));
#endif

        ////////////////////////////////////////////
        // Reconstruct f at the cell faces
        // This calculates s.R and s.L for the Right and Left
        // face values on this cell

        // Reconstruct f at the cell faces
        Stencil1D s;
        s.c = f(i, j, k);
        s.m = f(i, j - 1, k);
        s.p = f(i, j + 1, k);

        cellboundary(s); // Calculate s.R and s.L

        ////////////////////////////////////////////
        // Reconstruct v at the cell faces
        Stencil1D sv;
        sv.c = v(i, j, k);
        sv.m = v(i, j - 1, k);
        sv.p = v(i, j + 1, k);

        cellboundary(sv); // Calculate sv.R and sv.L

        ////////////////////////////////////////////
        // Right boundary

        BoutReal flux;

        if (mesh->lastY(i) && (j == mesh->yend) && !mesh->periodicY(i)) {
          // Last point in domain

          // Calculate velocity at right boundary (y+1/2)
          BoutReal vpar = 0.5 * (v(i, j, k) + v(i, j + 1, k));

          BoutReal bndryval = 0.5 * (s.c + s.p);
          if (fixflux) {
            // Use mid-point to be consistent with boundary conditions
            flux = bndryval * vpar;
          } else {
            // Add flux due to difference in boundary values
            flux = s.R * vpar + wave_speed(i, j, k) * (s.R - bndryval);
          }
        } else {

          // Maximum wave speed in the two cells
          BoutReal amax = BOUTMAX(wave_speed(i, j, k), wave_speed(i, j + 1, k),
                                  fabs(v(i, j, k)), fabs(v(i, j + 1, k)));

          flux = s.R * 0.5 * (sv.R + amax);
        }

        result(i, j, k) += flux * flux_factor_rc;
        result(i, j + 1, k) -= flux * flux_factor_rp;

        ////////////////////////////////////////////
        // Calculate at left boundary

        if (mesh->firstY(i) && (j == mesh->ystart) && !mesh->periodicY(i)) {
          // First point in domain
          BoutReal bndryval = 0.5 * (s.c + s.m);
          BoutReal vpar = 0.5 * (v(i, j, k) + v(i, j - 1, k));
          if (fixflux) {
            // Use mid-point to be consistent with boundary conditions
            flux = bndryval * vpar;
          } else {
            // Add flux due to difference in boundary values
            flux = s.L * vpar - wave_speed(i, j, k) * (s.L - bndryval);
          }
        } else {

          // Maximum wave speed in the two cells
          BoutReal amax = BOUTMAX(wave_speed(i, j, k), wave_speed(i, j - 1, k),
                                  fabs(v(i, j, k)), fabs(v(i, j - 1, k)));

          flux = s.L * 0.5 * (sv.L - amax);
        }

        result(i, j, k) -= flux * flux_factor_lc;
        result(i, j - 1, k) += flux * flux_factor_lm;
      }
    }
  }
  return are_unaligned ? fromFieldAligned(result, "RGN_NOBNDRY") : result;
}

} // namespace FV

#endif //  __DIV_OPS_H__
