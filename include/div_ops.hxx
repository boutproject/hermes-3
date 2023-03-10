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

#include <field3d.hxx>
#include <vector3d.hxx>
#include <bout/fv_ops.hxx>

/*!
 * Diffusion in index space
 * 
 * Similar to using Div_par_diffusion(SQ(mesh->dy)*mesh->g_22, f)
 *
 * @param[in] The field to be differentiated
 * @param[in] bndry_flux  Are fluxes through the boundary calculated?
 */
const Field3D Div_par_diffusion_index(const Field3D &f, bool bndry_flux=true);

const Field3D Div_n_bxGrad_f_B_XPPM(const Field3D &n, const Field3D &f, bool bndry_flux=true, bool poloidal=false, bool positive=false);

const Field3D Div_Perp_Lap_FV_Index(const Field3D &a, const Field3D &f, bool xflux);


// 4th-order flux conserving term, in index space
const Field3D D4DX4_FV_Index(const Field3D &f, bool bndry_flux=false);

// Div ( k * Grad(f) )
const Field2D Laplace_FV(const Field2D &k, const Field2D &f);

/// Perpendicular diffusion including X and Y directions
const Field3D Div_a_Grad_perp_upwind(const Field3D& a, const Field3D& f);

namespace FV {
  template<typename CellEdges = MC>
  const Field3D Div_par_fvv(const Field3D &f_in, const Field3D &v_in,
                            const Field3D &wave_speed_in, bool fixflux=true) {

    ASSERT1(areFieldsCompatible(f_in, v_in));
    ASSERT1(areFieldsCompatible(f_in, wave_speed_in));

    Mesh* mesh = f_in.getMesh();

    CellEdges cellboundary;

    /// Ensure that f, v and wave_speed are field aligned
    Field3D f = toFieldAligned(f_in, "RGN_NOX");
    Field3D v = toFieldAligned(v_in, "RGN_NOX");
    Field3D wave_speed = toFieldAligned(wave_speed_in, "RGN_NOX");

    Coordinates *coord = f_in.getCoordinates();

    Field3D result{zeroFrom(f)};
    
    // Only need one guard cell, so no need to communicate fluxes
    // Instead calculate in guard cells to preserve fluxes
    int ys = mesh->ystart-1;
    int ye = mesh->yend+1;

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
        BoutReal common_factor = (coord->J(i, j) + coord->J(i, j + 1)) /
          (sqrt(coord->g_22(i, j)) + sqrt(coord->g_22(i, j + 1)));
        
        BoutReal flux_factor_rc = common_factor / (coord->dy(i, j) * coord->J(i, j));
        BoutReal flux_factor_rp = common_factor / (coord->dy(i, j + 1) * coord->J(i, j + 1));

        // For left cell boundaries
        common_factor = (coord->J(i, j) + coord->J(i, j - 1)) /
          (sqrt(coord->g_22(i, j)) + sqrt(coord->g_22(i, j - 1)));

        BoutReal flux_factor_lc = common_factor / (coord->dy(i, j) * coord->J(i, j));
        BoutReal flux_factor_lm = common_factor / (coord->dy(i, j - 1) * coord->J(i, j - 1));
        
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
              flux = s.R * vpar * sv.R + wave_speed(i, j, k) * (s.R * sv.R - bndryval * vpar);
            }

          } else {

            // Maximum wave speed in the two cells
            BoutReal amax = BOUTMAX(wave_speed(i, j, k), wave_speed(i, j + 1, k));

            if (vpar > amax) {
              // Supersonic flow out of this cell
              flux = s.R * vpar * sv.R;
            } else if (vpar < -amax) {
              // Supersonic flow into this cell
              flux = 0.0;
            } else {
              // Subsonic flow, so a mix of right and left fluxes
              flux = s.R * 0.5 * (vpar + amax) * sv.R;
            }
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
              flux = s.L * vpar * sv.L - wave_speed(i, j, k) * (s.L * sv.L - bndryval * vpar);
            }
          } else {
            
            // Maximum wave speed in the two cells
            BoutReal amax = BOUTMAX(wave_speed(i, j, k), wave_speed(i, j - 1, k));

            if (vpar < -amax) {
              // Supersonic out of this cell
              flux = s.L * vpar * sv.L;
            } else if (vpar > amax) {
              // Supersonic into this cell
              flux = 0.0;
            } else {
              flux = s.L * 0.5 * (vpar - amax) * sv.L;
            }
          }
          
          result(i, j, k) -= flux * flux_factor_lc;
          result(i, j - 1, k) += flux * flux_factor_lm;
          
        }
      }
    }
    return fromFieldAligned(result, "RGN_NOBNDRY");
  }
  
}

#endif //  __DIV_OPS_H__
