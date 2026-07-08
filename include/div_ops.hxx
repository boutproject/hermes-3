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

#ifndef DIV_OPS_H
#define DIV_OPS_H

#include <bout/assert.hxx>
#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/coordinates.hxx>
#include <bout/field.hxx>
#include <bout/field3d.hxx>
#include <bout/fv_ops_impl.hxx>
#include <bout/output_bout_types.hxx>
#include <bout/region.hxx>
#include <bout/utils.hxx>
#include <bout/vector3d.hxx>

#include <cmath>

/*!
 * Diffusion in index space
 *
 * Similar to using Div_par_diffusion(SQ(mesh->dy)*mesh->g_22, f)
 *
 * @param[in] The field to be differentiated
 * @param[in] bndry_flux  Are fluxes through the boundary calculated?
 */
Field3D Div_par_diffusion_index(const Field3D& f, bool bndry_flux = true);

Field3D Div_n_bxGrad_f_B_XPPM(const Field3D& n, const Field3D& f, bool bndry_flux = true,
                              bool poloidal = false, bool positive = false);

/// This version has an extra coefficient 'g' that is linearly interpolated
/// onto cell faces
Field3D Div_n_g_bxGrad_f_B_XZ(const Field3D& n, const Field3D& g, const Field3D& f,
                              bool bndry_flux = true);

Field3D Div_Perp_Lap_FV_Index(const Field3D& a, const Field3D& f);

Field3D Div_Z_FV_Index(const Field3D& a, const Field3D& f);

// 4th-order flux conserving term, in index space
Field3D D4DX4_FV_Index(const Field3D& f, bool bndry_flux = false);
Field3D D4DZ4_Index(const Field3D& f);

// Div ( k * Grad(f) )
Field2D Laplace_FV(const Field2D& k, const Field2D& f);

/// Perpendicular diffusion including X and Y directions
/// Takes Div_a_Grad_perp from BOUT++ and adds flows
Field3D Div_a_Grad_perp_flows(const Field3D& a, const Field3D& f, Field3D& flux_xlow,
                              Field3D& flux_ylow);
/// Same but with upwinding
/// WARNING: Causes checkerboarding in neutral_mixed integrated test
Field3D Div_a_Grad_perp_upwind(const Field3D& a, const Field3D& f);
/// Same but with upwinding and flows
/// WARNING: Causes checkerboarding in neutral_mixed integrated test
Field3D Div_a_Grad_perp_upwind_flows(const Field3D& a, const Field3D& f,
                                     Field3D& flux_xlow, Field3D& flux_ylow);

/*!
 * Div ( a Grad_perp(f) ) -- ∇⊥ ( a ⋅ ∇⊥ f) -- Vorticity
 *
 * This version includes corrections for non-orthogonal meshes
 * in which the g12 and g13 components can be non-zero
 * i.e. X-Y, X-Z and Y-Z coordinates can all be non-orthogonal.
 */
Field3D Div_a_Grad_perp_nonorthog(const Field3D& a, const Field3D& x, Field3D& flux_xlow,
                                  Field3D& flux_ylow);

#endif //  DIV_OPS_H
