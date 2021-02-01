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

const Field3D Div_a_Laplace_perp_upwind(const Field3D& a, const Field3D& f);

#endif //  __DIV_OPS_H__
