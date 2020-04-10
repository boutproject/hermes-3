/*

    Copyright B.Dudson, J.Leddy, University of York, 2016-2020
              email: benjamin.dudson@york.ac.uk

    This file is part of Hermes-3 (Hot ion, multifluid)

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

#include "hermes-3.hxx"

#include <bout/constants.hxx>

int Hermes::init(bool restarting) {
  auto &options = Options::root()["hermes"];

  // For forward compatibility with multiple meshes
  MeshMap meshes = {{"mesh", bout::globals::mesh}};

  // Choose normalisations
  Tnorm = options["Tnorm"].doc("Reference temperature [eV]").withDefault(100.);
  Nnorm = options["Nnorm"].doc("Reference density [m^-3]").withDefault(1e19);
  Bnorm = options["Bnorm"].doc("Reference magnetic field [T]").withDefault(1.0);

  Cs0 = sqrt(SI::qe * Tnorm / SI::Mp); // Reference sound speed [m/s]
  Omega_ci = SI::qe * Bnorm / SI::Mp;  // Ion cyclotron frequency [1/s]
  rho_s0 = Cs0 / Omega_ci;             // Length scale [m]

  SAVE_ONCE(Tnorm, Nnorm, Bnorm); // Save normalisations
  SAVE_ONCE(Cs0, Omega_ci, rho_s0);

  // Put normalisation quantities into an Options to use later
  units["inv_meters_cubed"] = Nnorm;
  units["eV"] = Tnorm;
  units["Tesla"] = Bnorm;
  units["seconds"] = 1./Omega_ci;
  units["meters"] = rho_s0;

  // Put into the options tree, so quantities can be normalised
  // when creating components
  options["units"] = units;

  // Tell the components if they are restarting
  options["restarting"] = restarting;
  
  // Create the components
  scheduler = ComponentScheduler::create(options, meshes);
  
  return 0;
}

int Hermes::rhs(BoutReal time) {
  Options state;

  state["time"] = time;
  state["units"] = units; 

  // Call all the components
  scheduler->transform(state);

  return 0;
}

/*!
 * Preconditioner. Solves the heat conduction
 *
 * @param[in] t  The simulation time
 * @param[in] gamma   Factor in front of the Jacobian in (I - gamma*J). Related
 * to timestep
 * @param[in] delta   Not used here
 */
int Hermes::precon(BoutReal t, BoutReal gamma, BoutReal delta) { return 0; }

// Standard main() function
BOUTMAIN(Hermes);
