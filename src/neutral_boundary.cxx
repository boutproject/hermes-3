#include "bout/mesh.hxx"
#include <bout/constants.hxx>
using bout::globals::mesh;

#include "../include/neutral_boundary.hxx"

NeutralBoundary::NeutralBoundary(std::string name, Options& alloptions, Solver* solver)
    : name(name) {
  AUTO_TRACE();

  auto& options = alloptions[name];
  const Options& units = alloptions["units"];
  Tnorm = units["eV"];

  diagnose = options["diagnose"].doc("Save additional diagnostics?").withDefault<bool>(false);
  lower_y = options["neutral_boundary_lower_y"].doc("Boundary on lower y?").withDefault<bool>(true);
  upper_y = options["neutral_boundary_upper_y"].doc("Boundary on upper y?").withDefault<bool>(true);
  sol = options["neutral_boundary_sol"].doc("Boundary on SOL?").withDefault<bool>(false);
  pfr = options["neutral_boundary_pfr"].doc("Boundary on PFR?").withDefault<bool>(false);

  target_energy_refl_factor =
        options["target_energy_refl_factor"]
            .doc("Fraction of energy retained by neutral particles after wall reflection at target")
            .withDefault<BoutReal>(0.75);

  sol_energy_refl_factor =
        options["sol_energy_refl_factor"]
            .doc("Fraction of energy retained by neutral particles after wall reflection at SOL")
            .withDefault<BoutReal>(0.75);

  pfr_energy_refl_factor =
        options["pfr_energy_refl_factor"]
            .doc("Fraction of energy retained by neutral particles after wall reflection at PFR")
            .withDefault<BoutReal>(0.75);

  target_fast_refl_fraction =
        options["target_fast_refl_fraction"]
            .doc("Fraction of neutrals that are undergoing fast reflection at the target")
            .withDefault<BoutReal>(0.8);

  sol_fast_refl_fraction =
        options["sol_fast_refl_fraction"]
            .doc("Fraction of neutrals that are undergoing fast reflection at the sol")
            .withDefault<BoutReal>(0.8);
  
  pfr_fast_refl_fraction =
        options["pfr_fast_refl_fraction"]
            .doc("Fraction of neutrals that are undergoing fast reflection at the pfr")
            .withDefault<BoutReal>(0.8);

}

void NeutralBoundary::transform(Options& state) {
  AUTO_TRACE();
  auto& species = state["species"][name];
  const BoutReal AA = get<BoutReal>(species["AA"]);

  Field3D Nn = toFieldAligned(GET_NOBOUNDARY(Field3D, species["density"]));
  Field3D Pn = toFieldAligned(GET_NOBOUNDARY(Field3D, species["pressure"]));
  Field3D Tn = toFieldAligned(GET_NOBOUNDARY(Field3D, species["temperature"]));

  Field3D Vn = IS_SET_NOBOUNDARY(species["velocity"])
                   ? toFieldAligned(getNoBoundary<Field3D>(species["velocity"]))
                   : zeroFrom(Nn);

  Field3D NVn = IS_SET_NOBOUNDARY(species["momentum"])
                    ? toFieldAligned(getNoBoundary<Field3D>(species["momentum"]))
                    : zeroFrom(Nn);

  // Get the energy source, or create if not set
  Field3D energy_source =
      species.isSet("energy_source")
          ? toFieldAligned(getNonFinal<Field3D>(species["energy_source"]))
          : zeroFrom(Nn);

  Coordinates* coord = mesh->getCoordinates();
  target_energy_source = 0;
  wall_energy_source = 0;

  // Targets
  if (lower_y) {
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        auto i = indexAt(Nn, r.ind, mesh->ystart, jz);
        auto im = i.ym();
        auto ip = i.yp();

        // Free boundary condition on log(Nn), log(Pn)
        Nn[im] = SQ(Nn[i]) / Nn[ip];
        Pn[im] = SQ(Pn[i]) / Pn[ip];
        Tn[im] = SQ(Tn[i]) / Tn[ip];

        // No-flow boundary condition
        Vn[im] = -Vn[i];
        NVn[im] = -NVn[i];

        // Calculate midpoint values at wall
        const BoutReal nnsheath = 0.5 * (Nn[im] + Nn[i]);
        const BoutReal tnsheath = 0.5 * (Tn[im] + Tn[i]);

        // Thermal speed
        const BoutReal v_th = 0.25 * sqrt( 8*tnsheath / (PI*AA) );   // Stangeby p.69 eqns. 2.21, 2.24

        // Calculate effective gamma from particle and energy reflection coefficients
        BoutReal target_gamma_heat = 1 - target_energy_refl_factor * target_fast_refl_fraction 
                                  -(1-target_fast_refl_fraction) * (3/Tnorm) / (2*tnsheath);  // D. Power thesis 2023

        // Heat flux (> 0)
        const BoutReal q = target_gamma_heat * nnsheath * tnsheath * v_th;
        // Multiply by cell area to get power
        BoutReal flux = q * (coord->J[i] + coord->J[im])
                        / (sqrt(coord->g_22[i]) + sqrt(coord->g_22[im]));

        // Divide by volume of cell to get energy loss rate (> 0)
        BoutReal power = flux / (coord->dy[i] * coord->J[i]);

        // Subtract from cell next to boundary
        energy_source[i] -= power;
        target_energy_source[i] -= power;
      }
    }
  }

  if (upper_y) {
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        auto i = indexAt(Nn, r.ind, mesh->yend, jz);
        auto im = i.ym();
        auto ip = i.yp();

        // Free boundary condition on log(Nn), log(Pn)
        // This is problematic when Nn, Pn or Tn are zero
        // Nn[ip] = SQ(Nn[i]) / Nn[im];
        // Pn[ip] = SQ(Pn[i]) / Pn[im];
        // Tn[ip] = SQ(Tn[i]) / Tn[im];

        // Neumann boundary condition: do not extrapolate, but 
        // assume the target value is same as the final cell centre.
        // Shouldn't affect results much and more resilient to positivity issues
        Nn[ip] = Nn[i];
        Pn[ip] = Pn[i];
        Tn[ip] = Tn[i];

        // No-flow boundary condition
        Vn[ip] = -Vn[i];
        NVn[ip] = -NVn[i];

        // Calculate midpoint values at wall
        const BoutReal nnsheath = 0.5 * (Nn[ip] + Nn[i]);
        const BoutReal tnsheath = 0.5 * (Tn[ip] + Tn[i]);

        // Thermal speed
        const BoutReal v_th = 0.25 * sqrt( 8*tnsheath / (PI*AA) );   // Stangeby p.69 eqns. 2.21, 2.24

        // Approach adapted from D. Power thesis 2023
        BoutReal T_FC = 3 / Tnorm; // Franck-Condon temp (hardcoded for now)

        // Outgoing neutral heat flux [W/m^2]
        BoutReal q = 
                      (1 - target_energy_refl_factor * target_fast_refl_fraction ) * nnsheath * tnsheath * v_th  // unreflected energy
                    - (1 - target_fast_refl_fraction) * T_FC * 0.5 * nnsheath * v_th;  // energy returning as FC

        // Multiply by cell area to get neutral heat flow [W]
        BoutReal flow = q * (coord->J[i] + coord->J[ip])
                        * ((coord->dx[i] * coord->dz[i]) / (sqrt(coord->g_22[i]) + sqrt(coord->g_22[ip])));

        // Divide by cell volume to get source [W/m^3]
        BoutReal cooling_source = flow / (coord->dx[i] * coord->dy[i] * coord->dz[i] * coord->J[i]);

        // Subtract from cell next to boundary
        energy_source[i] -= cooling_source;
        target_energy_source[i] -= cooling_source;

      }
    }
  }

    // SOL edge
  if (sol) {
    if(mesh->lastX()){  // Only do this for the processor which has the edge region
      for(int iy=0; iy < mesh->LocalNy ; iy++){
        for(int iz=0; iz < mesh->LocalNz; iz++){

          auto i = indexAt(Nn, mesh->xend, iy, iz);  // Final domain cell
          auto ig = indexAt(Nn, mesh->xend+1, iy, iz);  // Guard cell
          
          // Calculate midpoint values at wall
          const BoutReal nnsheath = 0.5 * (Nn[ig] + Nn[i]);
          const BoutReal tnsheath = 0.5 * (Tn[ig] + Tn[i]);

          // Thermal speed of static Maxwellian in one direction
          const BoutReal v_th = 0.25 * sqrt( 8*tnsheath / (PI*AA) );   // Stangeby p.69 eqns. 2.21, 2.24

          // Calculate effective gamma from particle and energy reflection coefficients
          BoutReal sol_gamma_heat = 1 - sol_energy_refl_factor * sol_fast_refl_fraction 
                                    -(1-sol_fast_refl_fraction) * (3/Tnorm) / (2*tnsheath);  // D. Power thesis 2023

          // Heat flux (> 0)
          const BoutReal q = sol_gamma_heat * nnsheath * tnsheath * v_th;

          // Multiply by radial cell area to get power
          // Expanded form of the calculation for clarity

          // Converts dy to poloidal length: dl = dy * sqrt(g22) = dy * h_theta
          BoutReal dpolsheath = 0.5*(coord->dy[i] + coord->dy[ig]) *  1/( 0.5*(sqrt(coord->g22[i]) + sqrt(coord->g22[ig])) );

          // Converts dz to toroidal length:  = dz*sqrt(g_33) = dz * R = 2piR
          BoutReal dtorsheath = 0.5*(coord->dz[i] + coord->dz[ig]) * 0.5*(sqrt(coord->g_33[i]) + sqrt(coord->g_33[ig]));

          BoutReal dasheath = dpolsheath * dtorsheath;  // [m^2]
          BoutReal flux = q * dasheath;  // [W]

          // Divide by volume of cell to get energy loss rate (> 0)
          BoutReal power = flux / (coord->J[i] * coord->dx[i] * coord->dy[i] * coord->dz[i]);   // [W m^-3]

          // Subtract from cell next to boundary
          energy_source[i] -= power;
          wall_energy_source[i] -= power;

        }
      }
    }
  }

  // PFR edge
  if (pfr) {
    if ((mesh->firstX()) and (!mesh->periodicY(mesh->xstart))) {  // do loop if inner edge and not periodic (i.e. PFR)
      for(int iy=0; iy < mesh->LocalNy ; iy++){
        for(int iz=0; iz < mesh->LocalNz; iz++){

          auto i = indexAt(Nn, mesh->xstart, iy, iz);  // Final domain cell
          auto ig = indexAt(Nn, mesh->xstart-1, iy, iz);  // Guard cell
          
          // Calculate midpoint values at wall
          const BoutReal nnsheath = 0.5 * (Nn[ig] + Nn[i]);
          const BoutReal tnsheath = 0.5 * (Tn[ig] + Tn[i]);

          // Thermal speed of static Maxwellian in one direction
          const BoutReal v_th = 0.25 * sqrt( 8*tnsheath / (PI*AA) );   // Stangeby p.69 eqns. 2.21, 2.24

          
          // Calculate effective gamma from particle and energy reflection coefficients
          BoutReal pfr_gamma_heat = 1 - pfr_energy_refl_factor * pfr_fast_refl_fraction 
                                    -(1-pfr_fast_refl_fraction) * (3/Tnorm) / (2*tnsheath);  // D. Power thesis 2023

          // Heat flux (> 0)
          const BoutReal q = pfr_gamma_heat * nnsheath * tnsheath * v_th;

          // Multiply by radial cell area to get power
          // Expanded form of the calculation for clarity

          // Converts dy to poloidal length: dl = dy * sqrt(g22) = dy * h_theta
          BoutReal dpolsheath = 0.5*(coord->dy[i] + coord->dy[ig]) *  1/( 0.5*(sqrt(coord->g22[i]) + sqrt(coord->g22[ig])) );

          // Converts dz to toroidal length:  = dz*sqrt(g_33) = dz * R = 2piR
          BoutReal dtorsheath = 0.5*(coord->dz[i] + coord->dz[ig]) * 0.5*(sqrt(coord->g_33[i]) + sqrt(coord->g_33[ig]));

          BoutReal dasheath = dpolsheath * dtorsheath;  // [m^2]
          BoutReal flux = q * dasheath;  // [W]

          // Divide by volume of cell to get energy loss rate (> 0)
          BoutReal power = flux / (coord->J[i] * coord->dx[i] * coord->dy[i] * coord->dz[i]);   // [W m^-3]

          // Subtract from cell next to boundary
          energy_source[i] -= power;
          wall_energy_source[i] -= power;

        }
      }
    }
  }

  // Set density, pressure and temperature, now with boundary conditions
  setBoundary(species["density"], fromFieldAligned(Nn));
  setBoundary(species["temperature"], fromFieldAligned(Tn));
  setBoundary(species["pressure"], fromFieldAligned(Pn));
  if (IS_SET_NOBOUNDARY(species["velocity"])) {
    setBoundary(species["velocity"], fromFieldAligned(Vn));
  }
  if (IS_SET_NOBOUNDARY(species["momentum"])) {
    setBoundary(species["momentum"], fromFieldAligned(NVn));
  }

  // Set energy source (negative in cell next to sheath)
  // Note: energy_source includes any sources previously set in other components
  set(species["energy_source"], fromFieldAligned(energy_source));
}

void NeutralBoundary::outputVars(Options& state) {
  
  AUTO_TRACE();
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

  if (diagnose) {

      AUTO_TRACE();

      // Save particle and energy source for the species created during recycling

      // Target recycling

      if ((sol) or (pfr)) {
        set_with_attrs(state[{std::string("E") + name + std::string("_wall_refl")}], wall_energy_source,
                        {{"time_dimension", "t"},
                        {"units", "W m^-3"},
                        {"conversion", Pnorm * Omega_ci},
                        {"standard_name", "energy source"},
                        {"long_name", std::string("Wall reflection energy source of ") + name},
                        {"source", "neutral_boundary"}});
      }

      set_with_attrs(state[{std::string("E") + name + std::string("_target_refl")}], target_energy_source,
                      {{"time_dimension", "t"},
                      {"units", "W m^-3"},
                      {"conversion", Pnorm * Omega_ci},
                      {"standard_name", "energy source"},
                      {"long_name", std::string("Wall reflection energy source of ") + name},
                      {"source", "neutral_boundary"}});
  }
}

