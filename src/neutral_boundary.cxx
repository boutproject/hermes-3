#include "bout/mesh.hxx"
#include <bout/constants.hxx>
using bout::globals::mesh;

#include "../include/neutral_boundary.hxx"

NeutralBoundary::NeutralBoundary(std::string name, Options& alloptions, Solver* solver)
    : name(name) {
  AUTO_TRACE();

  auto& options = alloptions[name];

  target_gamma_heat = options["target_gamma_heat"]
                   .doc("Neutral boundary heat transmission coefficient on target")
                   .withDefault(0.0);

  sol_gamma_heat = options["sol_gamma_heat"]
                   .doc("Neutral boundary heat transmission coefficient in SOL")
                   .withDefault(0.0);

  pfr_gamma_heat = options["pfr_gamma_heat"]
                   .doc("Neutral boundary heat transmission coefficient in PFR")
                   .withDefault(0.0);

  diagnose = options["diagnose"].doc("Save additional diagnostics?").withDefault<bool>(false);

  lower_y = options["neutral_boundary_lower_y"].doc("Boundary on lower y?").withDefault<bool>(true);
  upper_y = options["neutral_boundary_upper_y"].doc("Boundary on upper y?").withDefault<bool>(true);
  sol = options["neutral_boundary_sol"].doc("Boundary on SOL?").withDefault<bool>(false);
  pfr = options["neutral_boundary_pfr"].doc("Boundary on PFR?").withDefault<bool>(false);
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
  energy_source_diagnostic = 0;

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
        const BoutReal v_th = sqrt(tnsheath / AA);

        // Heat flux (> 0)
        const BoutReal q = target_gamma_heat * nnsheath * tnsheath * v_th;
        // Multiply by cell area to get power
        BoutReal flux = q * (coord->J[i] + coord->J[im])
                        / (sqrt(coord->g_22[i]) + sqrt(coord->g_22[im]));

        // Divide by volume of cell to get energy loss rate (> 0)
        BoutReal power = flux / (coord->dy[i] * coord->J[i]);

        // Subtract from cell next to boundary
        energy_source[i] -= power;
        energy_source_diagnostic[i] -= power;
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
        Nn[ip] = SQ(Nn[i]) / Nn[im];
        Pn[ip] = SQ(Pn[i]) / Pn[im];
        Tn[ip] = SQ(Tn[i]) / Tn[im];

        // No-flow boundary condition
        Vn[ip] = -Vn[i];
        NVn[ip] = -NVn[i];

        // Calculate midpoint values at wall
        const BoutReal nnsheath = 0.5 * (Nn[ip] + Nn[i]);
        const BoutReal tnsheath = 0.5 * (Tn[ip] + Tn[i]);

        // Thermal speed
        const BoutReal v_th = sqrt(tnsheath / AA);

        // Heat flux (> 0)
        const BoutReal q = target_gamma_heat * nnsheath * tnsheath * v_th;
        // Multiply by cell area to get power
        BoutReal flux = q * (coord->J[i] + coord->J[ip])
                        / (sqrt(coord->g_22[i]) + sqrt(coord->g_22[ip]));

        // Divide by volume of cell to get energy loss rate (> 0)
        BoutReal power = flux / (coord->dy[i] * coord->J[i]);

        // Subtract from cell next to boundary
        energy_source[i] -= power;
        energy_source_diagnostic[i] -= power;
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

          // Thermal speed in one direction only (?)
          const BoutReal v_th = sqrt(tnsheath / AA);

          // Heat flux (> 0)
          const BoutReal q = sol_gamma_heat * nnsheath * tnsheath * v_th;

          // Multiply by cell area to get power
          BoutReal flux = q * (coord->dy[i] + coord->dy[ig]) * (coord->dz[i] + coord->dz[ig])
                          / (sqrt(coord->g22[i]) + sqrt(coord->g22[ig]));

          // Divide by volume of cell to get energy loss rate (> 0)
          BoutReal power = flux / (coord->J[i] * coord->dx[i] * coord->dy[i] * coord->dz[i]);

          // Subtract from cell next to boundary
          energy_source[i] -= power;
          energy_source_diagnostic[i] -= power;

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

          // Thermal speed in one direction only (?)
          const BoutReal v_th = sqrt(tnsheath / AA);

          // Heat flux (> 0)
          const BoutReal q = pfr_gamma_heat * nnsheath * tnsheath * v_th;

          // Multiply by cell area to get power
          BoutReal flux = q * (coord->dy[i] + coord->dy[ig]) * (coord->dz[i] + coord->dz[ig])
                          / (sqrt(coord->g22[i]) + sqrt(coord->g22[ig]));

          // Divide by volume of cell to get energy loss rate (> 0)
          BoutReal power = flux / (coord->J[i] * coord->dx[i] * coord->dy[i] * coord->dz[i]);

          // Subtract from cell next to boundary
          energy_source[i] -= power;
          energy_source_diagnostic[i] -= power;

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
        set_with_attrs(state[{std::string("E") + name + std::string("_wall_refl")}], energy_source_diagnostic,
                        {{"time_dimension", "t"},
                        {"units", "W m^-3"},
                        {"conversion", Pnorm * Omega_ci},
                        {"standard_name", "energy source"},
                        {"long_name", std::string("Wall reflection energy source of ") + name},
                        {"source", "neutral_boundary"}});
  }
}

