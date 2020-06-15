#include "../include/sheath_boundary.hxx"

#include "bout/constants.hxx"
#include "bout/mesh.hxx"
using bout::globals::mesh;

BoutReal clip(BoutReal value, BoutReal min, BoutReal max) {
  if (value < min)
    return min;
  if (value > max)
    return max;
  return value;
}

BoutReal floor(BoutReal value, BoutReal min) {
  if (value < min)
    return min;
  return value;
}

SheathBoundary::SheathBoundary(std::string name, Options &alloptions, Solver *) {
  const Options& units = alloptions["units"];
  
  Options &options = alloptions[name];

  Ge = options["secondary_electron_coef"]
           .doc("Effective secondary electron emission coefficient")
           .withDefault(0.0);

  sin_alpha = options["sin_alpha"]
                  .doc("Sin of the angle between magnetic field line and wall surface. "
                       "Should be between 0 and 1")
                  .withDefault(1.0);

  if ((sin_alpha < 0.0) or (sin_alpha > 1.0)) {
    throw BoutException("Range of sin_alpha must be between 0 and 1");
  }

  lower_y = options["lower_y"].doc("Boundary on lower y?").withDefault<bool>(true);
  upper_y = options["upper_y"].doc("Boundary on upper y?").withDefault<bool>(true);
}

void SheathBoundary::transform(Options &state) {
  Options& allspecies = state["species"];
  Options& electrons = allspecies["e"];
  
  // Need electron properties
  const Field3D Ne = get<Field3D>(electrons["density"]);
  const Field3D Te = get<Field3D>(electrons["temperature"]);

  // Electrostatic potential
  // If phi is set, use free boundary condition
  // If phi not set, calculate assuming zero current
  Field3D phi;
  if (state.isSection("fields") and state["fields"].isSet("phi")) {
    phi = get<Field3D>(state["fields"]["phi"]);
  } else {
    // Calculate potential phi assuming zero current
    // Note: This is equation (22) in Tskhakaya 2005, with I = 0
    
    // Need to sum  s_i Z_i C_i over all ion species
    //
    // To avoid looking up species for every grid point, this
    // loops over the boundaries once per species.
    Field3D ion_sum = 0.0;

    // Iterate through charged ion species
    for (auto& kv : allspecies.getChildren()) {
      const Options& species = kv.second;
      
      if ((kv.first == "e") or !species["charge"].isSet("charge")) {
        continue; // Skip electrons and non-charged ions
      }
      
      const Field3D Ni = get<Field3D>(species["density"]);
      const Field3D Ti = get<Field3D>(species["temperature"]);
      const BoutReal Mi = get<BoutReal>(species["mass"]);
      const BoutReal Zi = get<BoutReal>(species["charge"]);

      const BoutReal adiabatic = 3./2; // Ratio of specific heats
      
      if (lower_y) {
        // Sum values, put result in mesh->ystart
        
        for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            // Free boundary extrapolate ion concentration
            const BoutReal s_i =
                clip(0.5
                         * (3. * Ni(r.ind, mesh->ystart, jz) / Ne(r.ind, mesh->ystart, jz)
                            - Ni(r.ind, mesh->ystart + 1, jz)
                                  / Ne(r.ind, mesh->ystart + 1, jz)),
                     0.0, 1.0); // Limit range to [0,1]

            BoutReal te = Te(r.ind, mesh->ystart, jz);
            BoutReal ti = Ti(r.ind, mesh->ystart, jz);
            
            // Equation (9) in Tskhakaya 2005
            BoutReal C_i_sq = clip(
                (adiabatic * ti
                 + Zi * s_i * te
                       * (Ne(r.ind, mesh->ystart + 1, jz) - Ne(r.ind, mesh->ystart, jz))
                       / (Ni(r.ind, mesh->ystart + 1, jz) - Ni(r.ind, mesh->ystart, jz)))
                    / Mi,
                0, 100); // Limit for e.g. Ni zero gradient
            
            ion_sum(r.ind, mesh->ystart, jz) += s_i * Zi * sqrt(C_i_sq);
          }
        }
      }

      if (upper_y) {
        // Sum values, put results in mesh->yend
        
        for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            const BoutReal s_i = clip(
                0.5
                    * (3. * Ni(r.ind, mesh->yend, jz) / Ne(r.ind, mesh->yend, jz)
                       - Ni(r.ind, mesh->yend - 1, jz) / Ne(r.ind, mesh->yend - 1, jz)),
                0.0, 1.0);

            BoutReal te = Te(r.ind, mesh->yend, jz);
            BoutReal ti = Ti(r.ind, mesh->yend, jz);

            BoutReal C_i_sq = clip(
                (adiabatic * ti
                 + Zi * s_i * te
                       * (Ne(r.ind, mesh->yend - 1, jz) - Ne(r.ind, mesh->yend, jz))
                       / (Ni(r.ind, mesh->yend - 1, jz) - Ni(r.ind, mesh->yend, jz)))
                    / Mi,
                0, 100); // Limit for e.g. Ni zero gradient
            
            ion_sum(r.ind, mesh->yend, jz) += s_i * Zi * sqrt(C_i_sq);
          }
        }
      }
    }

    phi.allocate();
    
    // ion_sum now contains  sum  s_i Z_i C_i over all ion species
    // at mesh->ystart and mesh->yend indices
    if (lower_y) {
      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          Ind3D i {r.ind, mesh->ystart, jz};
          phi[i] =
              Te[i]
              * log(sqrt(Te[i] * SI::Mp / (SI::Me * TWOPI)) * (1. - Ge) / ion_sum[i]);

          phi[i.ym()] = phi[i]; // Constant into sheath
        }
      }
    }

    if (upper_y) {
      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          Ind3D i {r.ind, mesh->yend, jz};
          phi[i] =
              Te[i]
              * log(sqrt(Te[i] * SI::Mp / (SI::Me * TWOPI)) * (1. - Ge) / ion_sum[i]);
          phi[i.yp()] = phi[i];
        }
      }
    }

    // Set the potential at the wall
    set(state["fields"]["phi"], phi);
  }

  // Iterate through all ions
  for (auto& kv : allspecies.getChildren()) {
    if (kv.first == "e") {
      continue; // Skip electrons
    }

    Options& species = allspecies[kv.first]; // Note: Need non-const

    // Ion charge
    const BoutReal Zi =
        species["charge"].isSet("charge") ? get<BoutReal>(species["charge"]) : 0.0;
    
    if (Zi == 0.0)
      continue; // Neutral -> skip
    
    const Field3D Ni = get<Field3D>(species["density"]);
    const Field3D Ti = get<Field3D>(species["temperature"]);
    const BoutReal Mi = get<BoutReal>(species["mass"]);
    
    if (lower_y) {
      
      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          // Free gradient electron density
          BoutReal nesheath = 0.5 * (3. * Ne(r.ind, mesh->ystart, jz) -
                                     Ne(r.ind, mesh->ystart + 1, jz));
          if (nesheath < 0.0)
            nesheath = 0.0;
          
          // Electron temperature at the sheath entrance
          BoutReal tesheath = floor(Te(r.ind, mesh->ystart, jz), 0.0);

          BoutReal phisheath;
          if (phi.isAllocated()) {
            // Free boundary potential
            phisheath =
                0.5
                * (3. * phi(r.ind, mesh->ystart, jz) - phi(r.ind, mesh->ystart + 1, jz));
            
            // Electron saturation at phi = 0
            if (phisheath < 0.0) {
              phisheath = 0.0;
            }
          } else {
            
          }
          
          // Electron sheath heat transmission
          
          // Ion sheath heat transmission coefficient
          
          
        }
      }
    }
  }
}
