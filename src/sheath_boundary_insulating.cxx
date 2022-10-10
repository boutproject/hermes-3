/// Implements an insulating sheath, so J = 0 but ions still flow
/// to the wall. Potential (if set) is linearly extrapolated into the boundary.

#include "../include/sheath_boundary_insulating.hxx"

#include <bout/output_bout_types.hxx>

#include "bout/constants.hxx"
#include "bout/mesh.hxx"
using bout::globals::mesh;

namespace {
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

Ind3D indexAt(const Field3D& f, int x, int y, int z) {
  int ny = f.getNy();
  int nz = f.getNz();
  return Ind3D{(x * ny + y) * nz + z, ny, nz};
}

/// Limited free gradient of log of a quantity
/// This ensures that the guard cell values remain positive
/// while also ensuring that the quantity never increases
///
///  fm  fc | fp
///         ^ boundary
///
/// exp( 2*log(fc) - log(fm) )
///
BoutReal limitFree(BoutReal fm, BoutReal fc) {
  if (fm < fc) {
    return fc; // Neumann rather than increasing into boundary
  }
  if (fm < 1e-10) {
    return fc; // Low / no density condition
  }
  BoutReal fp = SQ(fc) / fm;
#if CHECKLEVEL >= 2
  if (!std::isfinite(fp)) {
    throw BoutException("SheathBoundary limitFree: {}, {} -> {}", fm, fc, fp);
  }
#endif

  return fp;
}

}

SheathBoundaryInsulating::SheathBoundaryInsulating(std::string name, Options &alloptions, Solver *) {
  AUTO_TRACE();

  Options &options = alloptions[name];

  Ge = options["secondary_electron_coef"]
           .doc("Effective secondary electron emission coefficient")
           .withDefault(0.0);

  if ((Ge < 0.0) or (Ge > 1.0)) {
    throw BoutException("Secondary electron emission must be between 0 and 1 ({:e})", Ge);
  }
  
  sin_alpha = options["sin_alpha"]
                  .doc("Sin of the angle between magnetic field line and wall surface. "
                       "Should be between 0 and 1")
                  .withDefault(1.0);

  if ((sin_alpha < 0.0) or (sin_alpha > 1.0)) {
    throw BoutException("Range of sin_alpha must be between 0 and 1");
  }

  lower_y = options["lower_y"].doc("Boundary on lower y?").withDefault<bool>(true);
  upper_y = options["upper_y"].doc("Boundary on upper y?").withDefault<bool>(true);

  gamma_e = options["gamma_e"]
    .doc("Electron sheath heat transmission coefficient")
    .withDefault(3.5);
}

void SheathBoundaryInsulating::transform(Options &state) {
  AUTO_TRACE();

  Options& allspecies = state["species"];
  Options& electrons = allspecies["e"];

  // Need electron properties
  // Not const because boundary conditions will be set
  Field3D Ne = toFieldAligned(floor(GET_NOBOUNDARY(Field3D, electrons["density"]), 0.0));
  Field3D Te = toFieldAligned(GET_NOBOUNDARY(Field3D, electrons["temperature"]));
  Field3D Pe = IS_SET_NOBOUNDARY(electrons["pressure"])
    ? toFieldAligned(getNoBoundary<Field3D>(electrons["pressure"]))
    : Te * Ne;

  // Ratio of specific heats
  const BoutReal electron_adiabatic =
      IS_SET(electrons["adiabatic"]) ? get<BoutReal>(electrons["adiabatic"]) : 5. / 3;

  // Mass, normalised to proton mass
  const BoutReal Me =
      IS_SET(electrons["AA"]) ? get<BoutReal>(electrons["AA"]) : SI::Me / SI::Mp;

  // This is for applying boundary conditions
  Field3D Ve = IS_SET_NOBOUNDARY(electrons["velocity"])
    ? toFieldAligned(getNoBoundary<Field3D>(electrons["velocity"]))
    : zeroFrom(Ne);

  Field3D NVe = IS_SET_NOBOUNDARY(electrons["momentum"])
    ? toFieldAligned(getNoBoundary<Field3D>(electrons["momentum"]))
    : zeroFrom(Ne);

  Coordinates *coord = mesh->getCoordinates();

  //////////////////////////////////////////////////////////////////
  // Electrostatic potential
  // If phi is set, set free boundary condition

  Field3D phi;
  if (IS_SET_NOBOUNDARY(state["fields"]["phi"])) {
    phi = toFieldAligned(getNoBoundary<Field3D>(state["fields"]["phi"]));

    // Free boundary potential linearly extrapolated
    if (lower_y) {
      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          auto i = indexAt(Ne, r.ind, mesh->ystart, jz);
          phi[i.ym()] = 2 * phi[i] - phi[i.yp()];
        }
      }
    }
    if (upper_y) {
      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          auto i = indexAt(Ne, r.ind, mesh->yend, jz);
          phi[i.yp()] = 2 * phi[i] - phi[i.ym()];
        }
      }
    }
    // Set the potential, including boundary conditions
    phi = fromFieldAligned(phi);
    setBoundary(state["fields"]["phi"], phi);
  }

  //////////////////////////////////////////////////////////////////
  // Electrons

  if (lower_y) {
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        auto i = indexAt(Ne, r.ind, mesh->ystart, jz);
        auto ip = i.yp();
        auto im = i.ym();

        // Free gradient of log electron density and temperature
        // Limited so that the values don't increase into the sheath
        // This ensures that the guard cell values remain positive
        // exp( 2*log(N[i]) - log(N[ip]) )

        Ne[im] = limitFree(Ne[ip], Ne[i]);
        Te[im] = limitFree(Te[ip], Te[i]);
        Pe[im] = limitFree(Pe[ip], Pe[i]);

        // Set zero flow through boundary
        // This will be modified when iterating over the ions

        Ve[im] =  - Ve[i];
        NVe[im] = - NVe[i];
      }
    }
  }
  if (upper_y) {
    // This is essentially the same as at the lower y boundary
    // except ystart -> yend, ip <-> im
    // 
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        auto i = indexAt(Ne, r.ind, mesh->yend, jz);
        auto ip = i.yp();
        auto im = i.ym();

        Ne[ip] = limitFree(Ne[im], Ne[i]);
        Te[ip] = limitFree(Te[im], Te[i]);
        Pe[ip] = limitFree(Pe[im], Pe[i]);

        Ve[ip] = - Ve[i];
        NVe[ip] = - NVe[i];
      }
    }
  }

  // Set electron density and temperature, now with boundary conditions
  setBoundary(electrons["density"], fromFieldAligned(Ne));
  setBoundary(electrons["temperature"], fromFieldAligned(Te));
  setBoundary(electrons["pressure"], fromFieldAligned(Pe));

  //////////////////////////////////////////////////////////////////
  // Iterate through all ions
  // Sum the ion currents into the wall to calculate the electron flow

  for (auto& kv : allspecies.getChildren()) {
    if (kv.first == "e") {
      continue; // Skip electrons
    }

    Options& species = allspecies[kv.first]; // Note: Need non-const

    // Ion charge
    const BoutReal Zi =
        IS_SET(species["charge"]) ? get<BoutReal>(species["charge"]) : 0.0;

    if (Zi == 0.0) {
      continue; // Neutral -> skip
    }

    // Characteristics of this species
    const BoutReal Mi = get<BoutReal>(species["AA"]);

    const BoutReal adiabatic = IS_SET(species["adiabatic"])
                                   ? get<BoutReal>(species["adiabatic"])
                                   : 5. / 3; // Ratio of specific heats (ideal gas)

    // Density and temperature boundary conditions will be imposed (free)
    Field3D Ni = toFieldAligned(floor(getNoBoundary<Field3D>(species["density"]), 0.0));
    Field3D Ti = toFieldAligned(getNoBoundary<Field3D>(species["temperature"]));
    Field3D Pi = species.isSet("pressure")
      ? toFieldAligned(getNoBoundary<Field3D>(species["pressure"]))
      : Ni * Ti;

    // Get the velocity and momentum
    // These will be modified at the boundaries
    // and then put back into the state
    Field3D Vi = species.isSet("velocity")
      ? toFieldAligned(getNoBoundary<Field3D>(species["velocity"]))
      : zeroFrom(Ni);
    Field3D NVi = species.isSet("momentum")
      ? toFieldAligned(getNoBoundary<Field3D>(species["momentum"]))
      : Mi * Ni * Vi;

    // Energy source will be modified in the domain
    Field3D energy_source = species.isSet("energy_source")
      ? toFieldAligned(getNonFinal<Field3D>(species["energy_source"]))
      : zeroFrom(Ni);

    if (lower_y) {
      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          auto i = indexAt(Ne, r.ind, mesh->ystart, jz);
          auto ip = i.yp();
          auto im = i.ym();

          // Free gradient of log electron density and temperature
          // This ensures that the guard cell values remain positive
          // exp( 2*log(N[i]) - log(N[ip]) )

          Ni[im] = limitFree(Ni[ip], Ni[i]);
          Ti[im] = limitFree(Ti[ip], Ti[i]);
          Pi[im] = limitFree(Pi[ip], Pi[i]);

          // Calculate sheath values at half-way points (cell edge)
          const BoutReal nesheath = 0.5 * (Ne[im] + Ne[i]);
          const BoutReal nisheath = 0.5 * (Ni[im] + Ni[i]);
          const BoutReal tesheath = floor(0.5 * (Te[im] + Te[i]), 1e-5);  // electron temperature
          const BoutReal tisheath = floor(0.5 * (Ti[im] + Ti[i]), 1e-5);  // ion temperature

          // Ion sheath heat transmission coefficient
          // Equation (22) in Tskhakaya 2005
          // with 
          //
          // 1 / (1 + ∂_{ln n_e} ln s_i = s_i ∂_z n_e / ∂_z n_i
          // (from comparing C_i^2 in eq. 9 with eq. 20
          //
          // 
          BoutReal s_i = clip(nisheath / floor(nesheath, 1e-10), 0, 1); // Concentration
          BoutReal grad_ne = Ne[i] - nesheath;
          BoutReal grad_ni = Ni[i] - nisheath;

          if (fabs(grad_ni) < 1e-3) {
            grad_ni = grad_ne = 1e-3; // Remove kinetic correction term
          }

          // Ion speed into sheath
          // Equation (9) in Tskhakaya 2005
          //
          BoutReal C_i_sq =
              clip((adiabatic * tisheath + Zi * s_i * tesheath * grad_ne / grad_ni) / Mi,
                   0, 100); // Limit for e.g. Ni zero gradient

          // Ion sheath heat transmission coefficient
          const BoutReal gamma_i = 2.5 + 0.5 * Mi * C_i_sq / tisheath;

          const BoutReal visheath = - sqrt(C_i_sq); // Negative -> into sheath

          // Set boundary conditions on flows
          Vi[im] = 2. * visheath - Vi[i];
          NVi[im] = 2. * Mi * nisheath * visheath - NVi[i];

          // Add electron flow to balance current
          Ve[im] += 2. * visheath * Zi;
          NVe[im] += 2. * Me * nisheath * visheath; 

          // Take into account the flow of energy due to fluid flow
          // This is additional energy flux through the sheath
          // Note: Here this is negative because visheath < 0
          BoutReal q =
              ((gamma_i - 1 - 1 / (adiabatic - 1)) * tisheath - 0.5 * Mi * C_i_sq)
              * nisheath * visheath;
          if (q > 0.0) {
            q = 0.0;
          }

          // Multiply by cell area to get power
          BoutReal flux = q * (coord->J[i] + coord->J[im])
                          / (sqrt(coord->g_22[i]) + sqrt(coord->g_22[im]));

          // Divide by volume of cell to get energy loss rate (< 0)
          BoutReal power = flux / (coord->dy[i] * coord->J[i]);
	  ASSERT1(std::isfinite(power));
          ASSERT2(power <= 0.0);

          energy_source[i] += power;
        }
      }
    }
    if (upper_y) {
      // Note: This is essentially the same as the lower boundary,
      // but with directions reversed e.g. ystart -> yend, ip <-> im
      // 
      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          auto i = indexAt(Ne, r.ind, mesh->yend, jz);
          auto ip = i.yp();
          auto im = i.ym();

          // Free gradient of log electron density and temperature
          // This ensures that the guard cell values remain positive
          // exp( 2*log(N[i]) - log(N[ip]) )

          Ni[ip] = limitFree(Ni[im], Ni[i]);
          Ti[ip] = limitFree(Ti[im], Ti[i]);
          Pi[ip] = limitFree(Pi[im], Pi[i]);

          // Calculate sheath values at half-way points (cell edge)
          const BoutReal nesheath = 0.5 * (Ne[ip] + Ne[i]);
          const BoutReal nisheath = 0.5 * (Ni[ip] + Ni[i]);
          const BoutReal tesheath = floor(0.5 * (Te[ip] + Te[i]), 1e-5);  // electron temperature
          const BoutReal tisheath = floor(0.5 * (Ti[ip] + Ti[i]), 1e-5);  // ion temperature

          // Ion sheath heat transmission coefficient
          //
          // 1 / (1 + ∂_{ln n_e} ln s_i = s_i * ∂n_e / (s_i * ∂n_e + ∂ n_i) 
          BoutReal s_i = (nesheath > 1e-5) ? nisheath / nesheath : 0.0; // Concentration
          BoutReal grad_ne = Ne[i] - nesheath;
          BoutReal grad_ni = Ni[i] - nisheath;

          if (fabs(grad_ni) < 1e-3) {
            grad_ni = grad_ne = 1e-3; // Remove kinetic correction term
          }

          // Ion speed into sheath
          // Equation (9) in Tskhakaya 2005
          //
          BoutReal C_i_sq =
              clip((adiabatic * tisheath + Zi * s_i * tesheath * grad_ne / grad_ni) / Mi,
                   0, 100); // Limit for e.g. Ni zero gradient

          const BoutReal gamma_i = 2.5 + 0.5 * Mi * C_i_sq / tisheath; // + Δγ 

          const BoutReal visheath = sqrt(C_i_sq); // Positive -> into sheath

          // Set boundary conditions on flows
          Vi[ip] = 2. * visheath - Vi[i];
          NVi[ip] = 2. * Mi * nisheath * visheath - NVi[i];

          // Add electron flow to balance current
          Ve[ip] += 2. * visheath * Zi;
          NVe[ip] += 2. * Me * nisheath * visheath; 

          // Take into account the flow of energy due to fluid flow
          // This is additional energy flux through the sheath
          // Note: Here this is positive because visheath > 0
          BoutReal q =
              ((gamma_i - 1 - 1 / (adiabatic - 1)) * tisheath - 0.5 * C_i_sq * Mi)
              * nisheath * visheath;

          if (q < 0.0) {
            q = 0.0;
          }

          // Multiply by cell area to get power
          BoutReal flux = q * (coord->J[i] + coord->J[ip])
                          / (sqrt(coord->g_22[i]) + sqrt(coord->g_22[ip]));

          // Divide by volume of cell to get energy loss rate (> 0)
          BoutReal power = flux / (coord->dy[i] * coord->J[i]);
          ASSERT1(std::isfinite(power));
          ASSERT2(power >= 0.0);

          energy_source[i] -= power; // Note: Sign negative because power > 0
        }
      }
    }

    // Finished boundary conditions for this species
    // Put the modified fields back into the state.
    setBoundary(species["density"], fromFieldAligned(Ni));
    setBoundary(species["temperature"], fromFieldAligned(Ti));
    setBoundary(species["pressure"], fromFieldAligned(Pi));

    if (species.isSet("velocity")) {
      setBoundary(species["velocity"], fromFieldAligned(Vi));
    }

    if (species.isSet("momentum")) {
      setBoundary(species["momentum"], fromFieldAligned(NVi));
    }

    // Additional loss of energy through sheath
    // Note: Already includes any previously set sources
    set(species["energy_source"], fromFieldAligned(energy_source));
  }

  //////////////////////////////////////////////////////////////////
  // Electrons
  // This time adding energy sink, having calculated flow
  //

  Field3D electron_energy_source = electrons.isSet("energy_source")
    ? toFieldAligned(getNonFinal<Field3D>(electrons["energy_source"]))
    : zeroFrom(Ne);

  if (lower_y) {
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        auto i = indexAt(Ne, r.ind, mesh->ystart, jz);
        auto ip = i.yp();
        auto im = i.ym();

        const BoutReal nesheath = 0.5 * (Ne[im] + Ne[i]);
        const BoutReal tesheath = 0.5 * (Te[im] + Te[i]);  // electron temperature
        // Electron velocity into sheath (< 0). Calculated from ion flow
        const BoutReal vesheath = 0.5 * (Ve[im] + Ve[i]);

        // Take into account the flow of energy due to fluid flow
        // This is additional energy flux through the sheath
        // Note: Here this is negative because vesheath < 0
        BoutReal q = ((gamma_e - 1 - 1 / (electron_adiabatic - 1)) * tesheath
                      - 0.5 * Me * SQ(vesheath))
                     * nesheath * vesheath;

        // Multiply by cell area to get power
        BoutReal flux = q * (coord->J[i] + coord->J[im])
                        / (sqrt(coord->g_22[i]) + sqrt(coord->g_22[im]));

        // Divide by volume of cell to get energy loss rate (< 0)
        BoutReal power = flux / (coord->dy[i] * coord->J[i]);

#if CHECKLEVEL >= 1
        if (!std::isfinite(power)) {
	  throw BoutException("Non-finite power at {} : Te {} Ne {} Ve {}", i, tesheath, nesheath, vesheath);
	}
#endif

        electron_energy_source[i] += power;
      }
    }
  }
  if (upper_y) {
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        auto i = indexAt(Ne, r.ind, mesh->yend, jz);
        auto ip = i.yp();
        auto im = i.ym();

        const BoutReal nesheath = 0.5 * (Ne[ip] + Ne[i]);
        const BoutReal tesheath = 0.5 * (Te[ip] + Te[i]);  // electron temperature
        const BoutReal vesheath = 0.5 * (Ve[im] + Ve[i]);  // From ion flow

        // Take into account the flow of energy due to fluid flow
        // This is additional energy flux through the sheath
        // Note: Here this is positive because vesheath > 0
        BoutReal q = ((gamma_e - 1 - 1 / (electron_adiabatic - 1)) * tesheath
                      - 0.5 * Me * SQ(vesheath))
                     * nesheath * vesheath;

        // Multiply by cell area to get power
        BoutReal flux = q * (coord->J[i] + coord->J[ip])
                        / (sqrt(coord->g_22[i]) + sqrt(coord->g_22[ip]));

        // Divide by volume of cell to get energy loss rate (> 0)
        BoutReal power = flux / (coord->dy[i] * coord->J[i]);
#if CHECKLEVEL >= 1
	if (!std::isfinite(power)) {
	  throw BoutException("Non-finite power {} at {} : Te {} Ne {} Ve {} => q {}, flux {}",
                              power, i, tesheath, nesheath, vesheath, q, flux);
        }
#endif
        electron_energy_source[i] -= power;
      }
    }
  }

  // Set energy source (negative in cell next to sheath)
  set(electrons["energy_source"], fromFieldAligned(electron_energy_source));

  if (IS_SET_NOBOUNDARY(electrons["velocity"])) {
    setBoundary(electrons["velocity"], fromFieldAligned(Ve));
  }
  if (IS_SET_NOBOUNDARY(electrons["momentum"])) {
    setBoundary(electrons["momentum"], fromFieldAligned(NVe));
  }
}
