#include "../include/sheath_boundary_simple.hxx"

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

} // namespace

SheathBoundarySimple::SheathBoundarySimple(std::string name, Options& alloptions,
                                           Solver*) {
  AUTO_TRACE();

  Options& options = alloptions[name];

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

  gamma_e = options["gamma_e"]
    .doc("Electron sheath heat transmission coefficient")
    .withDefault(3.5);

  gamma_i = options["gamma_i"]
    .doc("Ion sheath heat transmission coefficient")
    .withDefault(3.5);

  sheath_ion_polytropic = options["sheath_ion_polytropic"]
         .doc("Ion polytropic coefficient in Bohm sound speed")
         .withDefault(1.0);

  lower_y = options["lower_y"].doc("Boundary on lower y?").withDefault<bool>(true);
  upper_y = options["upper_y"].doc("Boundary on upper y?").withDefault<bool>(true);

  always_set_phi =
      options["always_set_phi"]
          .doc("Always set phi field? Default is to only modify if already set")
          .withDefault<bool>(false);

  const Options& units = alloptions["units"];
  const BoutReal Tnorm = units["eV"];

  // Read wall voltage, convert to normalised units
  wall_potential = options["wall_potential"]
                       .doc("Voltage of the wall [Volts]")
                       .withDefault(Field3D(0.0))
                   / Tnorm;
  // Convert to field aligned coordinates
  wall_potential = toFieldAligned(wall_potential);
}

void SheathBoundarySimple::transform(Options& state) {
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

  Coordinates* coord = mesh->getCoordinates();

  //////////////////////////////////////////////////////////////////
  // Electrostatic potential
  // If phi is set, use free boundary condition
  // If phi not set, calculate assuming zero current
  Field3D phi;
  if (state.isSection("fields") and state["fields"].isSet("phi")) {
    phi = toFieldAligned(getNoBoundary<Field3D>(state["fields"]["phi"]));
  } else {
    // Calculate potential phi assuming zero current

    // Need to sum  n_i Z_i C_i over all ion species
    //
    // To avoid looking up species for every grid point, this
    // loops over the boundaries once per species.
    Field3D ion_sum = 0.0;

    // Iterate through charged ion species
    for (auto& kv : allspecies.getChildren()) {
      const Options& species = kv.second;

      if ((kv.first == "e") or !species.isSet("charge")
          or (get<BoutReal>(species["charge"]) == 0.0)) {
        continue; // Skip electrons and non-charged ions
      }

      const Field3D Ni = getNoBoundary<Field3D>(species["density"]);
      const Field3D Ti = getNoBoundary<Field3D>(species["temperature"]);
      const BoutReal Mi = getNoBoundary<BoutReal>(species["AA"]);
      const BoutReal Zi = getNoBoundary<BoutReal>(species["charge"]);

      if (lower_y) {
        // Sum values, put result in mesh->ystart

        for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            auto i = indexAt(Ni, r.ind, mesh->ystart, jz);
            auto ip = i.yp();

            // Free gradient of log density and temperature
            // This ensures that the guard cell values remain positive
            // exp( 2*log(N[i]) - log(N[ip]) )

            const BoutReal Ni_im = limitFree(Ni[ip], Ni[i]);
            const BoutReal Ti_im = limitFree(Ti[ip], Ti[i]);
            const BoutReal Te_im = limitFree(Te[ip], Te[i]);

            // Calculate sheath values at half-way points (cell edge)
            const BoutReal nisheath = 0.5 * (Ni_im + Ni[i]);
            const BoutReal tesheath =
                floor(0.5 * (Te_im + Te[i]), 1e-5); // electron temperature
            const BoutReal tisheath =
                floor(0.5 * (Ti_im + Ti[i]), 1e-5); // ion temperature

            // Sound speed squared
            BoutReal C_i_sq = (sheath_ion_polytropic * tisheath + Zi * tesheath) / Mi;

            ion_sum[i] += Zi * nisheath * sqrt(C_i_sq);
          }
        }
      }

      if (upper_y) {
        // Sum values, put results in mesh->yend

        for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            auto i = indexAt(Ni, r.ind, mesh->yend, jz);
            auto im = i.ym();

            const BoutReal Ni_ip = limitFree(Ni[im], Ni[i]);
            const BoutReal Ti_ip = limitFree(Ti[im], Ti[i]);
            const BoutReal Te_ip = limitFree(Te[im], Te[i]);

            // Calculate sheath values at half-way points (cell edge)
            const BoutReal nisheath = 0.5 * (Ni_ip + Ni[i]);
            const BoutReal tesheath =
                floor(0.5 * (Te_ip + Te[i]), 1e-5); // electron temperature
            const BoutReal tisheath =
                floor(0.5 * (Ti_ip + Ti[i]), 1e-5); // ion temperature

            BoutReal C_i_sq = (sheath_ion_polytropic * tisheath + Zi * tesheath) / Mi;

            ion_sum[i] += Zi * nisheath * sqrt(C_i_sq);
          }
        }
      }
    }

    phi.allocate();

    // ion_sum now contains the ion current, sum Z_i n_i C_i over all ion species
    // at mesh->ystart and mesh->yend indices
    if (lower_y) {
      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          auto i = indexAt(phi, r.ind, mesh->ystart, jz);
          auto ip = i.yp();

          const BoutReal Ne_im = limitFree(Ne[ip], Ne[i]);
          const BoutReal Te_im = limitFree(Te[ip], Te[i]);

          // Calculate sheath values at half-way points (cell edge)
          const BoutReal nesheath = 0.5 * (Ne_im + Ne[i]);
          const BoutReal tesheath = floor(0.5 * (Te_im + Te[i]), 1e-5);

          phi[i] =
              tesheath
              * log(sqrt(tesheath / (Me * TWOPI)) * (1. - Ge) * nesheath / ion_sum[i]);

          const BoutReal phi_wall = wall_potential[i];
          phi[i] += phi_wall; // Add bias potential

          phi[i.yp()] = phi[i.ym()] = phi[i]; // Constant into sheath
        }
      }
    }

    if (upper_y) {
      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          auto i = indexAt(phi, r.ind, mesh->yend, jz);
          auto im = i.ym();

          const BoutReal Ne_ip = limitFree(Ne[im], Ne[i]);
          const BoutReal Te_ip = limitFree(Te[im], Te[i]);

          // Calculate sheath values at half-way points (cell edge)
          const BoutReal nesheath = 0.5 * (Ne_ip + Ne[i]);
          const BoutReal tesheath = floor(0.5 * (Te_ip + Te[i]), 1e-5);

          phi[i] =
              tesheath
              * log(sqrt(tesheath / (Me * TWOPI)) * (1. - Ge) * nesheath / ion_sum[i]);

          const BoutReal phi_wall = wall_potential[i];
          phi[i] += phi_wall; // Add bias potential

          phi[i.yp()] = phi[i.ym()] = phi[i];
        }
      }
    }
  }

  //////////////////////////////////////////////////////////////////
  // Electrons

  Field3D electron_energy_source = electrons.isSet("energy_source")
    ? toFieldAligned(getNonFinal<Field3D>(electrons["energy_source"]))
    : zeroFrom(Ne);

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

        // Free boundary potential linearly extrapolated
        phi[im] = 2 * phi[i] - phi[ip];

        const BoutReal nesheath = 0.5 * (Ne[im] + Ne[i]);
        const BoutReal tesheath = 0.5 * (Te[im] + Te[i]); // electron temperature
        const BoutReal phi_wall = wall_potential[i];
        const BoutReal phisheath =
            floor(0.5 * (phi[im] + phi[i]), phi_wall); // Electron saturation at phi = phi_wall

        // Electron velocity into sheath (< 0)
        BoutReal vesheath =
	  -sqrt(tesheath / (TWOPI * Me)) * (1. - Ge) * exp(-(phisheath - phi_wall) / floor(tesheath, 1e-5));

        Ve[im] = 2 * vesheath - Ve[i];
        NVe[im] = 2 * Me * nesheath * vesheath - NVe[i];

        // Take into account the flow of energy due to fluid flow
        // This is additional energy flux through the sheath
        // Note: Here this is negative because vesheath < 0
        BoutReal q = ((gamma_e - 2.5) * tesheath
                      - 0.5 * Me * SQ(vesheath))
                     * nesheath * vesheath;

        // Multiply by cell area to get power
        BoutReal flux = q * (coord->J[i] + coord->J[im])
                        / (sqrt(coord->g_22[i]) + sqrt(coord->g_22[im]));

        // Divide by volume of cell to get energy loss rate (< 0)
        BoutReal power = flux / (coord->dy[i] * coord->J[i]);

        electron_energy_source[i] += power;
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

        // Free gradient of log electron density and temperature
        // This ensures that the guard cell values remain positive
        // exp( 2*log(N[i]) - log(N[ip]) )

        Ne[ip] = limitFree(Ne[im], Ne[i]);
        Te[ip] = limitFree(Te[im], Te[i]);
        Pe[ip] = limitFree(Pe[im], Pe[i]);

        // Free boundary potential linearly extrapolated.
        phi[ip] = 2 * phi[i] - phi[im];

        const BoutReal nesheath = 0.5 * (Ne[ip] + Ne[i]);
        const BoutReal tesheath = 0.5 * (Te[ip] + Te[i]); // electron temperature
        const BoutReal phi_wall = wall_potential[i];
        const BoutReal phisheath =
            floor(0.5 * (phi[ip] + phi[i]), phi_wall); // Electron saturation at phi = phi_wall

        // Electron velocity into sheath (> 0)
        BoutReal vesheath =
	  sqrt(tesheath / (TWOPI * Me)) * (1. - Ge) * exp(-(phisheath - phi_wall) / floor(tesheath, 1e-5));

        Ve[ip] = 2 * vesheath - Ve[i];
        NVe[ip] = 2. * Me * nesheath * vesheath - NVe[i];

        // Take into account the flow of energy due to fluid flow
        // This is additional energy flux through the sheath
        // Note: Here this is positive because vesheath > 0
        BoutReal q = ((gamma_e - 2.5) * tesheath
                      - 0.5 * Me * SQ(vesheath))
                     * nesheath * vesheath;

        // Multiply by cell area to get power
        BoutReal flux = q * (coord->J[i] + coord->J[ip])
                        / (sqrt(coord->g_22[i]) + sqrt(coord->g_22[ip]));

        // Divide by volume of cell to get energy loss rate (> 0)
        BoutReal power = flux / (coord->dy[i] * coord->J[i]);

        electron_energy_source[i] -= power;
      }
    }
  }

  // Set electron density and temperature, now with boundary conditions
  setBoundary(electrons["density"], fromFieldAligned(Ne));
  setBoundary(electrons["temperature"], fromFieldAligned(Te));
  setBoundary(electrons["pressure"], fromFieldAligned(Pe));

  // Set energy source (negative in cell next to sheath)
  // Note: electron_energy_source includes any sources previously set in other components
  set(electrons["energy_source"], fromFieldAligned(electron_energy_source));

  if (IS_SET_NOBOUNDARY(electrons["velocity"])) {
    setBoundary(electrons["velocity"], fromFieldAligned(Ve));
  }
  if (IS_SET_NOBOUNDARY(electrons["momentum"])) {
    setBoundary(electrons["momentum"], fromFieldAligned(NVe));
  }

  if (always_set_phi or (state.isSection("fields") and state["fields"].isSet("phi"))) {
    // Set the potential, including boundary conditions
    setBoundary(state["fields"]["phi"], fromFieldAligned(phi));
  }

  //////////////////////////////////////////////////////////////////
  // Iterate through all ions
  for (auto& kv : allspecies.getChildren()) {
    if (kv.first == "e") {
      continue; // Skip electrons
    }

    Options& species = allspecies[kv.first]; // Note: Need non-const

    // Ion charge
    const BoutReal Zi = species.isSet("charge") ? get<BoutReal>(species["charge"]) : 0.0;

    if (Zi == 0.0) {
      continue; // Neutral -> skip
    }

    // Characteristics of this species
    const BoutReal Mi = get<BoutReal>(species["AA"]);

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
          const BoutReal tesheath =
              floor(0.5 * (Te[im] + Te[i]), 1e-5); // electron temperature
          const BoutReal tisheath =
              floor(0.5 * (Ti[im] + Ti[i]), 1e-5); // ion temperature

          // Ion speed into sheath
          BoutReal C_i_sq = (sheath_ion_polytropic * tisheath + Zi * tesheath) / Mi;

          BoutReal visheath = -sqrt(C_i_sq); // Negative -> into sheath

          if (Vi[i] < visheath) {
            visheath = Vi[i];
          }

          // Set boundary conditions on flows
          Vi[im] = 2. * visheath - Vi[i];
          NVi[im] = 2. * Mi * nisheath * visheath - NVi[i];

          // Take into account the flow of energy due to fluid flow
          // This is additional energy flux through the sheath
          // Note: Here this is negative because visheath < 0
          BoutReal q =
              ((gamma_i - 2.5) * tisheath - 0.5 * Mi * C_i_sq)
              * nisheath * visheath;

          // Multiply by cell area to get power
          BoutReal flux = q * (coord->J[i] + coord->J[im])
                          / (sqrt(coord->g_22[i]) + sqrt(coord->g_22[im]));

          // Divide by volume of cell to get energy loss rate (< 0)
          BoutReal power = flux / (coord->dy[i] * coord->J[i]);

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
          const BoutReal tesheath =
              floor(0.5 * (Te[ip] + Te[i]), 1e-5); // electron temperature
          const BoutReal tisheath =
              floor(0.5 * (Ti[ip] + Ti[i]), 1e-5); // ion temperature

          // Ion speed into sheath
          BoutReal C_i_sq = (sheath_ion_polytropic * tisheath + Zi * tesheath) / Mi;

          BoutReal visheath = sqrt(C_i_sq); // Positive -> into sheath

          if (Vi[i] > visheath) {
            visheath = Vi[i];
          }

          // Set boundary conditions on flows
          Vi[ip] = 2. * visheath - Vi[i];
          NVi[ip] = 2. * Mi * nisheath * visheath - NVi[i];

          // Take into account the flow of energy due to fluid flow
          // This is additional energy flux through the sheath
          // Note: Here this is positive because visheath > 0
          BoutReal q =
              ((gamma_i - 2.5) * tisheath - 0.5 * C_i_sq * Mi)
              * nisheath * visheath;

          // Multiply by cell area to get power
          BoutReal flux = q * (coord->J[i] + coord->J[ip])
                          / (sqrt(coord->g_22[i]) + sqrt(coord->g_22[ip]));

          // Divide by volume of cell to get energy loss rate (> 0)
          BoutReal power = flux / (coord->dy[i] * coord->J[i]);
          ASSERT2(std::isfinite(power));

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
    // Note: energy_source already includes previously set values
    set(species["energy_source"], fromFieldAligned(energy_source));
  }
}
