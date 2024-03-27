#include "../include/sheath_boundary_parallel.hxx"

#include <bout/output_bout_types.hxx>

#include "bout/constants.hxx"
#include "bout/mesh.hxx"

#include "bout/parallel_boundary_region.hxx"

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
    throw BoutException("SheathBoundaryParallel limitFree: {}, {} -> {}", fm, fc, fp);
  }
#endif

  return fp;
}

BoutReal limitFree(const Field3D& f, const BoundaryRegionParIter& pnt) {
  if (pnt.valid() > 0) {
    return limitFree(pnt.yprev(f), f[pnt.ind()]);
  }
  return f[pnt.ind()];
}

BoutReal limitFree(const Field3D& f, const BoundaryRegionIter& pnt) {
  return limitFree(pnt.yprev(f), f[pnt.ind()]);
}
}

SheathBoundaryParallel::SheathBoundaryParallel(std::string name, Options &alloptions, Solver *) {
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

  bool lower_y = options["lower_y"].doc("Boundary on lower y?").withDefault<bool>(true);
  bool upper_y = options["upper_y"].doc("Boundary on upper y?").withDefault<bool>(true);
  bool outer_x = options["outer_x"].doc("Boundary on inner y?").withDefault<bool>(true);
  bool inner_x = options["inner_x"].doc("Boundary on outer y?").withDefault<bool>(false);
  if (wall_potential.isFci()) {
    if (outer_x) {
      for (auto& bndry : mesh->getBoundariesPar(BoundaryParType::xout)) {
        boundary_regions_par.push_back(bndry);
      }
    }
    if (inner_x) {
      for (auto& bndry : mesh->getBoundariesPar(BoundaryParType::xin)) {
        boundary_regions_par.push_back(bndry);
      }
    }
  } else {
    for (auto& bndry : mesh->getBoundaries()) {
      if (upper_y && bndry->location == BndryLoc::yup) {
        boundary_regions.push_back(bndry);
      }
      if (lower_y && bndry->location == BndryLoc::ydown) {
        boundary_regions.push_back(bndry);
      }
    }
  }

  // Note: wall potential at the last cell before the boundary is used,
  // not the value at the boundary half-way between cells. This is due
  // to how twist-shift boundary conditions and non-aligned inputs are
  // treated; using the cell boundary gives incorrect results.

  floor_potential = options["floor_potential"]
                        .doc("Apply a floor to wall potential when calculating Ve?")
                        .withDefault<bool>(true);
}

void SheathBoundaryParallel::transform(Options &state) {
  AUTO_TRACE();

  Options& allspecies = state["species"];
  Options& electrons = allspecies["e"];

  // Need electron properties
  // Not const because boundary conditions will be set
  Field3D Ne = (floor(GET_NOBOUNDARY(Field3D, electrons["density"]), 0.0));
  Field3D Te = (GET_NOBOUNDARY(Field3D, electrons["temperature"]));
  Field3D Pe = IS_SET_NOBOUNDARY(electrons["pressure"])
    ? (getNoBoundary<Field3D>(electrons["pressure"]))
    : Te * Ne;

  // Ratio of specific heats
  const BoutReal electron_adiabatic =
      IS_SET(electrons["adiabatic"]) ? get<BoutReal>(electrons["adiabatic"]) : 5. / 3;

  // Mass, normalised to proton mass
  const BoutReal Me =
      IS_SET(electrons["AA"]) ? get<BoutReal>(electrons["AA"]) : SI::Me / SI::Mp;

  // This is for applying boundary conditions
  Field3D Ve = IS_SET_NOBOUNDARY(electrons["velocity"])
    ? (getNoBoundary<Field3D>(electrons["velocity"]))
    : zeroFrom(Ne);

  Field3D NVe = IS_SET_NOBOUNDARY(electrons["momentum"])
    ? (getNoBoundary<Field3D>(electrons["momentum"]))
    : zeroFrom(Ne);

  Coordinates *coord = mesh->getCoordinates();

  //////////////////////////////////////////////////////////////////
  // Electrostatic potential
  // If phi is set, use free boundary condition
  // If phi not set, calculate assuming zero current
  Field3D phi;
  if (IS_SET_NOBOUNDARY(state["fields"]["phi"])) {
    phi = (getNoBoundary<Field3D>(state["fields"]["phi"]));
  } else {
    // Calculate potential phi assuming zero current
    // Note: This is equation (22) in Tskhakaya 2005, with I = 0

    // Need to sum  s_i Z_i C_i over all ion species
    //
    // To avoid looking up species for every grid point, this
    // loops over the boundaries once per species.
    Field3D ion_sum {zeroFrom(Ne)};
    phi = emptyFrom(Ne); // So phi is field aligned

    // Iterate through charged ion species
    for (auto& kv : allspecies.getChildren()) {
      Options& species = allspecies[kv.first];

      if ((kv.first == "e") or !IS_SET(species["charge"])
          or (get<BoutReal>(species["charge"]) == 0.0)) {
        continue; // Skip electrons and non-charged ions
      }

      const Field3D Ni = (floor(GET_NOBOUNDARY(Field3D, species["density"]), 0.0));
      const Field3D Ti = (GET_NOBOUNDARY(Field3D, species["temperature"]));
      const BoutReal Mi = GET_NOBOUNDARY(BoutReal, species["AA"]);
      const BoutReal Zi = GET_NOBOUNDARY(BoutReal, species["charge"]);

      const BoutReal adiabatic = IS_SET(species["adiabatic"])
                                     ? get<BoutReal>(species["adiabatic"])
                                     : 5. / 3; // Ratio of specific heats (ideal gas)

      for (auto* region : boundary_regions) {
        for (auto& pnt : region) {
          const auto& i = pnt.ind();
          BoutReal s_i =
              clip(0.5 * pnt.extrapolate_next_o2([&, Ni, Ne](int yoffset, Ind3D ind) {
                return Ni.ynext(yoffset)[ind] / Ne.ynext(yoffset)[ind];
              }),
                   0.0, 1.0);

          if (!std::isfinite(s_i)) {
            s_i = 1.0;
          }
          BoutReal te = Te[i];
          BoutReal ti = Ti[i];

          // Equation (9) in Tskhakaya 2005
          BoutReal grad_ne = pnt.extrapolate_grad_o2(Ne);
          BoutReal grad_ni = pnt.extrapolate_grad_o2(Ni);

          // Note: Needed to get past initial conditions, perhaps
          // transients but this shouldn't happen in steady state
          if (fabs(grad_ni) < 1e-3) {
            grad_ni = grad_ne = 1e-3; // Remove kinetic correction term
          }

          BoutReal C_i_sq =
              clip((adiabatic * ti + Zi * s_i * te * grad_ne / grad_ni) / Mi, 0,
                   100); // Limit for e.g. Ni zero gradient

          // Note: Vzi = C_i * sin(α)
          ion_sum[pnt.ind()] += s_i * Zi * sin_alpha * sqrt(C_i_sq);
        }
      }
    }

    phi.allocate();

    // ion_sum now contains  sum  s_i Z_i C_i over all ion species
    // at mesh->ystart and mesh->yend indices
    for (const auto& region: boundary_regions) {
      for (const auto& pnt: *region) {
	auto i = pnt.ind();

	if (Te[i] <= 0.0) {
	  phi[i] = 0.0;
	} else {
	  phi[i] = Te[i] * log(sqrt(Te[i] / (Me * TWOPI)) * (1. - Ge) / ion_sum[i]);
	}

	const BoutReal phi_wall = wall_potential[i];
	phi[i] += phi_wall; // Add bias potential

	phi.ynext(1)[i.yp()] = phi.ynext(-1)[i.ym()] = phi[i]; // Constant into sheath
      }
    }
  }

  //////////////////////////////////////////////////////////////////
  // Electrons

  Field3D electron_energy_source = electrons.isSet("energy_source")
    ? (getNonFinal<Field3D>(electrons["energy_source"]))
    : zeroFrom(Ne);

  for (const auto& region: boundary_regions) {
    for (const auto& pnt: *region) {
      auto i = pnt.ind();

      // Free gradient of log electron density and temperature
      // Limited so that the values don't increase into the sheath
      // This ensures that the guard cell values remain positive
      // exp( 2*log(N[i]) - log(N[ip]) )

      pnt.ynext(Ne) = limitFree(Ne, pnt);
      pnt.ynext(Te) = limitFree(Te, pnt);
      pnt.ynext(Pe) = limitFree(Pe, pnt);

      // Free boundary potential linearly extrapolated
      pnt.ynext(phi) = phi[i] + pnt.extrapolate_grad_o2(phi);

      const BoutReal nesheath = pnt.interpolate_sheath(Ne);
      const BoutReal tesheath = pnt.interpolate_sheath(Te);  // electron temperature
      const BoutReal phi_wall = wall_potential[i];

      const BoutReal phisheath = floor_potential ? floor(
            pnt.interpolate_sheath(phi), phi_wall) // Electron saturation at phi = phi_wall
	    : pnt.interpolate_sheath(phi);

      // Electron sheath heat transmission
      const BoutReal gamma_e = floor(2 / (1. - Ge) + (phisheath - phi_wall) / floor(tesheath, 1e-5), 0.0);

      // Electron velocity into sheath (< 0)
      const BoutReal vesheath = (tesheath < 1e-10) ?
          0.0 :
          pnt.dir * sqrt(tesheath / (TWOPI * Me)) * (1. - Ge) * exp(-(phisheath - phi_wall) / tesheath);

      pnt.dirichlet_o2(Ve, vesheath);
      pnt.dirichlet_o2(NVe, Me * nesheath * vesheath);

      // Take into account the flow of energy due to fluid flow
      // This is additional energy flux through the sheath
      // Note: sign depends on sign of vesheath
      BoutReal q = ((gamma_e - 1 - 1 / (electron_adiabatic - 1)) * tesheath
                      - 0.5 * Me * SQ(vesheath))
                     * nesheath * vesheath;

      BoutReal flux;
      // Multiply by cell area to get power
      if (! Ve.isFci()) {
	const auto ip = i.yp(pnt.dir);
        flux = q * (coord->J[i] + coord->J[ip])
                        / (sqrt(coord->g_22[i]) + sqrt(coord->g_22[ip]));
      } else {
	flux = q * coord->J[i] / sqrt(coord->g_22[i]);
      }
      // Divide by volume of cell to get energy loss rate (sign depending on vesheath)
      BoutReal power = flux / (coord->dy[i] * coord->J[i]);

#if CHECKLEVEL >= 1
      if (!std::isfinite(power)) {
	throw BoutException("Non-finite power {} at {} : Te {} Ne {} Ve {} phi {}, {} => q {}, flux {}",
			    power, i, tesheath, nesheath, vesheath, phi[i], phi[im], q, flux);
      }
#endif

      electron_energy_source[i] -= pnt.dir * power;
    }
  }

  // Set electron density and temperature, now with boundary conditions
  setBoundary(electrons["density"], (Ne));
  setBoundary(electrons["temperature"], (Te));
  setBoundary(electrons["pressure"], (Pe));

  // Add energy source (negative in cell next to sheath)
  // Note: already includes previously set sources
  set(electrons["energy_source"], (electron_energy_source));

  if (IS_SET_NOBOUNDARY(electrons["velocity"])) {
    setBoundary(electrons["velocity"], (Ve));
  }
  if (IS_SET_NOBOUNDARY(electrons["momentum"])) {
    setBoundary(electrons["momentum"], (NVe));
  }

  if (always_set_phi or IS_SET_NOBOUNDARY(state["fields"]["phi"])) {
    // Set the potential, including boundary conditions
    phi = (phi);
    //output.write("-> phi {}\n", phi(10, mesh->yend+1, 0));
    setBoundary(state["fields"]["phi"], phi);
  }

  //////////////////////////////////////////////////////////////////
  // Iterate through all ions
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
    Field3D Ni = (floor(getNoBoundary<Field3D>(species["density"]), 0.0));
    Field3D Ti = (getNoBoundary<Field3D>(species["temperature"]));
    Field3D Pi = species.isSet("pressure")
      ? (getNoBoundary<Field3D>(species["pressure"]))
      : Ni * Ti;

    // Get the velocity and momentum
    // These will be modified at the boundaries
    // and then put back into the state
    Field3D Vi = species.isSet("velocity")
      ? (getNoBoundary<Field3D>(species["velocity"]))
      : zeroFrom(Ni);
    Field3D NVi = species.isSet("momentum")
      ? (getNoBoundary<Field3D>(species["momentum"]))
      : Mi * Ni * Vi;

    // Energy source will be modified in the domain
    Field3D energy_source = species.isSet("energy_source")
      ? (getNonFinal<Field3D>(species["energy_source"]))
      : zeroFrom(Ni);

    for (auto& region : boundary_regions) {
      for (const auto& pnt: *region) {
	
	auto i = pnt.ind();

	// Free gradient of log electron density and temperature
	// This ensures that the guard cell values remain positive
	// exp( 2*log(N[i]) - log(N[ip]) )

	pnt.ynext(Ni) = limitFree(Ni, pnt);
	pnt.ynext(Ti) = limitFree(Ti, pnt);
	pnt.ynext(Pi) = limitFree(Pi, pnt);

	// Calculate sheath values at half-way points (cell edge)
	const BoutReal nesheath = pnt.extrapolate_sheath_o2(Ne);
	const BoutReal nisheath = pnt.extrapolate_sheath_o2(Ni);
	const BoutReal tesheath = floor(pnt.extrapolate_sheath_o2(Te), 1e-5);  // electron temperature
	const BoutReal tisheath = floor(pnt.extrapolate_sheath_o2(Ti), 1e-5);  // ion temperature

	// Ion sheath heat transmission coefficient
	//
	// 1 / (1 + ∂_{ln n_e} ln s_i = s_i * ∂n_e / (s_i * ∂n_e + ∂ n_i) 
	BoutReal s_i = (nesheath > 1e-5) ? nisheath / nesheath : 0.0; // Concentration
	BoutReal grad_ne = pnt.extrapolate_grad_o2(Ne);
	BoutReal grad_ni = pnt.extrapolate_grad_o2(Ni);

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
	pnt.dirichlet_o2(Vi, visheath);
	pnt.dirichlet_o2(NVi, Mi * nisheath * visheath);

	// Take into account the flow of energy due to fluid flow
	// This is additional energy flux through the sheath
	// Note: Sign depends on sign of visheath
	BoutReal q =
	  ((gamma_i - 1 - 1 / (adiabatic - 1)) * tisheath - 0.5 * C_i_sq * Mi)
	  * nisheath * visheath;

	if (q * pnt.dir < 0.0) {
	  q = 0.0;
	}

	// Multiply by cell area to get power
	BoutReal flux;
	if (! Ti.isFci()) {
	  const auto ip = i.yp(pnt.dir);
	  flux = q * (coord->J[i] + coord->J[ip])
                          / (sqrt(coord->g_22[i]) + sqrt(coord->g_22[ip]));
	} else {
	  flux = q * coord->J[i] / sqrt(coord->g_22[i]);
	}

	// Divide by volume of cell to get energy loss rate
	BoutReal power = flux / (coord->dy[i] * coord->J[i]);
	ASSERT1(std::isfinite(power));
	ASSERT2(power * pnt.dir >= 0.0);

	energy_source[i] -= power * pnt.dir; // Note: Sign negative because power * direction > 0
      }
    }

    // Finished boundary conditions for this species
    // Put the modified fields back into the state.
    setBoundary(species["density"], (Ni));
    setBoundary(species["temperature"], (Ti));
    setBoundary(species["pressure"], (Pi));

    if (species.isSet("velocity")) {
      setBoundary(species["velocity"], (Vi));
    }

    if (species.isSet("momentum")) {
      setBoundary(species["momentum"], (NVi));
    }

    // Additional loss of energy through sheath
    // Note: Already includes previously set sources
    set(species["energy_source"], (energy_source));
  }
}
