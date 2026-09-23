#include <algorithm>

#include "../include/hermes_utils.hxx"
#include "../include/sheath_boundary_fci.hxx"
#include <bout/yboundary_regions.hxx>

#include "bout/constants.hxx"
#include "bout/mesh.hxx"
using bout::globals::mesh;

namespace {

/// Limited free gradient of log of a quantity
/// This ensures that the guard cell values remain positive
/// while also ensuring that the quantity never increases
///
///  fm  fc | fp
///         ^ boundary
///
/// exp( 2*log(fc) - log(fm) )
/// Mode 0: default (exponential extrapolation if decreases, Neumann if increases)
/// Mode 1: always exponential extrapolation
/// Mode 2: always linear extrapolation

BoutReal limitFree(BoutReal fm, BoutReal fc, BoutReal mode) {
  if ((fm < fc) && (mode == 0)) {
    return fc; // Neumann rather than increasing into boundary
  }
  if (fm < 1e-10) {
    return fc; // Low / no density condition
  }

  BoutReal fp = 0;
  if ((mode == 0) || (mode == 1)) {
    fp = SQ(fc) / fm; // Exponential
  } else if (mode == 2) {
    fp = 2.0 * fc - fm; // Linear
  } else {
    throw BoutException("Unknown boundary mode");
  }

  return fp; // Extrapolation

#if CHECKLEVEL >= 2
  if (!std::isfinite(fp)) {
    throw BoutException("SheathBoundary limitFree: {}, {} -> {}", fm, fc, fp);
  }
#endif

  return fp;
}


} // namespace

SheathBoundaryFci::SheathBoundaryFci(std::string name, Options& alloptions, Solver*)
    : NamedComponent(name, {
                               readIfSet("species:e:{e_whole_domain}"),
                               writeBoundary("species:e:{e_boundary}"),
                               readWrite("species:e:energy_source"),
                               readWrite("species:e:energy_flow_ylow"),
                               writeBoundaryIfSet("species:e:{e_optional}"),
                               writeBoundaryReadInteriorIfSet("species:e:pressure"),
                               readIfSet("species:{all_species}:charge"),
                               readOnly("species:{ions}:AA"),
                               readWrite("species:{ions}:energy_source"),
                               readWrite("species:{ions}:energy_flow_ylow"),
                               writeBoundary("species:{ions}:{ion_boundary}"),
                               writeBoundaryReadInteriorIfSet("species:{ions}:pressure"),
                               writeBoundaryIfSet("species:{ions}:{ion_optional}"),
                           }) {

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

  gamma_i =
      options["gamma_i"].doc("Ion sheath heat transmission coefficient").withDefault(3.5);

  sheath_ion_polytropic = options["sheath_ion_polytropic"]
                              .doc("Ion polytropic coefficient in Bohm sound speed")
                              .withDefault(1.0);

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
  wall_potential = wall_potential;

  no_flow = options["no_flow"]
                .doc("Set zero particle flow, keeping energy flow")
                .withDefault<bool>(false);

  density_boundary_mode =
      options["density_boundary_mode"]
          .doc("BC mode: 0=LimitFree, 1=ExponentialFree, 2=LinearFree")
          .withDefault<BoutReal>(1);

  pressure_boundary_mode =
      options["pressure_boundary_mode"]
          .doc("BC mode: 0=LimitFree, 1=ExponentialFree, 2=LinearFree")
          .withDefault<BoutReal>(1);

  temperature_boundary_mode =
      options["temperature_boundary_mode"]
          .doc("BC mode: 0=LimitFree, 1=ExponentialFree, 2=LinearFree")
          .withDefault<BoutReal>(1);

  diagnose = options["diagnose"]
                 .doc("Save additional output diagnostics")
                 .withDefault<bool>(false);

  substitutePermissions("e_whole_domain", {"AA", "charge"});
  substitutePermissions("e_boundary", {"density", "temperature"});
  substitutePermissions("e_optional", {"velocity", "momentum"});
  substitutePermissions("ion_boundary", {"density", "temperature"});
  substitutePermissions("ion_optional", {"velocity", "momentum"});
  setPermissions(always_set_phi ? writeBoundaryReadInteriorIfSet("fields:phi")
                                : writeBoundaryIfSet("fields:phi"));

  if (!mesh->isFci()) {
        throw BoutException("Using the Fci sheath variant while not using Fci. Please use sheath_boundary or sheath_boundary_simple instead!");
  }
  
}

void SheathBoundaryFci::transform_impl(GuardedOptions& state) {

  GuardedOptions allspecies = state["species"];
  GuardedOptions electrons = allspecies["e"];

  // Need electron properties
  // Not const because boundary conditions will be set
  Field3DParallel Ne =floor(GET_NOBOUNDARY(Field3D, electrons["density"]), 0.0);
  Field3DParallel Te = GET_NOBOUNDARY(Field3D, electrons["temperature"]);
  Field3DParallel Pe = IS_SET_NOBOUNDARY(electrons["pressure"])
                   ? getNoBoundary<Field3D>(electrons["pressure"])
                   : Te * Ne;
  
  
  // Mass, normalised to proton mass
  const BoutReal Me =
      IS_SET(electrons["AA"]) ? get<BoutReal>(electrons["AA"]) : SI::Me / SI::Mp;

  // This is for applying boundary conditions
  Field3DParallel Ve = IS_SET_NOBOUNDARY(electrons["velocity"])
                   ? getNoBoundary<Field3D>(electrons["velocity"])
                   : zeroFrom(Ne);

  Field3DParallel NVe = IS_SET_NOBOUNDARY(electrons["momentum"])
                    ? getNoBoundary<Field3D>(electrons["momentum"])
                    : zeroFrom(Ne);

  ASSERT2(Ne.hasParallelSlices());
  ASSERT2(Te.hasParallelSlices());
  ASSERT2(Pe.hasParallelSlices());
  ASSERT2(Ve.hasParallelSlices());
  ASSERT2(NVe.hasParallelSlices());


  
  Coordinates* coord = mesh->getCoordinates();

  //////////////////////////////////////////////////////////////////
  // Electrostatic potential
  // If phi is set, use free boundary condition
  // If phi not set, calculate assuming zero current
  Field3DParallel phi;
  if (IS_SET_NOBOUNDARY(state["fields"]["phi"])) {
    phi = getNoBoundary<Field3D>(state["fields"]["phi"]);
  } else {
    // Calculate potential phi assuming zero current

    // Need to sum  n_i Z_i C_i over all ion species
    //
    // To avoid looking up species for every grid point, this
    // loops over the boundaries once per species.
    ion_sum = 0.0;

    // Iterate through charged ion species
    for (auto& kv : allspecies.getChildren()) {
      const GuardedOptions species = kv.second;

      if ((kv.first == "e") or !species.isSet("charge")
          or (get<BoutReal>(species["charge"]) == 0.0)) {
        continue; // Skip electrons and non-charged ions
      }

      const Field3DParallel Ni = getNoBoundary<Field3D>(species["density"]);
      const Field3DParallel Ti = getNoBoundary<Field3D>(species["temperature"]);
      const BoutReal Mi = getNoBoundary<BoutReal>(species["AA"]);
      const BoutReal Zi = getNoBoundary<BoutReal>(species["charge"]);
      Field3DParallel Vi = species.isSet("velocity")
                       ? getNoBoundary<Field3D>(species["velocity"])
                       : zeroFrom(Ni);

      mesh->getCoordinates()->getYBoundary().iter([&](auto& pnt) {
	const auto& i = pnt.ind();
	const BoutReal Ni_im = limitFree(pnt.prev(Ni), pnt.current(Ni), density_boundary_mode);
	const BoutReal Ti_im = limitFree(pnt.prev(Ti), pnt.current(Ti), temperature_boundary_mode);
	const BoutReal Te_im = limitFree(pnt.prev(Te), pnt.current(Te), temperature_boundary_mode);

	const BoutReal nisheath = 0.5 * (Ni_im + pnt.current(Ni));
	const BoutReal tesheath =
	  floor(0.5 * (Te_im + pnt.current(Te)), 1e-5); // electron temperature
	const BoutReal tisheath =
	  floor(0.5 * (Ti_im + pnt.current(Ti)), 1e-5); // ion temperature

	// Sound speed squared
	BoutReal C_i_sq = (sheath_ion_polytropic * tisheath + Zi * tesheath) / Mi;

	BoutReal visheath;
	if (pnt.dir() > 0) {
	  visheath = std::max(pnt.current(Vi), pnt.dir() * sqrt(C_i_sq));
	} else {
	  visheath = std::min(pnt.current(Vi), pnt.dir() * sqrt(C_i_sq));
	}

	pnt.current(ion_sum) += pnt.dir() * Zi * nisheath * visheath;


      }); // End of yboundary.iter
      

    } // End of loop over all species

    phi.allocate();
    phi = 0.0;
    // ion_sum now contains the ion current, sum Z_i n_i C_i over all ion species
    // at mesh->ystart and mesh->yend indices


    mesh->getCoordinates()->getYBoundary().iter([&](auto& pnt) {
      const auto& i = pnt.ind();

      const BoutReal Ne_im = limitFree(pnt.prev(Ne), pnt.current(Ne), density_boundary_mode);
      const BoutReal Te_im = limitFree(pnt.prev(Te), pnt.current(Te), temperature_boundary_mode);

      const BoutReal nesheath = 0.5 * (Ne_im + pnt.current(Ne));
      const BoutReal tesheath = floor(0.5 * (Te_im + pnt.current(Te)), 1e-5);

      pnt.current(phi) = tesheath
	* log(sqrt(tesheath / (Me * TWOPI)) * (1. - Ge) * floor(nesheath, 1e-5)
	      / floor(pnt.current(ion_sum), 1e-5));

      const BoutReal phi_wall = wall_potential[i];
      pnt.current(phi) += phi_wall;
      pnt.next(phi) = pnt.current(phi);
	
    }); // End of yboundary.iter
    
    
    
    
  } // End of loop that calculates the phi if it is not set

  // Field to capture total sheath heat flux for diagnostics
  Field3D electron_sheath_power_ylow = zeroFrom(Ne);

  //////////////////////////////////////////////////////////////////
  // Electrons

  Field3D electron_energy_source =
      electrons.isSet("energy_source")
          ? getNonFinal<Field3D>(electrons["energy_source"])
          : zeroFrom(Ne);

  hflux_e = zeroFrom(electron_energy_source); // sheath heat flux for diagnostics

  mesh->getCoordinates()->getYBoundary().iter([&](auto& pnt) {
    const auto& i = pnt.ind();

    pnt.next(Ne) = limitFree(pnt.prev(Ne), pnt.current(Ne), density_boundary_mode);
    pnt.next(Te) = limitFree(pnt.prev(Te), pnt.current(Te), temperature_boundary_mode);
    pnt.next(Pe) = limitFree(pnt.prev(Pe), pnt.current(Pe), pressure_boundary_mode);

    pnt.next(phi) = 2.0 * pnt.current(phi) - pnt.prev(phi);

    const BoutReal nesheath = 0.5 * (pnt.next(Ne) + pnt.current(Ne));
    const BoutReal tesheath = 0.5 * (pnt.next(Te) + pnt.current(Te));
    const BoutReal phi_wall = wall_potential[i];
    const BoutReal phisheath = floor(
				     0.5 * (pnt.next(phi) + pnt.current(phi)), phi_wall); // Electron saturation at phi = phi_wall

    BoutReal vesheath = pnt.dir() * sqrt(tesheath / (TWOPI * Me)) * (1. - Ge)
                            * exp(-(phisheath - phi_wall) / floor(tesheath, 1e-5));

    pnt.next(Ve) = 2.0 * vesheath - pnt.current(Ve);
    pnt.next(NVe) = 2.0 * Me * nesheath * vesheath - pnt.current(NVe);

    if (abs(pnt.offset()) == 1) { // Only subtract flux when the cell is actually in direct contact with the sheath

      BoutReal q = gamma_e * tesheath * nesheath * vesheath;
      
      q -= (2.5 * tesheath + 0.5 * Me * SQ(vesheath)) * nesheath * vesheath;
      
      BoutReal flux;
      if (pnt.dir() > 0.0) {
	flux =  q * coord->cell_area_yhigh()[i];
      } else {
	flux =  q * coord->cell_area_ylow()[i];
      }
      
      // Divide by volume of cell to get energy loss rate (< 0)
      const BoutReal power = flux / coord->cell_volume()[i];

      electron_energy_source[i] -= pnt.dir() * power;
    }

  }); // End of yboundary.iter


  // Set electron density and temperature, now with boundary conditions
  // Note: Clear parallel slices because they do not contain boundary conditions.
  setBoundary(electrons["density"], Field3D{Ne});
  setBoundary(electrons["temperature"], Field3D{Te});
  setBoundary(electrons["pressure"], Field3D{Pe});

  // Set energy source (negative in cell next to sheath)
  // Note: electron_energy_source includes any sources previously set in other components
  set(electrons["energy_source"], electron_energy_source);


  if (IS_SET_NOBOUNDARY(electrons["velocity"])) {
    setBoundary(electrons["velocity"], Field3D{Ve});
  }
  if (IS_SET_NOBOUNDARY(electrons["momentum"])) {
    setBoundary(electrons["momentum"], Field3D{NVe});
  }

  if (always_set_phi or (state.isSection("fields") and state["fields"].isSet("phi"))) {
    // Set the potential, including boundary conditions
    setBoundary(state["fields"]["phi"], Field3D{phi});
  }

  //////////////////////////////////////////////////////////////////
  // Iterate through all ions
  for (auto& kv : allspecies.getChildren()) {
    if (kv.first == "e") {
      continue; // Skip electrons
    }

    GuardedOptions species = allspecies[kv.first]; // Note: Need non-const

    // Ion charge
    const BoutReal Zi = species.isSet("charge") ? get<BoutReal>(species["charge"]) : 0.0;

    if (Zi == 0.0) {
      continue; // Neutral -> skip
    }

    // Characteristics of this species
    const BoutReal Mi = get<BoutReal>(species["AA"]);

    // Density and temperature boundary conditions will be imposed (free)
    Field3DParallel Ni = (floor(getNoBoundary<Field3D>(species["density"]), 0.0));
    Field3DParallel Ti = getNoBoundary<Field3D>(species["temperature"]);
    Field3DParallel Pi = species.isSet("pressure")
                     ? getNoBoundary<Field3D>(species["pressure"])
                     : Ni * Ti;

    // Get the velocity and momentum
    // These will be modified at the boundaries
    // and then put back into the state
    Field3DParallel Vi = species.isSet("velocity")
                     ? getNoBoundary<Field3D>(species["velocity"])
                     : zeroFrom(Ni);
    Field3DParallel NVi = species.isSet("momentum")
                      ? getNoBoundary<Field3D>(species["momentum"])
                      : Mi * Ni * Vi;

    // Energy source will be modified in the domain
    Field3D energy_source =
        species.isSet("energy_source")
            ? getNonFinal<Field3D>(species["energy_source"])
            : zeroFrom(Ni);

    // Initialise sheath ion heat flux. This will be created for each species
    // saved in diagnostics struct and then destroyed and re-created for next species


    mesh->getCoordinates()->getYBoundary().iter([&](auto& pnt) {
      const auto& i = pnt.ind();

      // Free gradient of log electron density and temperature
      // This ensures that the guard cell values remain positive
      // exp( 2*log(N[i]) - log(N[ip]) )
      
      pnt.next(Ni) = limitFree(pnt.prev(Ni), pnt.current(Ni), density_boundary_mode);
      pnt.next(Ti) = limitFree(pnt.prev(Ti), pnt.current(Ti), temperature_boundary_mode);
      pnt.next(Pi) = limitFree(pnt.prev(Pi), pnt.current(Pi), pressure_boundary_mode);

      // Calculate sheath values at half-way points (cell edge)
      const BoutReal nisheath = 0.5 * (pnt.current(Ni) + pnt.next(Ni));
      const BoutReal tesheath = floor(pnt.current(Te) + pnt.next(Te), 1e-5);
      const BoutReal tisheath =	floor(pnt.current(Ti) + pnt.next(Ti), 1e-5);

      // Ion speed into sheath
      BoutReal C_i_sq = (sheath_ion_polytropic * tisheath + Zi * tesheath) / Mi;

      BoutReal visheath;
      if (pnt.dir() > 0) {
	visheath = std::max(pnt.current(Vi), pnt.dir() * sqrt(C_i_sq));
      } else {
	visheath = std::min(pnt.current(Vi), pnt.dir() * sqrt(C_i_sq));
      }

      pnt.next(Vi) = 2.0 * visheath - pnt.current(Vi);
      pnt.next(NVi) = 2.0 * Mi * nisheath * visheath - pnt.current(NVi);

      if (abs(pnt.offset()) == 1) { // Only subtract flux when the cell is actually in direct contact with the sheath
	
	BoutReal q = gamma_i * tisheath * nisheath * visheath;
	q -= (2.5 * tisheath + 0.5 * Mi * SQ(visheath)) * nisheath * visheath;
	
	BoutReal flux;
	if (pnt.dir() > 0.0) {
	  flux =  q * coord->cell_area_yhigh()[i];
	} else {
	  flux =  q * coord->cell_area_ylow()[i];
	}

	// Divide by volume of cell to get energy loss rate (< 0)
	
	const BoutReal power = flux / coord->cell_volume()[i];
	energy_source[i] -= pnt.dir() * power;
  
      }
      
    }); // End of yboundary.iter

    
    // Finished boundary conditions for this species
    // Put the modified fields back into the state.

    setBoundary(species["density"], Field3D{Ni});
    setBoundary(species["temperature"], Field3D{Ti});
    setBoundary(species["pressure"], Field3D{Pi});

    if (species.isSet("velocity")) {
      setBoundary(species["velocity"], Field3D{Vi});
    }

    if (species.isSet("momentum")) {
      setBoundary(species["momentum"], Field3D{NVi});
    }

    // Additional loss of energy through sheath
    // Note: energy_source already includes previously set values
    set(species["energy_source"], energy_source);

  }
}

void SheathBoundaryFci::outputVars(Options& state) {
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

  if (diagnose) {
    /// Iterate through the first species in each collision pair
    const std::map<std::string, Options>& level1 = diagnostics.getChildren();
    for (auto s1 = std::begin(level1); s1 != std::end(level1); ++s1) {
      auto species_name = s1->first;
      const Options& section = diagnostics[species_name];

      set_with_attrs(state[{"E" + species_name + "_sheath"}],
                     getNonFinal<Field3D>(section["energy_source"]),
                     {{"time_dimension", "t"},
                      {"units", "W / m^3"},
                      {"conversion", Pnorm * Omega_ci},
                      {"standard_name", "energy source"},
                      {"long_name", species_name + " sheath energy source"},
                      {"source", "sheath_boundary_fci"}});

    }
  }
}
