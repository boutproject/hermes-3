#include "../include/recycling_fci.hxx"

#include <bout/utils.hxx> // for trim, strsplit
#include "../include/hermes_utils.hxx"  // For indexAt
#include "../include/hermes_utils.hxx"  // For indexAt
#include <bout/coordinates.hxx>
#include <bout/mesh.hxx>
#include <bout/constants.hxx>
#include <bout/yboundary_regions.hxx>

using bout::globals::mesh;

RecyclingFCI::RecyclingFCI(std::string name, Options& alloptions, Solver*) {
  AUTO_TRACE();

  const Options& units = alloptions["units"];
  const BoutReal Tnorm = units["eV"];

  Options& options = alloptions[name];

  // init parallel bc iterator
  yboundary.init(options);
 
  auto species_list = strsplit(options["species"]
                                   .doc("Comma-separated list of species to recycle")
                                   .as<std::string>(),
                               ',');
  
  // Neutral pump
  // Mark cells as having a pump by setting the Field2D is_pump to 1 in the grid file
  // Works only on SOL and PFR edges, where it locally modifies the recycle multiplier to the pump albedo
  is_pump = 0.0;
  for (const auto& species : species_list) {
    std::string from = trim(species, " \t\r()"); // The species name in the list

    if (from.empty())
      continue; // Missing

    // Get the options for this species
    Options& from_options = alloptions[from];
    std::string to = from_options["recycle_as"]
                         .doc("Name of the species to recycle into")
                         .as<std::string>();

    diagnose =
      from_options["diagnose"].doc("Save additional diagnostics?").withDefault<bool>(false);
    
    BoutReal target_recycle_multiplier =
        from_options["target_recycle_multiplier"]
            .doc("Multiply the target recycled flux by this factor. Should be >=0 and <= 1")
            .withDefault<BoutReal>(1.0);

    BoutReal sol_recycle_multiplier =
        from_options["sol_recycle_multiplier"]
            .doc("Multiply the sol recycled flux by this factor. Should be >=0 and <= 1")
            .withDefault<BoutReal>(1.0);

    BoutReal pfr_recycle_multiplier =
        from_options["pfr_recycle_multiplier"]
            .doc("Multiply the pfr recycled flux by this factor. Should be >=0 and <= 1")
            .withDefault<BoutReal>(1.0);

    BoutReal pump_recycle_multiplier =
      from_options["pump_recycle_multiplier"]
          .doc("Multiply the pump boundary recycling flux by this factor (like albedo). Should be >=0 and <= 1")
          .withDefault<BoutReal>(1.0);

    BoutReal target_recycle_energy = from_options["target_recycle_energy"]
                                  .doc("Fixed energy of the recycled particles at target [eV]")
                                  .withDefault<BoutReal>(3.0)
                              / Tnorm; // Normalise from eV

    BoutReal sol_recycle_energy = from_options["sol_recycle_energy"]
                                  .doc("Fixed energy of the recycled particles at sol [eV]")
                                  .withDefault<BoutReal>(3.0)
                              / Tnorm; // Normalise from eV

    BoutReal pfr_recycle_energy = from_options["pfr_recycle_energy"]
                                  .doc("Fixed energy of the recycled particles at pfr [eV]")
                                  .withDefault<BoutReal>(3.0)
                              / Tnorm; // Normalise from eV

    BoutReal target_fast_recycle_fraction =
        from_options["target_fast_recycle_fraction"]
            .doc("Fraction of ions undergoing fast reflection at target")
            .withDefault<BoutReal>(0);

    BoutReal pfr_fast_recycle_fraction =
        from_options["pfr_fast_recycle_fraction"]
            .doc("Fraction of ions undergoing fast reflection at pfr")
            .withDefault<BoutReal>(0);

    BoutReal sol_fast_recycle_fraction =
        from_options["sol_fast_recycle_fraction"]
            .doc("Fraction of ions undergoing fast reflection at sol")
            .withDefault<BoutReal>(0);

    BoutReal target_fast_recycle_energy_factor =
        from_options["target_fast_recycle_energy_factor"]
            .doc("Fraction of energy retained by fast recycled neutrals at target")
            .withDefault<BoutReal>(0);

    BoutReal sol_fast_recycle_energy_factor =
        from_options["sol_fast_recycle_energy_factor"]
            .doc("Fraction of energy retained by fast recycled neutrals at sol")
            .withDefault<BoutReal>(0);

    BoutReal pfr_fast_recycle_energy_factor =
        from_options["pfr_fast_recycle_energy_factor"]
            .doc("Fraction of energy retained by fast recycled neutrals at pfr")
            .withDefault<BoutReal>(0);

    if ((target_recycle_multiplier < 0.0) or (target_recycle_multiplier > 1.0)
    or (sol_recycle_multiplier < 0.0) or (sol_recycle_multiplier > 1.0)
    or (pfr_recycle_multiplier < 0.0) or (pfr_recycle_multiplier > 1.0)
    or (pump_recycle_multiplier < 0.0) or (pump_recycle_multiplier > 1.0)) {
      throw BoutException("All recycle multipliers must be betweeen 0 and 1");
    }

    // Populate recycling channel vector
    channels.push_back({
      from, to, 
      target_recycle_multiplier, sol_recycle_multiplier, pfr_recycle_multiplier, pump_recycle_multiplier,
      target_recycle_energy, sol_recycle_energy, pfr_recycle_energy,
      target_fast_recycle_fraction, pfr_fast_recycle_fraction, sol_fast_recycle_fraction,
      target_fast_recycle_energy_factor, sol_fast_recycle_energy_factor, pfr_fast_recycle_energy_factor});

    // Boolean flags for enabling recycling in different regions
    target_recycle = from_options["target_recycle"]
                   .doc("Recycling in the targets?")
                   .withDefault<bool>(false);

    sol_recycle = from_options["sol_recycle"]
                   .doc("Recycling in the SOL edge?")
                   .withDefault<bool>(false);

    pfr_recycle = from_options["pfr_recycle"]
                   .doc("Recycling in the PFR edge?")
                   .withDefault<bool>(false);

    neutral_pump = from_options["neutral_pump"]
                   .doc("Neutral pump enabled? Note, need location in grid file")
                   .withDefault<bool>(false);                 
  }


  // Get metric tensor components
  Coordinates* coord = mesh->getCoordinates();

  // set parallel slices of grid spacing (assumed uniform)
  // FCI needs the yup/down fields for grid spacing.
  auto& dx = coord->dx;
  auto& dy = coord->dy;
  auto& dz = coord->dz;
  
}

void RecyclingFCI::transform(Options& state) {
  AUTO_TRACE();

  // Get metric tensor components
  Coordinates* coord = mesh->getCoordinates();

  
  for (auto& channel : channels) {
    const Options& species_from = state["species"][channel.from];

    const Field3D N = get<Field3D>(species_from["density"]);
    const Field3D V = get<Field3D>(species_from["velocity"]); // Parallel flow velocity
    const Field3D T = get<Field3D>(species_from["temperature"]); // Ion temperature
    const BoutReal AAf = get<BoutReal>(species_from["AA"]);

    
    Options& species_to = state["species"][channel.to];
    const Field3D Nn = get<Field3D>(species_to["density"]);
    const Field3D Pn = get<Field3D>(species_to["pressure"]);
    const Field3D Tn = get<Field3D>(species_to["temperature"]);
    const BoutReal AAn = get<BoutReal>(species_to["AA"]);

    bool dissipative = species_from["dissipative"].withDefault(true); 
    
    // Recycling particle and energy sources will be added to these global sources 
    // which are then passed to the density and pressure equations
    density_source = species_to.isSet("density_source")
                                 ? getNonFinal<Field3D>(species_to["density_source"])
                                 : 0.0;
    energy_source = species_to.isSet("energy_source")
                                ? getNonFinal<Field3D>(species_to["energy_source"])
                                : 0.0;

    // Recycling at the divertor target plates
    if (target_recycle) {


      channel.target_recycle_density_source = 0;
      channel.target_recycle_energy_source = 0;

      Field3D fastest_wave;
      if (dissipative) {
	if (state.isSet("fastest_wave")) {
	  fastest_wave = get<Field3D>(state["fastest_wave"]);
	} else {
	  fastest_wave = sqrt(T / AAf);
	}
      }
      
      // Y boundaries
      yboundary.iter_pnts([&](auto& pnt) {
	// BoutReal flux = pnt.dir * pnt.interpolate_sheath_o1(N) * pnt.interpolate_sheath_o1(V);
	BoutReal flux = 0.0;
	if (dissipative) {
	  const BoutReal amax = BOUTMAX(pnt.ythis(fastest_wave),
					fabs(pnt.ythis(V)),
					fabs(pnt.ynext(V)),
					fabs(pnt.yprev(V)));
	  flux = pnt.dir * (0.5 * (pnt.ythis(N) * (pnt.ythis(V) + amax) + pnt.ynext(N) * (pnt.ynext(N) - amax)));
	} else {
	  flux = pnt.dir * (0.25 * (pnt.ythis(N) + pnt.ynext(N)) * (pnt.ythis(V) + pnt.ynext(V)));
	}
	    
	if (flux < 0.0) {
	  flux = 0.0;
	}
	
	BoutReal area = (pnt.ythis(coord->J) + pnt.ynext(coord->J)) / (pnt.ythis(coord->dy) * (sqrt(pnt.ythis(coord->g_22)) + sqrt(pnt.ynext(coord->g_22))));

	BoutReal flow = channel.target_multiplier * flux * area;

	BoutReal volume  = pnt.ythis(coord->J);

	// Calculate sources in the final cell [m^-3 s^-1]
	if (pnt.abs_offset() == 1){
	  pnt.ythis(channel.target_recycle_density_source) += flow / volume;    // For diagnostic
	  pnt.ythis(density_source) += flow / volume;         // For use in solver
	}
        
        BoutReal recycle_energy_flow = flow  * channel.target_energy;   // Thermal recycling par
	

	// Divide heat flow in [W] by cell volume to get source in [m^-3 s^-1]
	if (pnt.abs_offset() == 1) {
	  pnt.ythis(channel.target_recycle_energy_source) += recycle_energy_flow / volume;
	  pnt.ythis(energy_source) += recycle_energy_flow / volume;
	}
	  
      }); // end yboundary.iter_pnts()

    }

    // Put the updated sources back into the state
    set<Field3D>(species_to["density_source"], density_source);
    set<Field3D>(species_to["energy_source"], energy_source);

  }
}

void RecyclingFCI::outputVars(Options& state) {
  AUTO_TRACE();

  if (neutral_pump) {
    // Save the pump mask as a time-independent field
    set_with_attrs(state[{std::string("is_pump")}], is_pump,
                   {{"standard_name", "neutral pump location"},
                    {"long_name", std::string("Neutral pump location")},
                    {"source", "recycling"}});
  }

  if (diagnose) {
    // Normalisations
    auto Nnorm = get<BoutReal>(state["Nnorm"]);
    auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
    auto Tnorm = get<BoutReal>(state["Tnorm"]);
    BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

    for (const auto& channel : channels) {
      // Save particle and energy source for the species created during recycling

      // Target recycling
      if (target_recycle) {
        set_with_attrs(state[{std::string("S") + channel.to + std::string("_target_recycle")}],
                       channel.target_recycle_density_source,
                       {{"time_dimension", "t"},
                        {"units", "m^-3 s^-1"},
                        {"conversion", Nnorm * Omega_ci},
                        {"standard_name", "particle source"},
                        {"long_name", std::string("Target recycling particle source of ") + channel.to},
                        {"source", "recycling"}});

        set_with_attrs(state[{std::string("E") + channel.to + std::string("_target_recycle")}],
                       channel.target_recycle_energy_source,
                       {{"time_dimension", "t"},
                        {"units", "W m^-3"},
                        {"conversion", Pnorm * Omega_ci},
                        {"standard_name", "energy source"},
                        {"long_name", std::string("Target recycling energy source of ") + channel.to},
                        {"source", "recycling"}});
      }

      // Wall recycling
      if ((sol_recycle) or (pfr_recycle)) {
        set_with_attrs(state[{std::string("S") + channel.to + std::string("_wall_recycle")}],
                       channel.wall_recycle_density_source,
                       {{"time_dimension", "t"},
                        {"units", "m^-3 s^-1"},
                        {"conversion", Nnorm * Omega_ci},
                        {"standard_name", "particle source"},
                        {"long_name", std::string("Wall recycling particle source of ") + channel.to},
                        {"source", "recycling"}});

        set_with_attrs(state[{std::string("E") + channel.to + std::string("_wall_recycle")}],
                       channel.wall_recycle_energy_source,
                       {{"time_dimension", "t"},
                        {"units", "W m^-3"},
                        {"conversion", Pnorm * Omega_ci},
                        {"standard_name", "energy source"},
                        {"long_name", std::string("Wall recycling energy source of ") + channel.to},
                        {"source", "recycling"}});
      }

      // Neutral pump
      if (neutral_pump) {
        set_with_attrs(state[{std::string("S") + channel.to + std::string("_pump")}],
                       channel.pump_density_source,
                       {{"time_dimension", "t"},
                        {"units", "m^-3 s^-1"},
                        {"conversion", Nnorm * Omega_ci},
                        {"standard_name", "particle source"},
                        {"long_name", std::string("Pump recycling particle source of ") + channel.to},
                        {"source", "recycling"}});

        set_with_attrs(state[{std::string("E") + channel.to + std::string("_pump")}],
                       channel.pump_energy_source,
                       {{"time_dimension", "t"},
                        {"units", "W m^-3"},
                        {"conversion", Pnorm * Omega_ci},
                        {"standard_name", "energy source"},
                        {"long_name", std::string("Pump recycling energy source of ") + channel.to},
                        {"source", "recycling"}});
      }

    }
  }
}
