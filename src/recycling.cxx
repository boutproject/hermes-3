
#include "../include/recycling.hxx"

#include <bout/utils.hxx> // for trim, strsplit
#include <bout/coordinates.hxx>
#include <bout/mesh.hxx>
#include <bout/constants.hxx>

using bout::globals::mesh;

Recycling::Recycling(std::string name, Options& alloptions, Solver*) {
  AUTO_TRACE();

  const Options& units = alloptions["units"];
  const BoutReal Tnorm = units["eV"];

  Options& options = alloptions[name];

  auto species_list = strsplit(options["species"]
                                   .doc("Comma-separated list of species to recycle")
                                   .as<std::string>(),
                               ',');

      
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
            .as<BoutReal>();

    BoutReal sol_recycle_multiplier =
        from_options["sol_recycle_multiplier"]
            .doc("Multiply the sol recycled flux by this factor. Should be >=0 and <= 1")
            .as<BoutReal>();

    BoutReal pfr_recycle_multiplier =
        from_options["pfr_recycle_multiplier"]
            .doc("Multiply the pfr recycled flux by this factor. Should be >=0 and <= 1")
            .as<BoutReal>();

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

    if ((target_recycle_multiplier < 0.0) or (target_recycle_multiplier > 1.0)
    or (sol_recycle_multiplier < 0.0) or (sol_recycle_multiplier > 1.0)
    or (pfr_recycle_multiplier < 0.0) or (pfr_recycle_multiplier > 1.0)) {
      throw BoutException("recycle_fraction must be betweeen 0 and 1");
    }

    // Populate recycling channel vector
    channels.push_back({
      from, to, 
      target_recycle_multiplier, sol_recycle_multiplier, pfr_recycle_multiplier,
      target_recycle_energy, sol_recycle_energy, pfr_recycle_energy});

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
  }
}

void Recycling::transform(Options& state) {
  AUTO_TRACE();

  // Get metric tensor components
  Coordinates* coord = mesh->getCoordinates();
  const Field2D& J = coord->J;
  const Field2D& dy = coord->dy;
  const Field2D& dx = coord->dx;
  const Field2D& dz = coord->dz;
  const Field2D& g_22 = coord->g_22;
  const Field2D& g11 = coord->g11;

  for (const auto& channel : channels) {
    const Options& species_from = state["species"][channel.from];

    const Field3D N = get<Field3D>(species_from["density"]);
    const Field3D V = get<Field3D>(species_from["velocity"]); // Parallel flow velocity

    Options& species_to = state["species"][channel.to];

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

      target_recycle_density_source = 0;
      target_recycle_energy_source = 0;

      // Lower Y boundary

      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          // Calculate flux through surface [normalised m^-2 s^-1],
          // should be positive since V < 0.0
          BoutReal flux =
              -0.5 * (N(r.ind, mesh->ystart, jz) + N(r.ind, mesh->ystart - 1, jz)) * 0.5
              * (V(r.ind, mesh->ystart, jz) + V(r.ind, mesh->ystart - 1, jz));

          if (flux < 0.0) {
            flux = 0.0;
          }

          // Flow of recycled species inwards
          BoutReal flow =
              channel.target_multiplier * flux
              * (J(r.ind, mesh->ystart) + J(r.ind, mesh->ystart - 1))
              / (sqrt(g_22(r.ind, mesh->ystart)) + sqrt(g_22(r.ind, mesh->ystart - 1)));

          // Add to density source
          target_recycle_density_source += flow / (J(r.ind, mesh->ystart) * dy(r.ind, mesh->ystart));
          density_source(r.ind, mesh->ystart, jz) += flow / (J(r.ind, mesh->ystart) * dy(r.ind, mesh->ystart));

          // energy of recycled particles
          target_recycle_energy_source += channel.target_energy * flow / (J(r.ind, mesh->ystart) * dy(r.ind, mesh->ystart));
          energy_source(r.ind, mesh->ystart, jz) += channel.target_energy * flow / (J(r.ind, mesh->ystart) * dy(r.ind, mesh->ystart));
        }
      }

      // Upper Y boundary

      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        // Calculate flux of ions into target from Ne and Vi boundary
        // This calculation is supposed to be consistent with the flow
        // of plasma from FV::Div_par(N, V)
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          // Flux through surface [normalised m^-2 s^-1], should be positive
          BoutReal flux = 0.5 * (N(r.ind, mesh->yend, jz) + N(r.ind, mesh->yend + 1, jz))
                          * 0.5 * (V(r.ind, mesh->yend, jz) + V(r.ind, mesh->yend + 1, jz));

          if (flux < 0.0) {
            flux = 0.0;
          }

          // Flow of neutrals inwards
          BoutReal flow =
              channel.target_multiplier * flux * (J(r.ind, mesh->yend) + J(r.ind, mesh->yend + 1))
              / (sqrt(g_22(r.ind, mesh->yend)) + sqrt(g_22(r.ind, mesh->yend + 1)));

          // Rate of change of neutrals in final cell
          // Add to density source
          target_recycle_density_source += flow / (J(r.ind, mesh->yend) * dy(r.ind, mesh->yend));
          density_source(r.ind, mesh->yend, jz) += flow / (J(r.ind, mesh->yend) * dy(r.ind, mesh->yend));

          target_recycle_energy_source += channel.target_energy * flow / (J(r.ind, mesh->yend) * dy(r.ind, mesh->yend));
          energy_source(r.ind, mesh->yend, jz) += channel.target_energy * flow / (J(r.ind, mesh->yend) * dy(r.ind, mesh->yend));
        }
      }
    }

    // Recycling at the SOL edge (2D/3D only)
    if (sol_recycle) {

      radial_particle_flow = get<Field3D>(species_from["particle_flow_xlow"]);
      radial_energy_flow = get<Field3D>(species_from["energy_flow_xlow"]);
      sol_recycle_density_source = 0;
      sol_recycle_energy_source = 0;

      if(mesh->lastX()){
        for(int iy=0; iy < mesh->LocalNy ; iy++){
          for(int iz=0; iz < mesh->LocalNz; iz++){

            // Volume of cell adjacent to wall which will receive source
            BoutReal volume = J(mesh->xend, iy) * dx(mesh->xend, iy)
                 * dy(mesh->xend, iy) * dz(mesh->xend, iy);

            // Flow of recycled species back from the edge
            // Edge = LHS flow of inner guard cells (mesh->xend-1)
            // TODO: Handle cases when flow is going into domain from edge
            BoutReal recycle_particle_flow = channel.sol_multiplier * radial_particle_flow(mesh->xend+1, iy, iz) * -1; 
            BoutReal recycle_energy_flow = channel.sol_multiplier * radial_energy_flow(mesh->xend, iy, iz) * -1 ;

            // Divide by volume to get source
            sol_recycle_density_source(mesh->xend, iy, iz) += recycle_particle_flow / volume;
            density_source(mesh->xend, iy, iz) += sol_recycle_density_source(mesh->xend, iy, iz);

            // For now, this is a fixed temperature
            sol_recycle_energy_source(mesh->xend, iy, iz) += channel.sol_energy * recycle_particle_flow / volume;
            energy_source(mesh->xend, iy, iz) += sol_recycle_energy_source(mesh->xend, iy, iz);

          

          }
        }
      }
    }

        // Recycling at the SOL edge (2D/3D only)
    if (pfr_recycle) {

      throw BoutException("Error: PFR recycling not implemented yet\n");
    }

    // Put the updated sources back into the state
    set<Field3D>(species_to["density_source"], density_source);
    set<Field3D>(species_to["energy_source"], energy_source);
  }
}

void Recycling::outputVars(Options& state) {

  AUTO_TRACE();
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

  if (diagnose) {

      for (const auto& channel : channels) {
        AUTO_TRACE();

        // Save particle and energy source for the species created during recycling

        if (target_recycle) {
          set_with_attrs(state[{std::string("S") + channel.to + std::string("_target_recycle")}], target_recycle_density_source,
                          {{"time_dimension", "t"},
                          {"units", "m^-3 s^-1"},
                          {"conversion", Nnorm * Omega_ci},
                          {"standard_name", "particle source"},
                          {"long_name", std::string("Target recycling particle source of ") + channel.to},
                          {"source", "recycling"}});
    
          set_with_attrs(state[{std::string("E") + channel.to + std::string("_target_recycle")}], target_recycle_energy_source,
                          {{"time_dimension", "t"},
                          {"units", "W m^-3"},
                          {"conversion", Pnorm * Omega_ci},
                          {"standard_name", "energy source"},
                          {"long_name", std::string("Target recycling energy source of ") + channel.to},
                          {"source", "recycling"}});
          }

        if (sol_recycle) {
          set_with_attrs(state[{std::string("S") + channel.to + std::string("_sol_recycle")}], sol_recycle_density_source,
                          {{"time_dimension", "t"},
                          {"units", "m^-3 s^-1"},
                          {"conversion", Nnorm * Omega_ci},
                          {"standard_name", "particle source"},
                          {"long_name", std::string("SOL recycling particle source of ") + channel.to},
                          {"source", "recycling"}});
    
          set_with_attrs(state[{std::string("E") + channel.to + std::string("_sol_recycle")}], sol_recycle_energy_source,
                          {{"time_dimension", "t"},
                          {"units", "W m^-3"},
                          {"conversion", Pnorm * Omega_ci},
                          {"standard_name", "energy source"},
                          {"long_name", std::string("SOL recycling energy source of ") + channel.to},
                          {"source", "recycling"}});
          }
      }

  }
}
