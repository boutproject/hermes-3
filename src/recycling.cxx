
#include "../include/recycling.hxx"

#include <bout/utils.hxx> // for trim, strsplit
#include <bout/coordinates.hxx>
#include <bout/mesh.hxx>

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

    BoutReal recycle_multiplier =
        from_options["recycle_multiplier"]
            .doc("Multiply the recycled flux by this factor. Should be >=0 and <= 1")
            .as<BoutReal>();

    BoutReal recycle_energy = from_options["recycle_energy"]
                                  .doc("Fixed energy of the recycled particles [eV]")
                                  .withDefault<BoutReal>(0.0)
                              / Tnorm; // Normalise from eV

    if ((recycle_multiplier < 0.0) or (recycle_multiplier > 1.0)) {
      throw BoutException("recycle_fraction must be betweeen 0 and 1");
    }
    channels.push_back({from, to, recycle_multiplier, recycle_energy});
  }
}

void Recycling::transform(Options& state) {
  AUTO_TRACE();

  // Get metric tensor components
  Coordinates* coord = mesh->getCoordinates();
  const Field3D& J = coord->J;
  const Field3D& dy = coord->dy;
  const Field3D& g_22 = coord->g_22;

  for (const auto& channel : channels) {
    const Options& species_from = state["species"][channel.from];

    const Field3D N = get<Field3D>(species_from["density"]);
    const Field3D V = get<Field3D>(species_from["velocity"]); // Parallel flow velocity

    Options& species_to = state["species"][channel.to];
    // Get the sources, so the values can be added
    Field3D density_source = species_to.isSet("density_source")
                                 ? getNonFinal<Field3D>(species_to["density_source"])
                                 : 0.0;
    Field3D energy_source = species_to.isSet("energy_source")
                                ? getNonFinal<Field3D>(species_to["energy_source"])
                                : 0.0;

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
            channel.multiplier * flux
	  * (J(r.ind, mesh->ystart, jz) + J(r.ind, mesh->ystart - 1, jz))
	  / (sqrt(g_22(r.ind, mesh->ystart, jz)) + sqrt(g_22(r.ind, mesh->ystart - 1, jz)));

        // Add to density source
        density_source(r.ind, mesh->ystart, jz) +=
	  flow / (J(r.ind, mesh->ystart, jz) * dy(r.ind, mesh->ystart, jz));

        // energy of recycled particles
        energy_source(r.ind, mesh->ystart, jz) +=
	  channel.energy * flow / (J(r.ind, mesh->ystart, jz) * dy(r.ind, mesh->ystart, jz));
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
	  channel.multiplier * flux * (J(r.ind, mesh->yend, jz) + J(r.ind, mesh->yend + 1, jz))
	  / (sqrt(g_22(r.ind, mesh->yend, jz)) + sqrt(g_22(r.ind, mesh->yend + 1, jz)));

        // Rate of change of neutrals in final cell
        // Add to density source
        density_source(r.ind, mesh->yend, jz) +=
	  flow / (J(r.ind, mesh->yend, jz) * dy(r.ind, mesh->yend, jz));

        energy_source(r.ind, mesh->yend, jz) +=
	  channel.energy * flow / (J(r.ind, mesh->yend, jz) * dy(r.ind, mesh->yend, jz));
      }
    }

    // Put the updated sources back into the state
    set<Field3D>(species_to["density_source"], density_source);
    set<Field3D>(species_to["energy_source"], energy_source);
  }
}
