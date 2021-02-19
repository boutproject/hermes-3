
#include <difops.hxx>

#include "../include/thermal_force.hxx"

void ThermalForce::transform(Options& state) {
  AUTO_TRACE();

  Options& allspecies = state["species"];

  if (electron_ion && allspecies.isSection("e")) {
    // Electron-ion collisions

    Options& electrons = allspecies["e"];
    const Field3D Te = get<Field3D>(electrons["temperature"]);
    const Field3D Grad_Te = Grad_par(Te);

    for (auto& kv : allspecies.getChildren()) {
      if (kv.first == "e") {
        continue; // Omit electron-electron
      }
      Options& species = allspecies[kv.first];

      if (!species.isSet("charge")) {
        continue; // Only considering charged particle interactions
      }

      const BoutReal Z = get<BoutReal>(species["charge"]);
      const Field3D nz = get<Field3D>(species["density"]);

      Field3D ion_force = nz * (0.71 * SQ(Z)) * Grad_Te;

      add(species["momentum_source"], ion_force);
      subtract(electrons["momentum_source"], ion_force);
    }
  }

  if (ion_ion) {
    // Iterate through other species
    // To avoid double counting, this needs to iterate over pairs
    // not including the same species. i.e. above the diagonal
    //
    // Iterators kv1 and kv2 over the species map
    //
    //               kv2 ->
    //             species1  species2  species3
    // kv1   species1               X         X
    //  ||   species2                         X
    //  \/   species3
    //
    const std::map<std::string, Options>& children = allspecies.getChildren();
    for (auto kv1 = std::begin(children); kv1 != std::end(children); ++kv1) {
      Options& species1 = allspecies[kv1->first];

      if (kv1->first == "e" or !species1.isSet("charge")) {
        continue; // Only considering charged particle interactions
      }

      // Copy the iterator, so we don't iterate over the
      // lower half of the matrix or the diagonal but start  the diagonal
      for (std::map<std::string, Options>::const_iterator kv2 = std::next(kv1);
           kv2 != std::end(children); ++kv2) {
        Options& species2 = allspecies[kv2->first];

        if (kv2->first == "e" or !species2.isSet("charge")) {
          continue; // Only considering charged particle interactions
        }

        // Now have two different ion species, species1 and species2
        // Only including one majority light species, and one trace heavy species

        Options *light, *heavy;
        if ((get<BoutReal>(species1["AA"]) < 4)
            and (get<BoutReal>(species2["AA"]) > 10)) {
          // species1 light, species2 heavy
          light = &species1;
          heavy = &species2;
        } else if (((get<BoutReal>(species1["AA"]) > 10)
                    and (get<BoutReal>(species2["AA"]) < 4))) {
          // species1 heavy, species2 light
          light = &species2;
          heavy = &species1;
        } else {
          // Ignore this combination
          if (first_time) {
            // Print warnings the first time
            output_warn.write(
                "ThermalForce: Not calculating thermal force between {} and {} species\n",
                kv1->first, kv2->first);
          }
          continue;
        }

        // Force on heavy (trace) ion due to light species
        // This follows Stangeby, page 298 and following

        const BoutReal mi = get<BoutReal>((*light)["AA"]);
        const Field3D Ti = get<Field3D>((*light)["temperature"]);

        const BoutReal mz = get<BoutReal>((*heavy)["AA"]);
        const BoutReal Z = get<BoutReal>((*heavy)["charge"]);
        const Field3D nz = get<Field3D>((*heavy)["density"]);

        if (Z == 0.0) {
          continue; // Check that the charge is not zero
        }

        const BoutReal mu = mz / (mi + mz);
        const BoutReal beta =
            3
            * (mu + 5 * sqrt(2) * SQ(Z) * (1.1 * pow(mu, 5. / 2) - 0.35 * pow(mu, 3. / 2))
               - 1)
            / (2.6 - 2 * mu + 5.4 * SQ(mu));

        const Field3D heavy_force = nz * beta * Grad_par(Ti);

        add((*heavy)["momentum_source"], heavy_force);
        subtract((*light)["momentum_source"], heavy_force);
      }
    }
  }
  first_time = false; // Don't print warnings every time
}
